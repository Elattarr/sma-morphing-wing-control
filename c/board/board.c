#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "board.h"
#include "board_cfg.h"

float board_laser_to_mm(const board_t *b)
{
    float v = (float)hal_adc_laser_counts() * (ADC_VREF / ADC_FULL);
    float mA = v / LASER_R_OHM * 1000.0f;
    float x = (mA - 4.0f) / 16.0f;
    return b->d_cold_mm + LASER_SIGN * (x - b->laser_x0) * LASER_RANGE_MM;
}

static float laser_x(void)
{
    float v = (float)hal_adc_laser_counts() * (ADC_VREF / ADC_FULL);
    return (v / LASER_R_OHM * 1000.0f - 4.0f) / 16.0f;
}

void board_init(board_t *b, const sma_wire_t *w, const sma_meff_t *M, const sma_cfg_t *cfg, const sma_safety_cfg_t *sc,
                float sigma_ext, float K_spring, float d_cold_mm)
{
    b->w = *w; b->M = *M; b->cfg = *cfg; b->sc = *sc;
    sma_ctl_init(&b->ctl, &b->w, &b->M, &b->cfg, (sma_real)sigma_ext, (sma_real)K_spring);
    sma_safety_init(&b->saf, &b->w, &b->M, b->cfg.Ts, &b->sc);
    fl_init(&b->fl);
    b->state = B_RELAY_OPEN; b->armed = 0; b->cal_from_ctl = 0; b->beat = 0; b->r_target = 0;
    b->laser_x0 = 0; b->d_cold_mm = d_cold_mm;
    b->T1 = b->T2 = b->d = b->I_cmd = b->I_perm = b->I_mean = b->I_rms = b->V_wire = b->R_est = 0; b->Vbus = 5.0f;
    b->tc_fault = 0; b->trip = SMA_S_OK; b->hw_trip = 0;
    hal_pwm_enable(0); hal_pwm_set_duty(0); hal_relay_set(0);
}

void board_set_target(board_t *b, float r)
{
    b->r_target = r < 0 ? 0 : (r > 1 ? 1 : r);
}

void board_rx_line(board_t *b, const char *line)
{
    if (line[0] == 'R' && line[1] == ' ') board_set_target(b, strtof(line + 2, NULL));
    else if (line[0] == 'X') b->state = B_OPENING;                /* stop: gates low, relay open (latches until the board is reset) */
}

void board_fast_isr(board_t *b)
{
    float D = fl_tick(&b->fl, hal_adc_current_counts(), hal_adc_vbus_counts());
    if (b->fl.latched_over) hal_pwm_enable(0);
    hal_pwm_set_duty(D);
}

void board_slow_tick(board_t *b)
{
    float T1 = 0, T2 = 0;
    b->beat++;
    int f1 = hal_tc_read(1, &T1), f2 = hal_tc_read(2, &T2);
    b->tc_fault = f1 | (f2 << 1);
    b->hw_trip = hal_hw_trip();
    fl_take_slow(&b->fl, &b->I_mean, &b->I_rms, &b->V_wire, &b->Vbus);

    if (!b->armed) {                                           /* arm with the wire cold and the heater off */
        if (b->tc_fault == 0) {
            sma_safety_arm(&b->saf, (sma_real)T1, 1, (sma_real)T2);
            b->laser_x0 = laser_x();
            b->armed = 1;
        }
        return;
    }
    if (b->tc_fault) { T1 = -999.0f; T2 = -999.0f; }          /* a thermocouple fault is a sensor_range trip in the monitor */
    float d = board_laser_to_mm(b);
    b->T1 = T1; b->T2 = T2; b->d = d;

    sma_real Ic = sma_ctl_step(&b->ctl, (sma_real)b->r_target, (sma_real)T1, (sma_real)d);
    if (!b->saf.stroke_valid)                                  /* heat-only, no homing: the monitor's stall check starts from the model stroke ... */
        sma_safety_set_cal(&b->saf, fabs(sma_strain(&b->ctl.st, 0) - sma_strain(&b->ctl.st, 1)) * b->w.L * (sma_real)1e3,
                           b->ctl.T_limit_model - (sma_real)7.5);
    if (b->ctl.cal_valid && !b->cal_from_ctl) {                /* ... and uses the learned one after the first stroke */
        b->cal_from_ctl = 1;
        sma_safety_set_cal(&b->saf, fabs(b->ctl.cal.eps_ext - b->ctl.cal.eps_con) * b->w.L * (sma_real)1e3, b->ctl.cal.Af);
    }
    int relay_ok = 1;
    sma_real Ip = sma_safety_step(&b->saf, Ic, (sma_real)b->I_mean, (sma_real)T1, 1, (sma_real)T2, 1, (sma_real)d,
                                  1, (sma_real)b->V_wire, 1, b->beat, &relay_ok);
    b->I_cmd = (float)Ic; b->I_perm = (float)Ip; b->trip = b->saf.reason;
    b->R_est = b->I_mean > 0.3f ? b->V_wire / b->I_mean : 0.0f;
    if (b->hw_trip || b->fl.latched_over || !relay_ok) b->state = B_OPENING;

    float r = b->ctl.r_f_valid ? (float)b->ctl.r_f : 0.0f;
    float R_ff = (float)sma_resistance(&b->w, (sma_real)r) + R_OTHER_OHM;
    switch (b->state) {
    case B_RELAY_OPEN:                                         /* gates low, then close the relay */
        fl_set_command(&b->fl, 0, R_ff, 0); hal_pwm_enable(0);
        hal_relay_set(1); b->state = B_RELAY_CLOSING;
        break;
    case B_RELAY_CLOSING:                                      /* one slow tick (50 ms) for the contacts to settle, then enable */
        fl_set_command(&b->fl, 0, R_ff, 0);
        hal_pwm_enable(1); b->state = B_RUN;
        break;
    case B_RUN:
        fl_set_command(&b->fl, b->I_perm, R_ff, 1);
        break;
    case B_OPENING:                                            /* gates first, relay one tick later: the contacts never break current */
    default:
        fl_set_command(&b->fl, 0, R_ff, 0); hal_pwm_enable(0); hal_pwm_set_duty(0);
        if (b->I_mean < 0.3f) hal_relay_set(0);
        break;
    }
    char line[160];
    snprintf(line, sizeof line, "%ld %.2f %.2f %.3f %.3f %.3f %.3f %.3f %d %d\n", b->beat, T1, T2, d, b->I_cmd, b->I_perm, b->I_mean, b->V_wire,
             (int)b->state, (int)b->trip);
    hal_uart_tx(line);
}
