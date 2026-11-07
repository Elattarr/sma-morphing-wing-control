/* PC simulation of the power stage and the sensors, behind the hal_* interface of board.h.  Used by tests/test_board.py (ctypes).
 * Power stage: average-value model of the half-bridge, one PWM period per fast tick:  L dI/dt = D_eff Vbus - I (R_wire + R_other),
 * the duty written in tick k acts in period k+1 (HRTIM compare preload), the ADC reads the period mean (sample in the middle of the
 * on-time), the dead time costs 2 t_dead f of duty.  Not modelled: switching ripple, thermal drift of RDS(on), noise coupling. */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "board.h"
#include "board_cfg.h"
#include "board_params.h"

static board_t B;
static double I_L, I_avg_period, R_w = 0.06, R_oth_true = 0.012, Vbus_true = 5.0, L_coil = 10e-6;
static double duty_next, duty_applied;
static int pwm_on, relay_closed, hw_trip_flag, tc_fail1, tc_fail2, arcs;
static double T1_in = 22, T2_in = 22, d_in = 6.5415, ina_gain_err = 0, ina_off_V = 0, vbus_gain_err = 0, noise_counts = 0;
static double acc_I, acc_I2, acc_D; static long acc_n;
static unsigned rng = 12345;
static char uart[4096];

static double nrand(void)                       /* ~N(0,1), cheap */
{
    double s = 0; int i;
    for (i = 0; i < 6; i++) { rng = rng * 1664525u + 1013904223u; s += (rng >> 8) / 16777216.0; }
    return (s - 3.0) * 1.4142;
}
static uint16_t quant(double v_adc)
{
    double c = v_adc / ADC_VREF * ADC_FULL + noise_counts * nrand();
    return (uint16_t)(c < 0 ? 0 : (c > ADC_FULL ? ADC_FULL : floor(c + 0.5)));
}

uint16_t hal_adc_current_counts(void)
{
    double v = I_avg_period * SHUNT_OHM * INA_GAIN * (1.0 + ina_gain_err) + ina_off_V + INA_OFFSET_V;
    return quant(v);
}
uint16_t hal_adc_vbus_counts(void) { return quant(Vbus_true * (1.0 + vbus_gain_err) * VBUS_DIV); }
uint16_t hal_adc_laser_counts(void)
{
    double x = (d_in - B.d_cold_mm) / (LASER_SIGN * LASER_RANGE_MM) + 0.01;           /* x - x0 = (d - d_cold)/(SIGN RANGE); wire cold at x0 = 0.01 of the range: the 6.54 mm stroke leaves only 0.1 mm of margin in a 6.64 mm laser */
    return quant((4.0 + 16.0 * x) * 1e-3 * LASER_R_OHM);
}
int hal_tc_read(int ch, float *T) { *T = (float)(ch == 1 ? T1_in : T2_in); return ch == 1 ? tc_fail1 : tc_fail2; }
void hal_pwm_set_duty(float D) { duty_next = D; }
void hal_pwm_enable(int on) { pwm_on = on; }
void hal_relay_set(int closed)
{
    if (!closed && relay_closed && I_L > 1.0) arcs++;       /* contacts opened under current */
    relay_closed = closed;
}
int hal_hw_trip(void) { return hw_trip_flag; }
void hal_uart_tx(const char *s) { strncat(uart, s, sizeof uart - strlen(uart) - 1); }

/* ------------------------------------------------------------ exports */
void sim_init(void)
{
    I_L = I_avg_period = 0; duty_next = duty_applied = 0; pwm_on = relay_closed = hw_trip_flag = tc_fail1 = tc_fail2 = arcs = 0;
    acc_I = acc_I2 = acc_D = 0; acc_n = 0; uart[0] = 0;
    R_w = 0.06; R_oth_true = 0.012; Vbus_true = 5.0; ina_gain_err = ina_off_V = vbus_gain_err = noise_counts = 0; T1_in = T2_in = 22.0; d_in = BOARD_D_COLD_MM;
    board_init(&B, &BOARD_WIRE, &BOARD_MEFF, &BOARD_CFG, &BOARD_SAFETY, BOARD_SIGMA, BOARD_K_SPRING, BOARD_D_COLD_MM);
}
void sim_set_sensors(double T1, double T2, double d_mm) { T1_in = T1; T2_in = T2; d_in = d_mm; }
void sim_set_plant(double R_wire, double Vbus, double R_other_true) { R_w = R_wire; Vbus_true = Vbus; R_oth_true = R_other_true; }
void sim_set_errors(double ina_gain, double ina_off, double vbus_gain, double noise) { ina_gain_err = ina_gain; ina_off_V = ina_off; vbus_gain_err = vbus_gain; noise_counts = noise; }
void sim_set_faults(int tc1, int tc2, int hw) { tc_fail1 = tc1; tc_fail2 = tc2; hw_trip_flag = hw; }

static void fast_ticks(int n, double *trace)
{
    const double Ts = TS_FAST, deadloss = 2.0 * DEAD_TIME_S * F_PWM_HZ;
    int k;
    for (k = 0; k < n; k++) {
        double D = (pwm_on && relay_closed) ? duty_applied : 0.0, Deff = D > 0 ? (D - deadloss > 0 ? D - deadloss : 0) : 0;
        double R = R_w + R_oth_true, V = Deff * Vbus_true, Iss = V / R, a = exp(-R * Ts / L_coil), x = R * Ts / L_coil;
        double I0 = relay_closed ? I_L : 0.0;
        I_avg_period = Iss + (I0 - Iss) * (1.0 - a) / x;
        I_L = Iss + (I0 - Iss) * a;
        if (I_L < 0) I_L = 0;                                   /* the synchronous Q2 would let it go negative; the wire has no source to push it */
        acc_I += I_avg_period; acc_I2 += I_avg_period * I_avg_period; acc_D += D; acc_n++;
        if (trace) trace[k] = I_avg_period;
        board_fast_isr(&B);
        duty_applied = duty_next;
    }
}

void sim_slow(double r_target)                                  /* one 50 ms slow tick, then 1000 fast ticks */
{
    acc_I = acc_I2 = acc_D = 0; acc_n = 0;
    board_set_target(&B, (float)r_target);
    board_slow_tick(&B);
    fast_ticks(FAST_PER_SLOW, 0);
}
void sim_rx(const char *line) { board_rx_line(&B, line); }
/* direct fast-loop test: slow loop bypassed, command written straight to the current loop */
void sim_direct(double I_cmd, double R_ff, int enable, int n, double *trace)
{
    relay_closed = 1; pwm_on = enable; fl_set_command(&B.fl, (float)I_cmd, (float)R_ff, enable);
    fast_ticks(n, trace);
}
/* out: 0 I_true_mean, 1 I_true_rms, 2 I_meas_mean(board), 3 V_wire_est, 4 R_est, 5 I_cmd, 6 I_perm, 7 state, 8 trip, 9 relay, 10 arcs, 11 latched_over,
 *      12 D_mean, 13 T1 board, 14 d board, 15 Vbus board, 16 r_f, 17 ctl fault, 18 hw_trip, 19 tc_fault */
void sim_get(double *o)
{
    double n = acc_n > 0 ? (double)acc_n : 1.0;
    o[0] = acc_I / n; o[1] = sqrt(acc_I2 / n); o[2] = B.I_mean; o[3] = B.V_wire; o[4] = B.R_est; o[5] = B.I_cmd; o[6] = B.I_perm;
    o[7] = B.state; o[8] = B.trip; o[9] = relay_closed; o[10] = arcs; o[11] = B.fl.latched_over; o[12] = acc_D / n; o[13] = B.T1; o[14] = B.d;
    o[15] = B.Vbus; o[16] = B.ctl.r_f; o[17] = B.ctl.fault; o[18] = B.hw_trip; o[19] = B.tc_fault;
}
const char *sim_uart(void) { return uart; }
