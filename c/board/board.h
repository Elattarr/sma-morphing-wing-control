/* Board layer: glues the portable controller / safety monitor (../sma.c) and the 20 kHz current loop (fast_loop.c) to the hardware.
 * Everything that touches a register is behind the hal_* functions below: hal_sim.c implements them for the PC test, hal_stm32g4.c
 * is the (untested) list of what the NUCLEO-G474RE must do.  NOT validated on hardware. */
#ifndef BOARD_H
#define BOARD_H

#include <stdint.h>
#include "../sma.h"
#include "fast_loop.h"

/* ---- what the target must provide ---- */
uint16_t hal_adc_current_counts(void);              /* INA240 output, sampled by HRTIM trigger in the middle of the Q1 on-time */
uint16_t hal_adc_vbus_counts(void);
uint16_t hal_adc_laser_counts(void);                /* slow, any time */
int      hal_tc_read(int channel, float *T_C);      /* MAX31856 over SPI; returns 0 if ok, non-zero on a fault bit or SPI error */
void     hal_pwm_set_duty(float D);                 /* Q1 duty for the next period (HRTIM compare preload); Q2 complementary with dead time */
void     hal_pwm_enable(int on);                    /* 0: both gates low (HRTIM outputs forced inactive) */
void     hal_relay_set(int closed);                 /* PB0 -> Q3 -> K1 coil */
int      hal_hw_trip(void);                         /* latched hardware over-current / over-temperature comparator (if built), 1 = tripped */
void     hal_uart_tx(const char *s);                /* LPUART1 to the laptop */

/* ---- board layer API ---- */
typedef enum { B_RELAY_OPEN = 0, B_RELAY_CLOSING, B_RUN, B_OPENING } board_state_t;

typedef struct {
    sma_wire_t w; sma_meff_t M; sma_cfg_t cfg; sma_safety_cfg_t sc;
    sma_ctl_t ctl; sma_safety_t saf;
    fl_t fl;
    board_state_t state;
    int armed, cal_from_ctl;
    long beat;
    float r_target;
    float laser_x0;                                  /* laser reading when armed (wire cold, extended) */
    float d_cold_mm;                                 /* displacement the controller expects at that moment */
    /* last slow-period values (for the log) */
    float T1, T2, d, I_cmd, I_perm, I_mean, I_rms, V_wire, Vbus, R_est;
    int tc_fault;
    sma_trip_t trip;
    int hw_trip;
} board_t;

void board_init(board_t *b, const sma_wire_t *w, const sma_meff_t *M, const sma_cfg_t *cfg, const sma_safety_cfg_t *sc,
                float sigma_ext, float K_spring, float d_cold_mm);
void board_set_target(board_t *b, float r);         /* from the laptop: 0..1 */
void board_fast_isr(board_t *b);                    /* call at 20 kHz from the ADC end-of-conversion interrupt */
void board_slow_tick(board_t *b);                   /* call every 50 ms from the main loop (flag set by a timer) */
void board_rx_line(board_t *b, const char *line);   /* "R 0.55" sets the target; "X" opens the relay (stop) */
float board_laser_to_mm(const board_t *b);

#endif
