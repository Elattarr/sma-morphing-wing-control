/* 20 kHz current loop: measured wire current -> PWM duty of the low-side MOSFET Q1 (Q2 is its complement with dead time, set in HRTIM).
 * Runs in the ADC end-of-conversion interrupt, sample taken in the middle of the Q1 on-time (HRTIM trigger) = the average current.
 * Pure arithmetic, no hardware access: the same file runs on the MCU and in the PC simulation. */
#ifndef FAST_LOOP_H
#define FAST_LOOP_H

#include <stdint.h>

typedef struct {
    float I_ref;          /* A, slew-limited command the loop follows */
    float I_cmd;          /* A, from the slow loop (already limited by the safety monitor) */
    float R_ff;           /* ohm, wire + other, for the feed-forward (slow loop updates it) */
    float integ;          /* duty */
    float D;              /* duty written last */
    int   enabled;        /* 0: gates low */
    int   over_ticks;
    int   latched_over;   /* software overcurrent latch */
    /* accumulators over one slow period */
    float s_I, s_V, s_I2, s_Iraw;
    int   n;
    float I_meas, Vbus;   /* last samples */
} fl_t;

void  fl_init(fl_t *f);
void  fl_set_command(fl_t *f, float I_cmd, float R_ff, int enable);   /* called by the slow loop */
/* one fast tick: ADC counts of the current and bus-voltage channels -> duty for the NEXT period */
float fl_tick(fl_t *f, uint16_t i_counts, uint16_t vbus_counts);
/* read and clear the slow-period accumulators: mean current, rms current, estimated wire voltage (see board_cfg.h) */
void  fl_take_slow(fl_t *f, float *I_mean, float *I_rms, float *V_wire, float *Vbus_mean);
float fl_counts_to_amps(uint16_t counts);
float fl_counts_to_vbus(uint16_t counts);

#endif
