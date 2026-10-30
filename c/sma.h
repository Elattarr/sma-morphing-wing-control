/* SMA wire controller, homing routine and independent safety monitor: portable C99, no malloc, no libm beyond sqrt/log/fabs.
 *
 * Direct port of src/sma_control/{spring,cycle_controller,safety}.py: same equations, same order of operations, so that replaying
 * logged inputs gives the Python outputs (tests/test_c_port.py compares them sample by sample).  Not validated on hardware.
 *
 * sma_real is double by default (bit-comparable with Python); compile with -DSMA_FLOAT for a single-precision build (MCUs with a
 * single-precision FPU).  Memory: dominated by the homing records, SMA_HEAT_MAX + SMA_COOL_MAX samples x 3 values; reduce with
 * SMA_REC_DECIM (keeps every n-th sample; the window temperatures then have n*Ts time resolution).
 */
#ifndef SMA_H
#define SMA_H

#include <stdint.h>

#ifdef SMA_FLOAT
typedef float sma_real;
#else
typedef double sma_real;
#endif

#ifndef SMA_HEAT_MAX
#define SMA_HEAT_MAX 12000   /* samples of the homing heating phase (600 s at 20 Hz) */
#endif
#ifndef SMA_COOL_MAX
#define SMA_COOL_MAX 170000  /* samples of the homing cooling phase (8500 s at 20 Hz) */
#endif
#ifndef SMA_REC_DECIM
#define SMA_REC_DECIM 1
#endif
#define SMA_RING_MAX 512     /* longest moving-window length in samples */

/* ------------------------------------------------------------------ constants of the wire (from the host / identification) */
typedef struct {            /* MEFF + elasticity (uniaxial reduction) */
    sma_real T0, R, beta, eps_L, h, m, E_A, E_M;
} sma_meff_t;

typedef struct {            /* geometry, thermal and electrical data of ONE wire */
    sma_real L;             /* active length [m] */
    sma_real A;             /* cross-section [m2] */
    sma_real S;             /* surface [m2] */
    sma_real mc;            /* m*c [J/K] */
    sma_real mL;            /* m*Lambda [J] */
    sma_real rho_eM, rho_eA;/* resistivity martensite / austenite [ohm m] */
    sma_real h_conv;        /* convection coefficient [W/m2/K] */
    sma_real Imax;          /* current limit [A] */
} sma_wire_t;

typedef struct {
    sma_real Ts, Kp, Ki, k_soft, dTdt_max, margin_up, margin_dn, T_safe_above_fs, Kr, Ki_r, Kr_dn, Ki_r_dn, tau_filt, Kp_hold, Ki_hold,
             drdt_up, dn_frac, Kaw, db_start, db_done, T_amb_est, T_limit_above_Af, T_abs_max, T_limit_expand_rate, learn_margin,
             hold_step, hold_drift, hold_below_As, sigma_unc;
    int overshoot_exit;
    /* heat-only (one-way actuator) mode, see cycle_controller.py: cfg_heat_only */
    int heat_only, creep_confirm;
    sma_real creep_rate, creep_detect, q_tau, approach_gain, approach_span, approach_min_rate, fast_detect, creep_pcap, settle_dr, settle_s,
             full_min, near_tol, drop_tol, reset_r, reset_still_s, reset_timeout;
} sma_cfg_t;

/* ------------------------------------------------------------------ bias-spring stress model */
typedef struct {
    sma_meff_t M;
    sma_real s0, K, eps_ref, c0, q_rec, q_fwd, a, kappa, eps_tr;
} sma_stress_t;

void     sma_stress_init(sma_stress_t *st, const sma_meff_t *M, sma_real sigma_ext, sma_real K);
sma_real sma_strain(const sma_stress_t *st, sma_real xi);
sma_real sma_sigma(const sma_stress_t *st, sma_real xi);
sma_real sma_xi_from_strain(const sma_stress_t *st, sma_real eps);
sma_real sma_T_rec(const sma_stress_t *st, sma_real xi);
sma_real sma_T_fwd(const sma_stress_t *st, sma_real xi);

/* ------------------------------------------------------------------ controller */
typedef enum { SMA_HOLD = 0, SMA_UP = 1, SMA_DOWN = 2 } sma_mode_t;

typedef struct {
    sma_real eps_ext, eps_con, As, Af, Tfs, Tfe, T_peak, t_heat, t_cool;
} sma_cal_t;

typedef enum { SMA_C_OK = 0, SMA_C_POSITION_DROPPED = 1, SMA_C_STALLED_SHORT = 2, SMA_C_RESET_TIMEOUT = 3 } sma_cfault_t;

typedef struct {
    sma_wire_t w;
    sma_cfg_t c;
    sma_stress_t st;
    sma_real As, Af, Tfs, Tfe, W, m_eff, T_limit, T_pre, T_safe, rate_dn, T_motion;
    sma_mode_t mode;
    int in_loop;
    int Tcmd_valid;
    sma_real Tcmd, I_a, I_b, rref, r_est, r_f, r_entry, rf_entry;
    int r_f_valid;
    sma_real T_limit_model, T_lim, T_hold, T_hold_max, T_peak, r_hold_ref;
    int cal_valid;
    sma_cal_t cal;
    /* heat-only state */
    sma_cfault_t fault;
    int d_first_set, resetting, saturated, r_q_valid, d_valid, n_detect, Tg_valid, T_start_valid, settle_valid;
    sma_real T_ref, r_goal, r_max, r_q, d_q, d_first, Tg, T_start_obs, settle_ref, t_settle, t_reset, t_rest;
} sma_ctl_t;

void     sma_ctl_init(sma_ctl_t *ct, const sma_wire_t *w, const sma_meff_t *M, const sma_cfg_t *cfg, sma_real sigma_ext, sma_real K_spring);
void     sma_ctl_apply_cal(sma_ctl_t *ct, const sma_cal_t *cal);
sma_real sma_ctl_r_from_disp(const sma_ctl_t *ct, sma_real d_mm);
sma_real sma_resistance(const sma_wire_t *w, sma_real r);
/* one control step: target r in [0,1], thermocouple [C], displacement [mm]; returns the commanded current [A] */
sma_real sma_ctl_step(sma_ctl_t *ct, sma_real r_target, sma_real T_meas, sma_real d_meas_mm);
/* heat-only: call when the wing is being returned externally; power stays off until r < reset_r for reset_still_s (fault 'reset_timeout' otherwise) */
void     sma_ctl_begin_reset(sma_ctl_t *ct);

/* ------------------------------------------------------------------ homing */
typedef enum { SMA_H_WAIT = 0, SMA_H_HEAT = 1, SMA_H_COOL = 2 } sma_hphase_t;
typedef enum { SMA_F_NONE = 0, SMA_F_NOT_COLD, SMA_F_NO_MOTION_ON_HEATING, SMA_F_NO_EXTENSION, SMA_F_PRELOAD_TOO_LOW, SMA_F_BUFFER_FULL } sma_hfault_t;

typedef struct { sma_real v[SMA_RING_MAX]; int n, cap, head; } sma_ring_t;

typedef struct {
    sma_ctl_t *ctl;
    sma_real ramp, t_wait, thr, confirm_s, min_stroke_frac, cool_still_s, cool_still_max_s, cool_T_margin, preload_margin_K,
             cold_no_ext_s, t_cool_max, t_wait_max, tau_d, tau_fast, T_cold_max, tau_rate, motion_s, cold;
    int Tf_valid;
    sma_real Tf, I_rate;
    sma_ring_t Tf_hist, wait_T, wait_d, hist, cool_df, cool_T;
    sma_real min_stroke_mm;
    sma_hphase_t phase;
    sma_real t;
    int done;
    sma_hfault_t fault;
    int cal_valid;
    sma_cal_t cal;
    int df_valid;
    sma_real df, dff, calm, still, d_e0, goal, d_c, t_cool0, T_peak;
    int n_heat, n_cool, k_dec;
    sma_real (*heat)[3];   /* records: t, T_meas, fast-filtered displacement */
    sma_real (*cool)[3];   /* records: t, T_meas, filtered displacement */
} sma_home_t;

/* heat and cool must point to arrays of SMA_HEAT_MAX / SMA_COOL_MAX rows */
void     sma_home_init(sma_home_t *h, sma_ctl_t *ctl, sma_real (*heat)[3], sma_real (*cool)[3]);
sma_real sma_home_step(sma_home_t *h, sma_real T_meas, sma_real d_meas_mm);

/* ------------------------------------------------------------------ independent safety monitor */
typedef enum { SMA_S_OK = 0, SMA_S_SENSOR_RANGE, SMA_S_OVER_TEMPERATURE, SMA_S_TC_DISAGREE, SMA_S_TC_READS_HIGH, SMA_S_TC_READS_LOW,
               SMA_S_STALL, SMA_S_CURRENT_MISMATCH, SMA_S_RESISTANCE, SMA_S_WATCHDOG, SMA_S_TC_FROZEN } sma_trip_t;

typedef struct {
    sma_real T_trip, T_disagree, T_disagree_dyn, h_min, h_max, dT_amb, band, T_lat_lo, stall_dT, t_stall, stall_frac, I_tol_abs, I_tol_rel,
             t_I_mismatch, R_tol, mc_min_frac, mc_max_frac, lat_max_frac, tc_lag_max, T_range_lo, T_range_hi, t_tc_freeze;
    int n_freeze;
} sma_safety_cfg_t;

typedef struct { sma_real hi, lo, lof, lat; } sma_obs_t;

typedef struct {
    sma_safety_cfg_t c;
    sma_obs_t obs[2];       /* envelope observers: small and large m c */
    sma_real mcs[2], mL_max;
    sma_real Ts, mc, mL, S, R_min, R_max, Imax;
    int tripped, armed;
    sma_trip_t reason;
    sma_real t, T_lo, T_lo_f, T_hi, Tamb0;
    sma_real I_prev, t_Imis, t_hot, d_hot0;
    int d_hot_valid, T_prev_valid, beat_valid, n_same;
    sma_real T1_prev, T2_prev, T1_frz, t_frz;
    int T1_frz_valid;
    long beat;
    int stroke_valid;
    sma_real stroke_mm, Af_max;
} sma_safety_t;

void sma_safety_init(sma_safety_t *s, const sma_wire_t *w, const sma_meff_t *M, sma_real Ts, const sma_safety_cfg_t *cfg);
void sma_safety_arm(sma_safety_t *s, sma_real T1, int have_T2, sma_real T2);
void sma_safety_set_cal(sma_safety_t *s, sma_real stroke_mm, sma_real Af_max);
/* optional inputs are passed with a have_* flag; returns the permitted current and writes the relay state (1 = closed) */
sma_real sma_safety_step(sma_safety_t *s, sma_real I_cmd, sma_real I_meas, sma_real T1, int have_T2, sma_real T2, int have_d, sma_real d,
                         int have_V, sma_real V, int have_beat, long beat, int *relay_closed);

#endif
