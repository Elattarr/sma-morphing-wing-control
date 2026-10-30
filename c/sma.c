/* Port of spring.py, cycle_controller.py (controller + homing) and safety.py. See sma.h. */
#include "sma.h"
#include <math.h>

#ifdef SMA_FLOAT
#define SQRT sqrtf
#define FABS fabsf
#define LOG logf
#else
#define SQRT sqrt
#define FABS fabs
#define LOG log
#endif

static sma_real maxr(sma_real a, sma_real b) { return a > b ? a : b; }
static sma_real minr(sma_real a, sma_real b) { return a < b ? a : b; }
static sma_real clipr(sma_real x, sma_real lo, sma_real hi) { return minr(maxr(x, lo), hi); }   /* np.clip order */

/* ================================================================== stress model (spring.py) */
void sma_stress_init(sma_stress_t *st, const sma_meff_t *M, sma_real sigma_ext, sma_real K)
{
    st->M = *M;
    st->s0 = sigma_ext;
    st->K = K;
    st->kappa = SQRT((sma_real)2.0 / (sma_real)27.0);
    st->eps_tr = SQRT((sma_real)2.0 / (sma_real)3.0) * M->eps_L;
    st->eps_ref = sigma_ext / M->E_M + st->eps_tr;
    st->c0 = sigma_ext + K * st->eps_ref;
    st->q_rec = M->R / (1 - st->kappa * M->m);
    st->q_fwd = M->R / (1 + st->kappa * M->m);
    st->a = SQRT((sma_real)2.0 / (sma_real)3.0);
}

static sma_real E_of(const sma_stress_t *st, sma_real xi) { return st->M.E_A + xi * (st->M.E_M - st->M.E_A); }

sma_real sma_strain(const sma_stress_t *st, sma_real xi)
{
    sma_real E = E_of(st, xi);
    return (st->c0 / E + st->eps_tr * xi) / (1 + st->K / E);
}

sma_real sma_sigma(const sma_stress_t *st, sma_real xi) { return st->s0 + st->K * (st->eps_ref - sma_strain(st, xi)); }

sma_real sma_xi_from_strain(const sma_stress_t *st, sma_real eps)
{
    sma_real xi = (sma_real)0.5;
    for (int i = 0; i < 8; i++)
        xi = minr(maxr((eps - sma_sigma(st, xi) / E_of(st, xi)) / st->eps_tr, 0), 1);
    return xi;
}

sma_real sma_T_rec(const sma_stress_t *st, sma_real xi)
{
    const sma_meff_t *M = &st->M;
    return M->T0 + (st->a * sma_sigma(st, xi) + st->q_rec - M->h * M->eps_L * xi) / M->beta;
}

sma_real sma_T_fwd(const sma_stress_t *st, sma_real xi)
{
    const sma_meff_t *M = &st->M;
    return M->T0 + (st->a * sma_sigma(st, xi) - st->q_fwd - M->h * M->eps_L * xi) / M->beta;
}

sma_real sma_resistance(const sma_wire_t *w, sma_real r)
{
    sma_real rho = w->rho_eM + r * (w->rho_eA - w->rho_eM);
    return rho * w->L / w->A;
}

/* ================================================================== controller (cycle_controller.py) */
void sma_ctl_init(sma_ctl_t *ct, const sma_wire_t *w, const sma_meff_t *M, const sma_cfg_t *cfg, sma_real sigma_ext, sma_real K_spring)
{
    ct->w = *w;
    ct->c = *cfg;
    sma_stress_init(&ct->st, M, sigma_ext, K_spring);
    ct->As = sma_T_rec(&ct->st, 1);
    ct->Af = sma_T_rec(&ct->st, 0);
    ct->Tfs = sma_T_fwd(&ct->st, 0);
    ct->Tfe = sma_T_fwd(&ct->st, 1);
    ct->W = ct->Af - ct->As;
    ct->m_eff = w->mL + w->mc * ct->W;
    sma_real a = SQRT((sma_real)2.0 / (sma_real)3.0) / M->beta;
    ct->T_limit = ct->Af + cfg->T_limit_above_Af + a * cfg->sigma_unc;
    ct->T_pre = ct->As - cfg->margin_up;
    ct->T_safe = ct->Tfs + cfg->T_safe_above_fs;
    sma_real mid = (sma_real)0.5 * (ct->Tfs + ct->Tfe);
    ct->rate_dn = cfg->dn_frac * w->h_conv * w->S * (mid - cfg->T_amb_est) / ct->m_eff;
    if (cfg->heat_only) ct->rate_dn = 0;
    ct->mode = SMA_HOLD;
    ct->in_loop = 0;
    ct->Tcmd_valid = 0;
    ct->Tcmd = 0;
    ct->I_a = ct->I_b = 0;
    ct->rref = ct->r_est = ct->r_entry = ct->rf_entry = 0;
    ct->r_f = 0;
    ct->r_f_valid = 0;
    ct->T_limit_model = ct->T_limit;
    ct->T_lim = ct->T_limit;
    ct->T_hold = ct->T_safe;
    ct->T_hold_max = maxr(ct->T_safe, ct->As - cfg->hold_below_As);
    ct->T_peak = (sma_real)-1e9;
    ct->r_hold_ref = 0;
    ct->cal_valid = 0;
    ct->T_motion = ct->Tfs + a * cfg->sigma_unc + (sma_real)2.0;
    ct->fault = SMA_C_OK;
    ct->d_first_set = ct->resetting = ct->saturated = ct->r_q_valid = ct->d_valid = ct->n_detect = ct->Tg_valid = ct->T_start_valid = ct->settle_valid = 0;
    ct->r_goal = -1; ct->T_ref = 0;
    ct->r_max = ct->r_q = ct->d_q = ct->d_first = ct->Tg = ct->T_start_obs = ct->settle_ref = ct->t_settle = ct->t_reset = ct->t_rest = 0;
}

static void ctl_enter(sma_ctl_t *ct, sma_mode_t mode)
{
    ct->r_entry = ct->r_est;
    ct->r_f = ct->r_est;
    ct->r_f_valid = 1;
    ct->rf_entry = ct->r_f;
    ct->r_hold_ref = ct->r_f;
    ct->T_peak = (sma_real)-1e9;
    ct->mode = mode;
    ct->in_loop = 0;
    ct->Tcmd_valid = 0;
    ct->Tg_valid = 0;
    ct->I_a = ct->I_b = 0;
    ct->n_detect = 0;
    ct->settle_valid = 0;
    ct->t_settle = 0;
}

void sma_ctl_apply_cal(sma_ctl_t *ct, const sma_cal_t *cal)
{
    sma_cfg_t *c = &ct->c;
    ct->cal = *cal;
    ct->cal_valid = 1;
    ct->As = cal->As; ct->Af = cal->Af; ct->Tfs = cal->Tfs; ct->Tfe = cal->Tfe;
    ct->W = maxr(ct->Af - ct->As, (sma_real)0.5);
    ct->m_eff = ct->w.mL + ct->w.mc * ct->W;
    sma_real mid = (sma_real)0.5 * (ct->Tfs + ct->Tfe);
    ct->rate_dn = c->dn_frac * ct->w.h_conv * ct->w.S * maxr(mid - c->T_amb_est, 1) / ct->m_eff;
    if (c->heat_only) ct->rate_dn = 0;
    ct->T_pre = ct->As - c->margin_up;
    ct->T_safe = ct->Tfs + c->T_safe_above_fs;
    ct->T_limit = minr(ct->Af + c->T_limit_above_Af, c->T_abs_max);
    ct->T_limit_model = ct->T_lim = ct->T_limit;
    ct->T_hold = ct->T_safe;
    ct->T_hold_max = maxr(ct->T_safe, ct->As - c->hold_below_As);
    ct->T_motion = ct->Tfs + (sma_real)2.0;
    ctl_enter(ct, SMA_HOLD);
}

sma_real sma_ctl_r_from_disp(const sma_ctl_t *ct, sma_real d_mm)
{
    sma_real eps = d_mm / (ct->w.L * (sma_real)1e3);
    if (ct->cal_valid) {
        sma_real m1 = sma_strain(&ct->st, 1), m0 = sma_strain(&ct->st, 0);
        eps = m1 + (eps - ct->cal.eps_ext) * (m0 - m1) / (ct->cal.eps_con - ct->cal.eps_ext);
    } else if (ct->c.heat_only && ct->d_first_set) {
        eps = eps - ct->d_first / (ct->w.L * (sma_real)1e3) + sma_strain(&ct->st, 1);       /* one-point zero at power-up */
    }
    return 1 - sma_xi_from_strain(&ct->st, eps);
}

static sma_real temp_loop(sma_ctl_t *ct, sma_real goal, sma_real T_meas, sma_real Q, sma_real P_max, sma_real Kp, sma_real Ki, sma_real rate_cap)
{
    const sma_cfg_t *c = &ct->c;
    sma_real dmax = (rate_cap < 0) ? c->dTdt_max : minr(rate_cap, c->dTdt_max);      /* rate_cap < 0: none */
    if (!ct->Tcmd_valid) { ct->Tcmd = T_meas; ct->Tcmd_valid = 1; }
    sma_real step = clipr(c->k_soft * (goal - ct->Tcmd), -dmax, dmax) * c->Ts;
    if (FABS(goal - ct->Tcmd) > dmax / c->k_soft)
        step = ((goal - ct->Tcmd) > 0 ? (sma_real)1 : (sma_real)-1) * dmax * c->Ts;
    ct->Tcmd += step;
    sma_real e = ct->Tcmd - T_meas;
    sma_real P_un = Q + Kp * e + ct->I_a;
    sma_real P = clipr(P_un, 0, P_max);
    ct->I_a += c->Ts * (Ki * e + c->Kaw * (P - P_un));
    return P;
}

static sma_real pos_loop(sma_ctl_t *ct, sma_real r_target, sma_real r_meas, sma_real Q, sma_real P_max, sma_real rate_up, sma_real Kr, sma_real Ki_r)
{
    const sma_cfg_t *c = &ct->c;
    if (!ct->in_loop) { ct->in_loop = 1; ct->rref = r_meas; ct->I_b = 0; }
    sma_real step = (rate_up != 0) ? clipr(r_target - ct->rref, -ct->rate_dn * c->Ts, rate_up * c->Ts) : (sma_real)0;
    ct->rref += step;
    sma_real e = ct->rref - r_meas;
    sma_real mL = ct->m_eff;
    sma_real P_un = Q + mL * (step / c->Ts + Kr * e + ct->I_b);
    sma_real P = clipr(P_un, 0, P_max);
    ct->I_b += c->Ts * (Ki_r * e + c->Kaw * (P - P_un) / mL);
    return P;
}

/* ---- heat-only helpers (ports of CycleController._stopped / begin_reset / _reset_step) */
static void ctl_stopped(sma_ctl_t *ct, sma_real r_target)
{
    const sma_cfg_t *c = &ct->c;
    ct->r_goal = r_target;
    if (r_target >= (sma_real)0.97 && ct->r_q >= c->full_min) {
        ct->saturated = 1;
        if (!ct->cal_valid && ct->d_valid) {
            sma_real L = ct->w.L * (sma_real)1e3;
            sma_cal_t cal;
            cal.eps_ext = ct->d_first / L;
            cal.eps_con = ct->d_q / L;
            cal.As = ct->T_start_valid ? ct->T_start_obs : ct->As;
            cal.Af = ct->T_ref;
            cal.Tfs = ct->Tfs; cal.Tfe = ct->Tfe; cal.T_peak = ct->T_peak; cal.t_heat = 0; cal.t_cool = 0;
            sma_ctl_apply_cal(ct, &cal);
            ct->r_est = ct->r_f = ct->r_q = ct->r_max = sma_ctl_r_from_disp(ct, ct->d_q);
        } else {
            ctl_enter(ct, SMA_HOLD);
        }
    } else if (r_target < (sma_real)0.97 && ct->r_q >= r_target - c->near_tol) {
        ctl_enter(ct, SMA_HOLD);
    } else {
        ct->fault = SMA_C_STALLED_SHORT;
        ctl_enter(ct, SMA_HOLD);
    }
}

void sma_ctl_begin_reset(sma_ctl_t *ct)
{
    if (ct->fault != SMA_C_OK && ct->fault != SMA_C_POSITION_DROPPED && ct->fault != SMA_C_RESET_TIMEOUT && ct->fault != SMA_C_STALLED_SHORT) return;
    ct->fault = SMA_C_OK; ct->resetting = 1; ct->t_reset = 0; ct->t_rest = 0;
}

static sma_real ctl_reset_step(sma_ctl_t *ct)
{
    const sma_cfg_t *c = &ct->c;
    ct->t_reset += c->Ts;
    ct->t_rest = (ct->r_f < c->reset_r) ? ct->t_rest + c->Ts : 0;
    if (ct->t_rest >= c->reset_still_s) {
        ct->resetting = 0; ct->saturated = 0; ct->r_max = ct->r_f; ct->r_goal = -1;
        ctl_enter(ct, SMA_HOLD);
    } else if (ct->t_reset > c->reset_timeout) {
        ct->fault = SMA_C_RESET_TIMEOUT; ct->resetting = 0;
    }
    return 0;
}

sma_real sma_ctl_step(sma_ctl_t *ct, sma_real r_target, sma_real T_meas, sma_real d_meas_mm)
{
    const sma_cfg_t *c = &ct->c;
    if (c->heat_only && !ct->d_first_set) { ct->d_first = d_meas_mm; ct->d_first_set = 1; }
    sma_real r = sma_ctl_r_from_disp(ct, d_meas_mm);
    ct->r_est = r;
    ct->r_f = ct->r_f_valid ? ct->r_f + c->Ts / c->tau_filt * (r - ct->r_f) : r;
    ct->r_f_valid = 1;
    sma_real R = sma_resistance(&ct->w, r);
    sma_real P_max = ct->w.Imax * ct->w.Imax * R;
    sma_real Q = ct->w.h_conv * ct->w.S * maxr(T_meas - c->T_amb_est, 0);
    ct->T_peak = maxr(ct->T_peak, T_meas);
    if (c->heat_only) {
        sma_real a = c->Ts / c->q_tau;
        ct->r_q = ct->r_q_valid ? ct->r_q + a * (r - ct->r_q) : r;
        ct->r_q_valid = 1;
        ct->d_q = ct->d_valid ? ct->d_q + a * (d_meas_mm - ct->d_q) : d_meas_mm;
        ct->d_valid = 1;
        if (ct->fault != SMA_C_OK) return 0;
        if (ct->resetting) return ctl_reset_step(ct);
        ct->r_max = maxr(ct->r_max, ct->r_f);
        if (ct->r_f < ct->r_max - c->drop_tol) { ct->fault = SMA_C_POSITION_DROPPED; return 0; }
    }
    if (ct->mode == SMA_HOLD) {
        if (r_target > ct->r_f + c->db_start && (!c->heat_only || r_target > ct->r_goal + (sma_real)1e-9)) ctl_enter(ct, SMA_UP);
        else if (r_target < ct->r_f - c->db_start && !c->heat_only) ctl_enter(ct, SMA_DOWN);
    }
    sma_real P = 0;
    sma_real tol = (r_target > (sma_real)0.03 && r_target < (sma_real)0.97) ? c->db_done : (sma_real)0.03;
    if (ct->mode == SMA_UP) {
        int moved, pre;
        if (c->heat_only && !ct->in_loop) {
            ct->n_detect = (ct->r_q > ct->r_entry + c->creep_detect) ? ct->n_detect + 1 : 0;
            moved = ct->n_detect >= c->creep_confirm || r > ct->r_entry + c->fast_detect * c->creep_detect;
            if (moved && !ct->T_start_valid) {
                ct->T_start_obs = T_meas; ct->T_start_valid = 1;
                ct->T_pre = minr(ct->T_pre, T_meas - c->margin_up);
            }
        } else {
            moved = r >= ct->r_entry + c->db_start;
        }
        pre = !ct->in_loop && !moved;
        if (c->heat_only && pre && (ct->Tg_valid || T_meas >= ct->T_pre - (sma_real)0.1)) {
            if (!ct->Tg_valid) { ct->Tg = T_meas; ct->Tg_valid = 1; ct->Tcmd = T_meas; ct->Tcmd_valid = 1; ct->I_a = 0; }
            ct->Tg += c->creep_rate * c->Ts;
            if (T_meas > ct->T_lim - 1) ct->T_lim = minr(ct->T_lim + c->T_limit_expand_rate * c->Ts, c->T_abs_max);
            P = minr(temp_loop(ct, ct->Tg, T_meas, Q, P_max, c->Kp, c->Ki, -1), Q + c->creep_pcap * ct->w.mc * c->creep_rate);
            if (T_meas > ct->T_lim) P = 0;
        } else if (pre && T_meas < ct->T_pre - (sma_real)0.1) {
            sma_real cap = c->heat_only ? maxr(c->approach_min_rate, c->approach_gain * (ct->T_pre + c->approach_span - T_meas)) : (sma_real)-1;
            P = temp_loop(ct, ct->T_pre, T_meas, Q, P_max, c->Kp, c->Ki, cap);
        } else {
            P = pos_loop(ct, r_target, r, Q, P_max, c->drdt_up, c->Kr, c->Ki_r);
            if (ct->in_loop && T_meas > ct->T_lim - 1 && r_target - ct->r_f > c->db_start)
                ct->T_lim = minr(ct->T_lim + c->T_limit_expand_rate * c->Ts, c->T_abs_max);
            if (T_meas > ct->T_lim) P = 0;
            if (c->heat_only && ct->in_loop && (P > 0 || T_meas > ct->T_pre)) {
                if (!ct->settle_valid || FABS(ct->r_q - ct->settle_ref) > c->settle_dr) {
                    ct->settle_ref = ct->r_q; ct->settle_valid = 1; ct->t_settle = 0; ct->T_ref = T_meas;
                } else {
                    ct->t_settle += c->Ts;
                }
                if (ct->t_settle >= c->settle_s && r_target - ct->r_f > c->db_done) {
                    ctl_stopped(ct, r_target);
                    return 0;
                }
            }
            if (FABS(ct->rref - r_target) < (sma_real)1e-9 && FABS(r_target - r) < tol) {
                ct->T_lim = maxr(ct->T_limit_model, minr(ct->T_peak + c->learn_margin, c->T_abs_max));
                ct->T_hold = minr(maxr(ct->T_hold, ct->T_safe + maxr(ct->T_peak - ct->Af, 0)), ct->T_hold_max);
                ct->r_goal = r_target;
                ctl_enter(ct, SMA_HOLD);
            } else if (c->overshoot_exit && ct->r_f > r_target + c->db_start) {
                ctl_enter(ct, SMA_HOLD);
            }
        }
    } else if (ct->mode == SMA_DOWN) {
        if (!ct->in_loop && (T_meas > ct->T_motion || (T_meas > ct->Tfs + c->margin_dn && ct->r_f > ct->rf_entry - (sma_real)0.5 * c->db_start))) {
            P = 0;
        } else {
            P = pos_loop(ct, r_target, ct->r_f, Q, P_max, c->drdt_up, c->Kr_dn, c->Ki_r_dn);
            if (FABS(ct->rref - r_target) < (sma_real)1e-9 && FABS(ct->r_f - r_target) < tol) {
                ctl_enter(ct, SMA_HOLD);
            } else if (c->overshoot_exit && ct->in_loop && ct->r_f < r_target - c->db_start) {
                ctl_enter(ct, SMA_HOLD);
            }
        }
    } else {
        sma_real d = ct->r_f - ct->r_hold_ref;
        if (d < -c->hold_drift) {
            ct->T_hold = minr(ct->T_hold + c->hold_step, ct->T_hold_max);
            ct->r_hold_ref = ct->r_f;
        } else if (d > c->hold_drift) {
            ct->T_hold = maxr(ct->T_hold - c->hold_step, ct->T_safe);
            ct->r_hold_ref = ct->r_f;
        }
        P = c->heat_only ? (sma_real)0 : temp_loop(ct, ct->T_hold, T_meas, Q, P_max, c->Kp_hold, c->Ki_hold, -1);
    }
    if (T_meas > ct->T_lim || T_meas > c->T_abs_max) P = 0;
    return SQRT(P / R);
}

/* ================================================================== homing */
static void ring_init(sma_ring_t *r, int cap) { r->n = 0; r->head = 0; r->cap = cap < SMA_RING_MAX ? cap : SMA_RING_MAX; }
static void ring_push(sma_ring_t *r, sma_real x)
{
    if (r->n < r->cap) { r->v[(r->head + r->n) % r->cap] = x; r->n++; }
    else { r->v[r->head] = x; r->head = (r->head + 1) % r->cap; }
}
static sma_real ring_get(const sma_ring_t *r, int i) { return r->v[(r->head + i) % r->cap]; }
static int ring_full(const sma_ring_t *r) { return r->n == r->cap; }

void sma_home_init(sma_home_t *h, sma_ctl_t *ctl, sma_real (*heat)[3], sma_real (*cool)[3])
{
    const sma_cfg_t *c = &ctl->c;
    h->ctl = ctl;
    h->ramp = (sma_real)1.5; h->t_wait = (sma_real)5.0; h->motion_s = (sma_real)2.0; h->thr = (sma_real)0.05; h->confirm_s = (sma_real)3.0;
    h->min_stroke_frac = (sma_real)0.3; h->cool_still_s = (sma_real)20.0; h->cool_still_max_s = (sma_real)90.0;
    h->cool_T_margin = (sma_real)4.0; h->preload_margin_K = (sma_real)2.0; h->cold_no_ext_s = (sma_real)60.0; h->cold = 0;
    h->t_cool_max = (sma_real)8000.0; h->t_wait_max = (sma_real)600.0;
    h->tau_d = (sma_real)0.5; h->tau_fast = (sma_real)0.15; h->T_cold_max = (sma_real)45.0; h->tau_rate = (sma_real)10.0;
    h->Tf_valid = 0; h->Tf = 0; h->I_rate = 0;
    ring_init(&h->Tf_hist, (int)(2.0 / c->Ts) + 1);
    ring_init(&h->wait_T, (int)(h->t_wait / c->Ts));
    ring_init(&h->wait_d, (int)(h->t_wait / c->Ts));
    ring_init(&h->hist, (int)(h->motion_s / c->Ts) + 1);
    ring_init(&h->cool_df, (int)(h->cool_still_s / c->Ts) + 1);
    ring_init(&h->cool_T, (int)(h->cool_still_s / c->Ts) + 1);
    h->min_stroke_mm = h->min_stroke_frac * ctl->st.eps_tr * ctl->w.L * (sma_real)1e3;
    h->phase = SMA_H_WAIT; h->t = 0; h->done = 0; h->fault = SMA_F_NONE; h->cal_valid = 0;
    h->df_valid = 0; h->df = h->dff = 0; h->calm = h->still = 0; h->d_e0 = h->goal = h->d_c = h->t_cool0 = 0;
    h->T_peak = (sma_real)-1e9;
    h->n_heat = h->n_cool = 0; h->k_dec = 0;
    h->heat = heat; h->cool = cool;
}

static int home_moving(sma_home_t *h, sma_real d)
{
    sma_real Ts = h->ctl->c.Ts;
    sma_real a = Ts / h->tau_d, af = Ts / h->tau_fast;
    if (!h->df_valid) { h->df = d; h->dff = d; h->df_valid = 1; }
    else { h->df = h->df + a * (d - h->df); h->dff = h->dff + af * (d - h->dff); }
    ring_push(&h->hist, h->df);
    return ring_full(&h->hist) && FABS(h->df - ring_get(&h->hist, 0)) > h->thr;
}

static sma_real mean_range(const sma_ring_t *r, int i0, int i1)
{
    sma_real s = 0;
    for (int i = i0; i < i1; i++) s += ring_get(r, i);
    return s / (sma_real)(i1 - i0);
}

static sma_real std_all(const sma_ring_t *r)
{
    sma_real m = mean_range(r, 0, r->n), s = 0;
    for (int i = 0; i < r->n; i++) { sma_real d = ring_get(r, i) - m; s += d * d; }
    return SQRT(s / (sma_real)r->n);
}

/* temperature at the first time the fraction (rec[i][2]-off)/den is >= x for n_persist samples in a row, minus delay (np.interp on time) */
static sma_real at_frac(sma_real (*rec)[3], int n, sma_real off, sma_real den, sma_real x, sma_real delay, int n_persist)
{
    int idx = n - 1, run = 0;
    for (int i = 0; i < n; i++) {
        sma_real f = (rec[i][2] - off) / den;
        run = (f >= x) ? run + 1 : 0;
        if (run >= n_persist) { idx = i - n_persist + 1; break; }
    }
    sma_real tq = rec[idx][0] - delay;
    if (tq <= rec[0][0]) return rec[0][1];
    if (tq >= rec[n - 1][0]) return rec[n - 1][1];
    int lo = 0, hi = n - 1;                      /* binary search for the bracketing samples */
    while (hi - lo > 1) { int mid = (lo + hi) / 2; if (rec[mid][0] <= tq) lo = mid; else hi = mid; }
    sma_real w = (tq - rec[lo][0]) / (rec[hi][0] - rec[lo][0]);
    return rec[lo][1] + w * (rec[hi][1] - rec[lo][1]);
}

static void home_finish(sma_home_t *h)
{
    const sma_cfg_t *c = &h->ctl->c;
    sma_real L = h->ctl->w.L * (sma_real)1e3;
    int n_persist = (int)(0.5 / c->Ts); if (n_persist < 1) n_persist = 1;
    sma_real d_e2 = h->cool[h->n_cool - 1][2];
    sma_cal_t *cal = &h->cal;
    cal->eps_ext = d_e2 / L;
    cal->eps_con = h->d_c / L;
    sma_real xh = minr(maxr((sma_real)0.005, h->thr / maxr(FABS(h->d_c - h->d_e0), (sma_real)1e-9)), (sma_real)0.25);
    sma_real xc = minr(maxr((sma_real)0.005, h->thr / maxr(FABS(d_e2 - h->d_c), (sma_real)1e-9)), (sma_real)0.25);
    cal->As = at_frac(h->heat, h->n_heat, h->d_e0, h->d_c - h->d_e0, xh, h->tau_fast, n_persist);
    cal->Af = at_frac(h->heat, h->n_heat, h->d_e0, h->d_c - h->d_e0, 1 - xh, h->tau_fast * LOG(1 / xh), n_persist);
    cal->Tfs = at_frac(h->cool, h->n_cool, h->d_c, d_e2 - h->d_c, xc, 0, n_persist);
    cal->Tfe = at_frac(h->cool, h->n_cool, h->d_c, d_e2 - h->d_c, 1 - xc, 0, n_persist);
    cal->T_peak = h->T_peak;
    cal->t_heat = (sma_real)h->n_heat * c->Ts;
    cal->t_cool = (sma_real)h->n_cool * c->Ts;
    h->done = 1;
    h->cal_valid = 1;
    if (cal->Tfe < c->T_amb_est + h->preload_margin_K) { h->fault = SMA_F_PRELOAD_TOO_LOW; h->cal_valid = 0; }
}

sma_real sma_home_step(sma_home_t *h, sma_real T_meas, sma_real d_meas_mm)
{
    sma_ctl_t *ct = h->ctl;
    sma_cfg_t *c = &ct->c;
    h->t += c->Ts;
    int moving = home_moving(h, d_meas_mm);
    h->T_peak = maxr(h->T_peak, T_meas);
    sma_real R = sma_resistance(&ct->w, sma_ctl_r_from_disp(ct, d_meas_mm));
    sma_real P_max = ct->w.Imax * ct->w.Imax * R;
    sma_real Q = ct->w.h_conv * ct->w.S * maxr(T_meas - c->T_amb_est, 0);
    if (h->phase == SMA_H_WAIT) {
        ring_push(&h->wait_T, T_meas);
        ring_push(&h->wait_d, d_meas_mm);
        int steady = 0;
        if (ring_full(&h->wait_T)) {
            int n = h->wait_T.n, hh = n / 2;
            sma_real se_d = maxr(std_all(&h->wait_d), (sma_real)1e-6) * SQRT((sma_real)2.0 / (sma_real)hh);
            steady = FABS(mean_range(&h->wait_T, hh, n) - mean_range(&h->wait_T, 0, hh)) < (sma_real)0.3 && T_meas < h->T_cold_max
                     && FABS(mean_range(&h->wait_d, hh, n) - mean_range(&h->wait_d, 0, hh)) < maxr((sma_real)0.02, (sma_real)4.0 * se_d);
        }
        if (steady) {
            c->T_amb_est = mean_range(&h->wait_T, 0, h->wait_T.n);
            sma_real s_f = std_all(&h->wait_d) * SQRT(c->Ts / ((sma_real)2.0 * h->tau_d));
            h->thr = maxr(h->thr, (sma_real)5.0 * SQRT((sma_real)2.0) * s_f);
            h->d_e0 = h->df; h->goal = T_meas; ct->Tcmd_valid = 0; h->phase = SMA_H_HEAT;
        } else if (h->t > h->t_wait_max) {
            h->fault = SMA_F_NOT_COLD; h->done = 1;
        }
        return 0;
    }
    if (h->phase == SMA_H_HEAT) {
        if (h->k_dec++ % SMA_REC_DECIM == 0) {
            if (h->n_heat >= SMA_HEAT_MAX) { h->fault = SMA_F_BUFFER_FULL; h->done = 1; return 0; }
            h->heat[h->n_heat][0] = h->t; h->heat[h->n_heat][1] = T_meas; h->heat[h->n_heat][2] = h->dff; h->n_heat++;
        }
        sma_real stroke = FABS(h->df - h->d_e0);
        h->still = moving ? 0 : h->still + c->Ts;
        if (T_meas >= c->T_abs_max - 3 || (stroke > h->min_stroke_mm && h->still >= h->confirm_s)) {
            if (stroke < h->min_stroke_mm) { h->fault = SMA_F_NO_MOTION_ON_HEATING; h->done = 1; return 0; }
            h->d_c = h->df; h->phase = SMA_H_COOL; h->still = 0; h->t_cool0 = h->t; h->k_dec = 0;
            return 0;
        }
        h->Tf = h->Tf_valid ? h->Tf + c->Ts / (sma_real)1.0 * (T_meas - h->Tf) : T_meas;
        h->Tf_valid = 1;
        ring_push(&h->Tf_hist, h->Tf);
        if (ring_full(&h->Tf_hist)) {
            sma_real rate = (ring_get(&h->Tf_hist, h->Tf_hist.n - 1) - ring_get(&h->Tf_hist, 0)) / (sma_real)2.0;
            h->I_rate += c->Ts * ct->w.mc / h->tau_rate * (h->ramp - rate);
        }
        sma_real P_un = Q + ct->w.mc * h->ramp + h->I_rate;
        sma_real P = clipr(P_un, 0, P_max);
        if (P != P_un) h->I_rate += P - P_un;
        return T_meas <= c->T_abs_max ? SQRT(P / R) : (sma_real)0;
    }
    /* cool */
    if (h->k_dec++ % SMA_REC_DECIM == 0) {
        if (h->n_cool >= SMA_COOL_MAX) { h->fault = SMA_F_BUFFER_FULL; h->done = 1; return 0; }
        h->cool[h->n_cool][0] = h->t; h->cool[h->n_cool][1] = T_meas; h->cool[h->n_cool][2] = h->df; h->n_cool++;
    }
    ring_push(&h->cool_df, h->df);
    ring_push(&h->cool_T, T_meas);
    int started = FABS(h->df - h->d_c) > h->thr;
    int drifting = ring_full(&h->cool_df) && FABS(h->df - ring_get(&h->cool_df, 0)) > h->thr;
    h->still = (moving || drifting || !started) ? 0 : h->still + c->Ts;
    h->cold = (!started && T_meas - c->T_amb_est < h->cool_T_margin) ? h->cold + c->Ts : 0;
    int settled = ring_full(&h->cool_T) && FABS(T_meas - ring_get(&h->cool_T, 0)) < (sma_real)0.2;
    if (started && h->still >= h->cool_still_s && (T_meas - c->T_amb_est < h->cool_T_margin || (h->still >= h->cool_still_max_s && settled)))
        home_finish(h);
    else if (h->cold >= h->cold_no_ext_s || h->t - h->t_cool0 > h->t_cool_max) { h->fault = SMA_F_NO_EXTENSION; h->done = 1; }
    return 0;
}

/* ================================================================== safety monitor (safety.py) */
void sma_safety_init(sma_safety_t *s, const sma_wire_t *w, const sma_meff_t *M, sma_real Ts, const sma_safety_cfg_t *cfg)
{
    (void)M;
    s->c = *cfg;
    s->Ts = Ts; s->mc = w->mc; s->mL = w->mL; s->S = w->S; s->Imax = w->Imax;
    s->mcs[0] = s->mc * cfg->mc_min_frac; s->mcs[1] = s->mc * cfg->mc_max_frac; s->mL_max = s->mL * cfg->lat_max_frac;
    s->R_min = (sma_real)1e30; s->R_max = (sma_real)-1e30;
    for (int i = 0; i <= 10; i++) {
        sma_real R = sma_resistance(w, (sma_real)i / (sma_real)10.0);
        s->R_min = minr(s->R_min, R); s->R_max = maxr(s->R_max, R);
    }
    s->tripped = 0; s->reason = SMA_S_OK; s->t = 0; s->armed = 0;
    s->T_lo = s->T_lo_f = s->T_hi = s->Tamb0 = 0;
    s->I_prev = 0; s->t_Imis = 0; s->t_hot = 0; s->d_hot0 = 0; s->d_hot_valid = 0;
    s->T_prev_valid = 0; s->beat_valid = 0; s->n_same = 0; s->T1_prev = s->T2_prev = 0; s->beat = 0;
    s->stroke_valid = 0; s->stroke_mm = s->Af_max = 0;
}

static void safety_trip(sma_safety_t *s, sma_trip_t why)
{
    if (!s->tripped) { s->tripped = 1; s->reason = why; }
}

void sma_safety_arm(sma_safety_t *s, sma_real T1, int have_T2, sma_real T2)
{
    const sma_safety_cfg_t *c = &s->c;
    if (!(c->T_range_lo < T1 && T1 < c->T_range_hi) || (have_T2 && !(c->T_range_lo < T2 && T2 < c->T_range_hi))) { safety_trip(s, SMA_S_SENSOR_RANGE); return; }
    if (have_T2 && FABS(T1 - T2) > c->T_disagree) { safety_trip(s, SMA_S_TC_DISAGREE); return; }
    s->Tamb0 = T1;
    s->T_lo = T1 - c->dT_amb; s->T_hi = T1 + c->dT_amb;
    for (int k = 0; k < 2; k++) { s->obs[k].hi = s->T_hi; s->obs[k].lo = s->T_lo; s->obs[k].lof = s->T_lo; s->obs[k].lat = s->mL_max; }
    s->T_lo_f = s->T_lo;
    s->armed = 1;
}

void sma_safety_set_cal(sma_safety_t *s, sma_real stroke_mm, sma_real Af_max)
{
    s->stroke_valid = 1; s->stroke_mm = stroke_mm; s->Af_max = Af_max;
}

sma_real sma_safety_step(sma_safety_t *s, sma_real I_cmd, sma_real I_meas, sma_real T1, int have_T2, sma_real T2, int have_d, sma_real d,
                         int have_V, sma_real V, int have_beat, long beat, int *relay_closed)
{
    const sma_safety_cfg_t *c = &s->c;
    sma_real dt = s->Ts;
    s->t += dt;
    if (!s->armed && !s->tripped) sma_safety_arm(s, T1, have_T2, T2);
    if (s->tripped) { *relay_closed = 0; return 0; }
    sma_real Tm = have_T2 ? maxr(T1, T2) : T1;

    if (!(c->T_range_lo < T1 && T1 < c->T_range_hi) || (have_T2 && !(c->T_range_lo < T2 && T2 < c->T_range_hi))) safety_trip(s, SMA_S_SENSOR_RANGE);
    if (Tm > c->T_trip) safety_trip(s, SMA_S_OVER_TEMPERATURE);

    if (have_T2) {
        int fast = 0;
        if (s->T_prev_valid) fast = maxr(FABS(T1 - s->T1_prev), FABS(T2 - s->T2_prev)) / dt > (sma_real)3.0;
        if (FABS(T1 - T2) > (fast ? c->T_disagree_dyn : c->T_disagree)) safety_trip(s, SMA_S_TC_DISAGREE);
    }
    s->T1_prev = T1; s->T2_prev = have_T2 ? T2 : T1; s->T_prev_valid = 1;

    sma_real P_hi = I_meas * I_meas * s->R_max, P_lo = I_meas * I_meas * s->R_min;
    sma_real Tamb_hi = s->Tamb0 + c->dT_amb, Tamb_lo = s->Tamb0 - c->dT_amb;
    for (int k = 0; k < 2; k++) {
        sma_obs_t *o = &s->obs[k];
        sma_real mc = s->mcs[k];
        o->hi += dt * (P_hi - c->h_min * s->S * (o->hi - Tamb_hi)) / mc;
        sma_real T_old = o->lo;
        sma_real T_new = T_old + dt * (P_lo - c->h_max * s->S * (T_old - Tamb_lo)) / mc;
        if (T_new > c->T_lat_lo && T_new > T_old && o->lat > 0) {
            sma_real dE = mc * (T_new - maxr(T_old, c->T_lat_lo));
            sma_real take = minr(dE, o->lat);
            o->lat -= take;
            T_new -= take / mc;
        } else if (T_new < c->T_lat_lo - 5) {
            o->lat = s->mL_max;
        }
        o->lo = T_new;
        o->lof += dt / c->tc_lag_max * (o->lo - o->lof);
    }
    s->T_hi = maxr(s->obs[0].hi, s->obs[1].hi);
    s->T_lo = minr(s->obs[0].lo, s->obs[1].lo);
    s->T_lo_f = minr(s->obs[0].lof, s->obs[1].lof);
    if (T1 > s->T_hi + c->band) safety_trip(s, SMA_S_TC_READS_HIGH);
    if (T1 < s->T_lo_f - c->band) safety_trip(s, SMA_S_TC_READS_LOW);

    if (s->stroke_valid && have_d) {
        if (T1 > s->Af_max + c->stall_dT) {
            if (!s->d_hot_valid) { s->d_hot0 = d; s->t_hot = 0; s->d_hot_valid = 1; }
            if (I_meas > (sma_real)0.3) s->t_hot += dt;
            if (FABS(d - s->d_hot0) > c->stall_frac * s->stroke_mm) { s->d_hot0 = d; s->t_hot = 0; }
            if (s->t_hot > c->t_stall) safety_trip(s, SMA_S_STALL);
        } else {
            s->d_hot_valid = 0; s->t_hot = 0;
        }
    }

    if (I_meas > s->I_prev * (1 + c->I_tol_rel) + c->I_tol_abs) {
        s->t_Imis += dt;
        if (s->t_Imis > c->t_I_mismatch) safety_trip(s, SMA_S_CURRENT_MISMATCH);
    } else {
        s->t_Imis = 0;
    }

    if (have_V && I_meas > (sma_real)0.3) {
        sma_real R = V / I_meas;
        if (!(s->R_min * (1 - c->R_tol) < R && R < s->R_max * (1 + c->R_tol))) safety_trip(s, SMA_S_RESISTANCE);
    }

    if (have_beat) {
        s->n_same = (s->beat_valid && beat == s->beat) ? s->n_same + 1 : 0;
        s->beat = beat; s->beat_valid = 1;
        if (s->n_same >= c->n_freeze) safety_trip(s, SMA_S_WATCHDOG);
    }

    if (s->tripped) { *relay_closed = 0; return 0; }
    s->I_prev = minr(maxr(I_cmd, 0), s->Imax);
    *relay_closed = 1;
    return s->I_prev;
}
