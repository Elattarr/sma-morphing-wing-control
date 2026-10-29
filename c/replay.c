/* Replay harness: reads a scenario written by scripts/export_replay.py (a logged closed-loop Python run: wire, controller and monitor
 * constants, then one row per sample with the inputs the Python code saw and the outputs it produced), runs the C controller / homing
 * and the C safety monitor on the SAME inputs and reports the largest difference.  Open loop on purpose: any deviation shows up at the
 * sample where it starts instead of being absorbed by feedback.
 *
 * usage: replay scenario.txt        output (one line, parsed by tests/test_c_port.py):
 *   rows N max_dI_ctl X max_dI_perm Y relay_mismatch K first_bad R home_done_row H cal as af tfs tfe eps_ext eps_con trip T fault F
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sma.h"

static sma_real heat_rec[SMA_HEAT_MAX][3];
static sma_real cool_rec[SMA_COOL_MAX][3];
static sma_home_t HOME;
static sma_ctl_t CTL;
static sma_safety_t SAF;

static int rd(FILE *f, double *x) { return fscanf(f, "%lf", x) == 1; }

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: replay scenario.txt\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("scenario"); return 2; }
    double v[64];
    sma_wire_t w; sma_meff_t M; sma_cfg_t cfg; sma_safety_cfg_t sc;
    sma_real sigma, K, arm_T1, arm_T2; int arm_have_T2;
    int i, n;
    char key[32];

    if (fscanf(f, "%31s", key) != 1) return 2;          /* wire */
    for (i = 0; i < 9; i++) rd(f, &v[i]);
    w.L = v[0]; w.A = v[1]; w.S = v[2]; w.mc = v[3]; w.mL = v[4]; w.rho_eM = v[5]; w.rho_eA = v[6]; w.h_conv = v[7]; w.Imax = v[8];
    if (fscanf(f, "%31s", key) != 1) return 2;          /* meff */
    for (i = 0; i < 8; i++) rd(f, &v[i]);
    M.T0 = v[0]; M.R = v[1]; M.beta = v[2]; M.eps_L = v[3]; M.h = v[4]; M.m = v[5]; M.E_A = v[6]; M.E_M = v[7];
    if (fscanf(f, "%31s", key) != 1) return 2;          /* cfg */
    for (i = 0; i < 30; i++) rd(f, &v[i]);
    cfg.Ts = v[0]; cfg.Kp = v[1]; cfg.Ki = v[2]; cfg.k_soft = v[3]; cfg.dTdt_max = v[4]; cfg.margin_up = v[5]; cfg.margin_dn = v[6];
    cfg.T_safe_above_fs = v[7]; cfg.Kr = v[8]; cfg.Ki_r = v[9]; cfg.Kr_dn = v[10]; cfg.Ki_r_dn = v[11]; cfg.tau_filt = v[12]; cfg.Kp_hold = v[13];
    cfg.Ki_hold = v[14]; cfg.drdt_up = v[15]; cfg.dn_frac = v[16]; cfg.Kaw = v[17]; cfg.db_start = v[18]; cfg.db_done = v[19]; cfg.T_amb_est = v[20];
    cfg.T_limit_above_Af = v[21]; cfg.T_abs_max = v[22]; cfg.T_limit_expand_rate = v[23]; cfg.learn_margin = v[24]; cfg.hold_step = v[25];
    cfg.hold_drift = v[26]; cfg.hold_below_As = v[27]; cfg.sigma_unc = v[28]; cfg.overshoot_exit = (int)v[29];
    if (fscanf(f, "%31s", key) != 1) return 2;          /* ctl */
    rd(f, &v[0]); rd(f, &v[1]); sigma = v[0]; K = v[1];
    if (fscanf(f, "%31s", key) != 1) return 2;          /* safety */
    for (i = 0; i < 22; i++) rd(f, &v[i]);
    sc.T_trip = v[0]; sc.T_disagree = v[1]; sc.T_disagree_dyn = v[2]; sc.h_min = v[3]; sc.h_max = v[4]; sc.dT_amb = v[5]; sc.band = v[6];
    sc.T_lat_lo = v[7]; sc.stall_dT = v[8]; sc.t_stall = v[9]; sc.stall_frac = v[10]; sc.I_tol_abs = v[11]; sc.I_tol_rel = v[12];
    sc.t_I_mismatch = v[13]; sc.R_tol = v[14]; sc.mc_min_frac = v[15]; sc.mc_max_frac = v[16]; sc.lat_max_frac = v[17]; sc.tc_lag_max = v[18];
    sc.T_range_lo = v[19]; sc.T_range_hi = v[20]; sc.n_freeze = (int)v[21];
    if (fscanf(f, "%31s", key) != 1) return 2;          /* arm */
    rd(f, &v[0]); rd(f, &v[1]); rd(f, &v[2]);
    arm_T1 = v[0]; arm_have_T2 = (int)v[1]; arm_T2 = v[2];
    if (fscanf(f, "%31s", key) != 1) return 2;          /* rows */
    rd(f, &v[0]); n = (int)v[0];

    sma_ctl_init(&CTL, &w, &M, &cfg, sigma, K);
    sma_home_init(&HOME, &CTL, heat_rec, cool_rec);
    sma_safety_init(&SAF, &w, &M, cfg.Ts, &sc);
    sma_safety_arm(&SAF, arm_T1, arm_have_T2, arm_T2);

    double max_dc = 0, max_dp = 0; int relay_bad = 0, first_bad = -1, home_done_row = -1, home_applied = 0;
    for (i = 0; i < n; i++) {
        double phase, tgt, T1, d, I_ctl, I_in, I_meas, T2, V, beat, I_perm, relay, hT2, hV, hd;
        if (!(rd(f, &phase) && rd(f, &tgt) && rd(f, &T1) && rd(f, &d) && rd(f, &I_ctl) && rd(f, &I_in) && rd(f, &I_meas) && rd(f, &hT2) && rd(f, &T2) &&
              rd(f, &hV) && rd(f, &V) && rd(f, &hd) && rd(f, &beat) && rd(f, &I_perm) && rd(f, &relay))) { fprintf(stderr, "short file at row %d\n", i); return 2; }
        sma_real Ic;
        if (phase < 0.5 || !HOME.done) {          /* a single-precision build may need a sample more than Python: keep homing until it is done */
            Ic = sma_home_step(&HOME, (sma_real)T1, (sma_real)d);
            if (HOME.done && home_done_row < 0) home_done_row = i;
        } else {
            if (!home_applied && HOME.done && HOME.cal_valid) {          /* what the host does when homing has finished */
                sma_ctl_apply_cal(&CTL, &HOME.cal);
                sma_safety_set_cal(&SAF, (HOME.cal.eps_ext - HOME.cal.eps_con) * w.L * (sma_real)1e3, HOME.cal.Af);
                home_applied = 1;
            }
            Ic = sma_ctl_step(&CTL, (sma_real)tgt, (sma_real)T1, (sma_real)d);
        }
        int relay_c = 0;
        sma_real Ip = sma_safety_step(&SAF, (sma_real)I_in, (sma_real)I_meas, (sma_real)T1, hT2 > 0.5, (sma_real)T2, hd > 0.5, (sma_real)d,
                                      hV > 0.5, (sma_real)V, 1, (long)beat, &relay_c);
        double dc = fabs((double)Ic - I_ctl), dp = fabs((double)Ip - I_perm);
        if (dc > max_dc) max_dc = dc;
        if (dp > max_dp) max_dp = dp;
        if (relay_c != (int)relay) relay_bad++;
        if (first_bad < 0 && (dc > 1e-6 || dp > 1e-6 || relay_c != (int)relay)) first_bad = i;
    }
    fclose(f);
    printf("rows %d max_dI_ctl %.3e max_dI_perm %.3e relay_mismatch %d first_bad %d home_done_row %d home_fault %d cal %d "
           "as %.9g af %.9g tfs %.9g tfe %.9g eps_ext %.12g eps_con %.12g trip %d\n",
           n, max_dc, max_dp, relay_bad, first_bad, home_done_row, (int)HOME.fault, HOME.cal_valid,
           (double)HOME.cal.As, (double)HOME.cal.Af, (double)HOME.cal.Tfs, (double)HOME.cal.Tfe, (double)HOME.cal.eps_ext, (double)HOME.cal.eps_con,
           (int)SAF.reason);
    return 0;
}
