#include <math.h>
#include "fast_loop.h"
#include "board_cfg.h"

float fl_counts_to_amps(uint16_t c)
{
    float v = (float)c * (ADC_VREF / ADC_FULL) - INA_OFFSET_V;
    return v / (INA_GAIN * SHUNT_OHM);
}

float fl_counts_to_vbus(uint16_t c)
{
    return (float)c * (ADC_VREF / ADC_FULL) / VBUS_DIV;
}

void fl_init(fl_t *f)
{
    f->I_cmd = 0; f->I_ref = 0; f->R_ff = 0.08f; f->integ = 0; f->D = 0; f->enabled = 0; f->over_ticks = 0; f->latched_over = 0;
    f->s_I = f->s_V = f->s_I2 = f->s_Iraw = 0; f->n = 0; f->I_meas = 0; f->Vbus = 5.0f;
}

void fl_set_command(fl_t *f, float I_cmd, float R_ff, int enable)
{
    f->I_cmd = I_cmd > 0 ? I_cmd : 0;
    f->R_ff = R_ff;
    f->enabled = enable && !f->latched_over;
}

float fl_tick(fl_t *f, uint16_t i_counts, uint16_t vbus_counts)
{
    float I = fl_counts_to_amps(i_counts);
    float Vb = fl_counts_to_vbus(vbus_counts);
    f->I_meas = I; f->Vbus = Vb;

    if (I > I_TRIP_A) { if (++f->over_ticks >= I_TRIP_TICKS) { f->latched_over = 1; f->enabled = 0; } } else f->over_ticks = 0;

    float D_in_effect = f->D;                              /* the duty that produced THIS sample (written one period ago) */
    float D = 0.0f;
    {                                                       /* command slew limit: a step of the command never becomes a step of the duty */
        float dI = f->I_cmd - f->I_ref;
        f->I_ref += dI > FL_SLEW_A ? FL_SLEW_A : (dI < -FL_SLEW_A ? -FL_SLEW_A : dI);
    }
    if (f->enabled && f->I_cmd > 0.0f && Vb > 1.0f) {
        float e = f->I_ref - I;
        float un = f->I_ref * f->R_ff / Vb + FL_KP * e + f->integ;      /* feed-forward + P + I */
        D = un < 0.0f ? 0.0f : (un > D_MAX ? D_MAX : un);
        if (!(un > D_MAX && e > 0.0f) && !(un < 0.0f && e < 0.0f))     /* anti-windup: no integration into a saturated output */
            f->integ += FL_KI * e;
        if (f->integ > D_MAX) f->integ = D_MAX;
        if (f->integ < -D_MAX) f->integ = -D_MAX;
    } else {
        f->integ = 0.0f;                                   /* off: no wind-up, no residual duty */
        f->I_ref = I > 0.0f && f->enabled ? I : 0.0f;       /* restart the ramp from the measured current */
    }
    f->D = D;

    /* wire voltage = mean(D Vbus) - I R_OTHER over the slow period (coil energy change over 50 ms is < 10 mV: ignored) */
    f->s_I += I; f->s_I2 += I * I; f->s_V += D_in_effect * Vb - I * R_OTHER_OHM; f->s_Iraw += Vb; f->n++;
    return D;
}

void fl_take_slow(fl_t *f, float *I_mean, float *I_rms, float *V_wire, float *Vbus_mean)
{
    float n = f->n > 0 ? (float)f->n : 1.0f;
    *I_mean = f->s_I / n;
    *I_rms = sqrtf(f->s_I2 / n);
    *V_wire = f->s_V / n;
    *Vbus_mean = f->s_Iraw / n;
    f->s_I = f->s_V = f->s_I2 = f->s_Iraw = 0; f->n = 0;
}
