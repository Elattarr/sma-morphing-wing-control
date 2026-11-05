/* Board constants of the drive in hardware/pin_table.md.  Every value is a DESIGN value, none is measured on a built board. */
#ifndef BOARD_CFG_H
#define BOARD_CFG_H

#define F_PWM_HZ        20000.0f
#define TS_FAST         (1.0f / F_PWM_HZ)          /* 50 us: one ADC sample and one duty update per PWM period */
#define TS_SLOW         0.05f                      /* 50 ms: thermocouples, laser, controller, safety monitor */
#define FAST_PER_SLOW   1000                       /* TS_SLOW / TS_FAST */

/* current path: shunt R1 -> INA240A2 -> ADC (PA0).  35 A -> 35 mV -> 1.75 V */
#define SHUNT_OHM       0.001f
#define INA_GAIN        50.0f
#define ADC_BITS        12
#define ADC_FULL        4095.0f
#define ADC_VREF        3.3f
#define INA_OFFSET_V    0.0f                       /* REF1/REF2 at ground: output is 0 V at 0 A, if the INA240 swings that low (NOT checked) */

/* bus voltage: divider into a second ADC channel (needed for the duty -> voltage estimate) */
#define VBUS_DIV        0.5f                       /* ADC volts per bus volt (2 x 10 k): 5 V -> 2.5 V */

/* wire-voltage estimate  V_wire = mean(D * Vbus) - I * R_OTHER.  R_OTHER = Q1 RDS(on) hot + shunt + coil winding + leads and clamps */
#define R_OTHER_OHM     0.012f

/* laser 4..20 mA over 150 ohm on the wire end */
#define LASER_R_OHM     150.0f
#define LASER_RANGE_MM  6.64f
#define LASER_SIGN      (-1.0f)                    /* d [mm] falls when the laser reading rises (mounting dependent: set on the bench) */

/* limits */
#define D_MAX           0.95f                      /* bootstrap capacitor needs the low side on for part of every period */
#define I_TRIP_A        40.0f                      /* software overcurrent latch (the fuse is 40 A; the real fast trip is hardware, see docs) */
#define I_TRIP_TICKS    3                          /* consecutive fast samples above I_TRIP_A */
#define DEAD_TIME_S     100e-9f                    /* HRTIM dead time, placeholder: set from QG and gate resistors */

#define FL_SLEW_A       1.0f                       /* A per fast tick: 35 A in 1.75 ms, much faster than the wire (seconds) */

/* current loop (duty = feed-forward + PI on the measured current); gains found in simulation, tests/test_board.py */
#define FL_KP           0.02f                      /* duty per ampere */
#define FL_KI           0.006f                     /* duty per ampere per fast tick */

#endif
