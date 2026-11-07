# portable C99 build of the controller, homing routine and safety monitor
CC      ?= gcc
CFLAGS  ?= -std=c99 -Wall -Wextra -pedantic -O2

all: replay replay_f

replay: sma.c sma.h replay.c
	$(CC) $(CFLAGS) sma.c replay.c -o $@ -lm

replay_f: sma.c sma.h replay.c           # single precision (MCU with a single-precision FPU)
	$(CC) $(CFLAGS) -DSMA_FLOAT sma.c replay.c -o $@ -lm

clean:
	rm -f libboard_sim.so
	rm -f replay replay_f

# board layer: 20 kHz current loop + glue, built here as a PC simulation library (power stage and sensors modelled in c/board/hal_sim.c)
board_sim: libboard_sim.so
libboard_sim.so: board/board.c board/fast_loop.c board/hal_sim.c sma.c sma.h board/*.h
	$(CC) $(CFLAGS) -DSMA_FLOAT -shared -fPIC board/board.c board/fast_loop.c board/hal_sim.c sma.c -o $@ -lm
