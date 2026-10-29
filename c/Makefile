# portable C99 build of the controller, homing routine and safety monitor
CC      ?= gcc
CFLAGS  ?= -std=c99 -Wall -Wextra -pedantic -O2

all: replay replay_f

replay: sma.c sma.h replay.c
	$(CC) $(CFLAGS) sma.c replay.c -o $@ -lm

replay_f: sma.c sma.h replay.c           # single precision (MCU with a single-precision FPU)
	$(CC) $(CFLAGS) -DSMA_FLOAT sma.c replay.c -o $@ -lm

clean:
	rm -f replay replay_f
