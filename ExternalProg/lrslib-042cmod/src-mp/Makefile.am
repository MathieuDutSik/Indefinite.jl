# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -O2
AM_LDFLAGS =

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS    = lrs_mp redund_mp
lrs_mp_SOURCES    = lrs.c lrslib.c lrsmp.c
redund_mp_SOURCES = redund.c lrslib.c lrsmp.c
AM_CPPFLAGS     = -DTIMES
