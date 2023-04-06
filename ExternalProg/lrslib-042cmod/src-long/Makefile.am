# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -O2
AM_LDFLAGS =

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS    = lrs_long redund_long
lrs_long_SOURCES    = lrs.c lrslib.c lrslong.c
redund_long_SOURCES = redund.c lrslib.c lrslong.c
AM_CPPFLAGS     = -DTIMES -DLONG
