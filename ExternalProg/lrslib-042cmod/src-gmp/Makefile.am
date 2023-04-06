# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -O2
AM_LDFLAGS =

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS    = glrs gredund
glrs_SOURCES    = lrs.c lrslib.c lrsgmp.c
gredund_SOURCES = redund.c lrslib.c lrsgmp.c
AM_CPPFLAGS     = -DTIMES -DGMP
