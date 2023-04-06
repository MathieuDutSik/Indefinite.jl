# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -O2
AM_LDFLAGS =

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = ISOM AUTO
ISOM_SOURCES = ISOM.c bachtools.c iotools.c isotools.c lattools.c mattools.c orbtools.c preproc.c sorttools.c types.h
AUTO_SOURCES = AUTO.c auttools.c bachtools.c iotools.c lattools.c mattools.c orbtools.c preproc.c sorttools.c types.h
