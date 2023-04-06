# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -O2
AM_LDFLAGS =

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = sv sv_gmp_read
sv_SOURCES = sv.c la_support.c lll_basis.c panic.c panic.h shvec.c shvec.h lll_basis.h la_support.h
sv_gmp_read_SOURCES = sv_gmp_read.c la_support.c lll_basis.c panic.c panic.h shvec.c shvec.h lll_basis.h la_support.h
