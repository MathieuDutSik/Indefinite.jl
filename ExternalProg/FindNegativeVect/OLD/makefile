all: FindNegativeVect

GMPPREFIX = /opt/gmp
GSLPREFIX = /opt/gsl

#C_INCLUDE = -I$(GMPPREFIX)/include -I$(GSLPREFIX)/include
#LIBDIR    = -L$(GMPPREFIX)/lib -L$(GSLPREFIX)/lib

C_INCLUDE = 
LIBDIR    = 
LINKAGE   = -lgsl -lgslcblas -lgmp -lm

C_OPTION  = -g -Wall

FindNegativeVect.o: FindNegativeVect.c
	gcc $(C_OPTION) $(C_INCLUDE) FindNegativeVect.c -c

FindNegativeVect: FindNegativeVect.o
	gcc FindNegativeVect.o $(LIBDIR) $(LINKAGE) -o FindNegativeVect

clean:
	rm -f FindNegativeVect *.o
