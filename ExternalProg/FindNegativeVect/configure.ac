#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FindNegativeVect, 1.0, mathieu.dutour@gmail.com)
AM_INIT_AUTOMAKE([])
AC_OUTPUT(Makefile src/Makefile)
AM_MAINTAINER_MODE([disable])
AC_CONFIG_SRCDIR([src/FindNegativeVect.c])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([gslcblas], [main])
AC_CHECK_LIB([gsl], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h])
AC_CHECK_HEADERS([gsl/gsl_vector.h gsl/gsl_matrix.h gsl/gsl_blas.h gsl/gsl_eigen.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT
