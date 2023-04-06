#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(SV, 0.2, mathieu.dutour@gmail.com)
AM_INIT_AUTOMAKE([-Wall])
AC_OUTPUT(Makefile src/Makefile)
AM_MAINTAINER_MODE([disable])
AC_CONFIG_SRCDIR([src/lll_basis.c])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lgmp':
AC_CHECK_LIB([gmp], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor sqrt])

AC_OUTPUT
