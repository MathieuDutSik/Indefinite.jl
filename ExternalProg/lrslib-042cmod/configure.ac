#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(lrslib-042cmod, 1.0, mathieu.dutour@gmail.com)
AM_INIT_AUTOMAKE([])
AC_OUTPUT(Makefile src-long/Makefile src-gmp/Makefile src-mp/Makefile)
AM_MAINTAINER_MODE([disable])
AC_CONFIG_SRCDIR([src-gmp/lrslib.c])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgmp':
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([gmp], [main])

AC_CHECK_FUNCS([malloc])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([alarm strchr])

AC_OUTPUT
