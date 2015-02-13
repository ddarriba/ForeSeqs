# Copyright (C) 2005 G. Hickey <hickey@cs.dal.ca>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([seqpred],
        [1.0],
        [Diego Darriba diego.darriba@h-its.org])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2])
AC_CONFIG_MACRO_DIR([m4])

# CFLAGS="$CFLAGS -I${PWD}/libpll/include"
# CXXFLAGS="$CXXFLAGS -I${PWD}/libpll/include"
# LDFLAGS="$LDFLAGS -L${PWD}/libpll/lib"

AX_EXT

LF_CONFIGURE_CC
LF_CONFIGURE_CXX
LF_HOST_TYPE
LF_SET_WARNINGS
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CXX
AC_CONFIG_FILES([
   Makefile
   m4/Makefile
   src/Makefile
])
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB([m], [main], [], [
  AC_MSG_ERROR([Could  not find MATH library])
])

# AC_CHECK_LIB([openblas], [cblas_daxpy], [], [
#   AC_MSG_ERROR([Could  not find OPENBLAS library])
# ])

# AC_CHECK_LIB([pthread], [pthread_create])
# AC_CHECK_LIB([pll-sse3-pthreads], [pllCreateInstance])
AC_CHECK_LIB([pll-avx], [pllCreateInstance], [], [
  AC_CHECK_LIB([pll-sse3], [pllCreateInstance], [], [
    AC_MSG_ERROR([Could  not find PLL library])
])])

AC_OUTPUT
