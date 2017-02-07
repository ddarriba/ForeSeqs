# Copyright (C) 2005 G. Hickey <hickey@cs.dal.ca>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([foreseqs],
        [1.0.3],
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
   man/Makefile
])
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB([m], [main], [], [
  AC_MSG_ERROR([Could  not find MATH library])
])

# AC_CHECK_LIB([openblas], [cblas_daxpy], [], [
#   AC_MSG_ERROR([Could  not find OPENBLAS library])
# ])

AC_CHECK_HEADERS([pll.h], [], [
    AC_CHECK_HEADERS([libpll/pll.h], [
      have_prefix=yes
    ], [
      AC_MSG_ERROR([pll.h not found])
    ])])

if test "x$have_prefix" = "xyes"; then
  AC_CHECK_HEADERS([libpll/pll_binary.h], [], [AC_MSG_ERROR([pll_binary.h not found])])
  AC_CHECK_HEADERS([libpll/pllmod_algorithm.h], [], [AC_MSG_ERROR([pllmod_algorithm.h not found])])
  AC_CHECK_HEADERS([libpll/pll_optimize.h], [], [AC_MSG_ERROR([pll_optimize.h not found])])
  AC_CHECK_HEADERS([libpll/pll_tree.h], [], [AC_MSG_ERROR([pll_tree.h not found])])
else
  AC_CHECK_HEADERS([pll_binary.h], [], [AC_MSG_ERROR([pll_binary.h not found])])
  AC_CHECK_HEADERS([pllmod_algorithm.h], [], [AC_MSG_ERROR([pllmod_algorithm.h not found])])
  AC_CHECK_HEADERS([pll_optimize.h], [], [AC_MSG_ERROR([pll_optimize.h not found])])
  AC_CHECK_HEADERS([pll_tree.h], [], [AC_MSG_ERROR([pll_tree.h not found])])
fi

AC_CHECK_LIB([pll], [pll_partition_create], [], [
  AC_MSG_ERROR([PLL library not found])
])

AC_OUTPUT
