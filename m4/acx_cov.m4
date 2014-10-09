dnl autoconf macro to check libcov 
dnl
AC_DEFUN([ACX_COV],
	[acx_cov_save_LIBS="$LIBS"
	AC_LANG(C++)
	COV_LIB=""
	AC_CACHE_CHECK([for library libcov], acx_cv_cov,
		[LIBS="$acx_cov_save_LIBS -lcov"
        AC_LINK_IFELSE(
        	AC_LANG_PROGRAM([[#include <cov/cov.h>]],
				[[covNode n; bool b=n.IsLeaf();]]),
			[acx_cv_cov=yes], [acx_cv_cov=no])])
	if test "$acx_cv_cov" = yes
	then
		COV_LIB="-lcov"
	fi
	LIBS="$acx_cov_save_LIBS"
	])dnl
