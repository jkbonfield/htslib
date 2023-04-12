# SYNOPSIS
#
#   AX_WITH_ZSTD([ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro checks whether Libzstd is installed and adds a
# --with-zstd=DIR option to override the search path.
#
#   The following output variables are set by this macro:
#
#     ZSTD_CPPFLAGS     Preprocessor flags for compiling
#     ZSTD_LDFLAGS      Linker flags for linking against the library
#     ZSTD_LIBS         Library list
#
#   The HAVE_LIBZSTD cpp variable will be defined in a working
#   libzstd was found.
#
# LICENSE
#
#   Copyright (C) 2023 Genome Research Ltd
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.
AC_DEFUN([AX_WITH_ZSTD],
[
  AC_ARG_WITH(zstd,
              AC_HELP_STRING([--with-zstd=DIR],[root DIR for zstd installation, e.g. /usr/local]),
              [_zstd_with=$withval],[_zstd_with="auto"])

  # Check if it's a working library
  _cppflags=$CPPFLAGS
  _ldflags=$LDFLAGS
  zstd_ok=no
  if test "x$_zstd_with" != "xno"
  then
    if test "$_zstd_with" != "yes"
    then
      if test -f "${_zstd_with}/include/zstd.h"
      then
        # Installation root
        ZSTD_CPPFLAGS="-I${_zstd_with}/include"
      elif test -f "${_zstd_with}/lib/zstd.h"
      then
        # Possible compiled source tree
        ZSTD_CPPFLAGS="-I${_zstd_with}/lib"
      else
        # Otherwise just try it and see
        ZSTD_CPPFLAGS="-I${_zstd_with}"
      fi
      if test -f "${_zstd_with}/lib/libzstd.a" -o -f "${_zstd_with}/lib/libzstd.so"
      then
        # Installation root or compiled source tree
        ZSTD_LDFLAGS="-L${_zstd_with}/lib"
      else
        ZSTD_LDFLAGS="-L${_zstd_with}"
      fi

    if test "$_zstd_with" != "yes" -a "$_zstd_with" != "auto"
    then
      # Explicit --with-zstd=DIR so add to RPATH too
      ZSTD_LDFLAGS="$ZSTD_LDFLAGS -Wl,-R$_zstd_with/lib"
    fi

    fi
    CPPFLAGS="$CPPFLAGS $ZSTD_CPPFLAGS"
    LDFLAGS="$LDFLAGS $ZSTD_LDFLAGS"
    AC_CHECK_LIB(zstd, ZSTD_createCStream,
        [AC_CHECK_HEADER(zstd.h,
                         [zstd_ok=yes ZSTD_LIBS="-lzstd"],
                         [zstd_ok=no])
        ])
    if test "$zstd_ok" != "yes"
    then
        AC_MSG_WARN("--with-zstd specified, but non functioning")
    fi

    # perform substitutions
    if test "$zstd_ok" = "yes"
    then
      AC_DEFINE(HAVE_ZSTD, 1,
                [Define to 1 if you have a functional libzstd.]) 
      AC_SUBST(ZSTD_CPPFLAGS)
      AC_SUBST(ZSTD_LDFLAGS)
      AC_SUBST(ZSTD_LIBS)
   else
      AC_MSG_WARN("No functioning zstd found")
      unset ZSTD_CPPFLAGS
      unset ZSTD_LDFLAGS
      unset ZSTD_LIBS
    fi

    CPPFLAGS=$_cppflags
    LDFLAGS=$_ldflags
    unset _cppflags
    unset _ldflags
  fi

  AH_TEMPLATE([HAVE_ZSTD], [Define if zstd is installed])
  AM_CONDITIONAL(HAVE_ZSTD, test "$zstd_ok" = "yes")

  # Execute the conditional expressions
  if test "$zstd_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$1],,:,[$1])
  else
     # This is the IF-NO path
     ifelse([$2],,:,[$2])
  fi

  # Tidy up
  unset zstd_ok
  unset _cppflags
  unset _ldflags
])
