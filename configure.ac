#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([komega], [1.1.0], [mkawamura@issp.u-tokyo.ac.jp])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADERS([config.h])

# Define required macros if not already defined
m4_pattern_allow([AC_MSG_FAILURE])
m4_pattern_allow([AC_MSG_WARN])

# Define macros if they don't exist
m4_ifndef([AC_MSG_FAILURE], [
  m4_define([AC_MSG_FAILURE], [AC_MSG_ERROR([$1])])
])
m4_ifndef([AC_MSG_WARN], [
  m4_define([AC_MSG_WARN], [AC_MSG_WARNING([$1])])
])

echo ""
echo "  Check compiler"
echo ""
AC_PROG_F77
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
                                 [compile with MPI (parallelization) support.
                                  If none is found, MPI is not used. Default: no])],,
            [with_mpi=no])
if test x"$with_mpi" != xno; then
  AC_PROG_FC
  FCFLAGS="$FCFLAGS -D__MPI"
else
  AC_PROG_FC
fi
echo ""
echo "  Check OpenMP"
echo ""
AC_ARG_WITH(openmp, [AS_HELP_STRING([--with-openmp],
                                 [compile with OpenMP (parallelization) support.
                                  If none is found, OpenMP is not used. Default: yes])],,
            [with_openmp=yes])
if test x"$with_openmp" != xno; then
  # Check for OpenMP support
  AC_LANG_PUSH([Fortran])
  AC_OPENMP
  if test x"$ac_cv_prog_c_openmp" = xyes; then
    FCFLAGS="$FCFLAGS $OPENMP_CFLAGS"
  else
    AC_MSG_WARN([OpenMP requested, but couldn't use OpenMP.])
  fi
  AC_LANG_POP([Fortran])
fi
echo ""
echo "  Check fortran module"
echo ""
# Set up Fortran module support
FC_MODINC="-I"
FC_MODEXT="mod"
AC_SUBST([FC_MODINC])
AC_SUBST([FC_MODEXT])
echo ""
echo "  Check library generation"
echo ""
AM_PROG_AR
LT_INIT
echo ""
echo "  Check BLAS and LAPACK"
echo ""
AX_BLAS
AX_LAPACK
#
# Shread-safe
#
AC_ARG_ENABLE(threadsafe, [AS_HELP_STRING([--enable-threadsafe],
                                 [Thread-safe mode (experimental).])],,
                           [enable_threadsafe=no])
if test x"$enable_threadsafe" != xno; then
  FCFLAGS="$FCFLAGS -D__KOMEGA_THREAD"
fi
#
# zdotc and adotu
#
AC_ARG_ENABLE(zdot, [AS_HELP_STRING([--disable-zdot],
                           [Do not use ZDOTC and ZDOTU in BLAS.])],,
                           [enable_zdot=yes])
if test x"$enable_zdot" != xyes; then
  FCFLAGS="$FCFLAGS -D__NO_ZDOT"
fi
echo ""
echo "  Output files"
echo ""
AC_CONFIG_FILES([
Makefile
src/Makefile
test/Makefile
app/src/Makefile
])

AC_OUTPUT

echo ""
echo "  SUMMARY"
echo ""
echo "             FC = ${FC}"
echo "        FCFLAGS = ${FCFLAGS}"
echo "    LAPACK_LIBS = ${LAPACK_LIBS}"
echo "      BLAS_LIBS = ${BLAS_LIBS}"
echo ""
echo "  Configure done."
echo ""
