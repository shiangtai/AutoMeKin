AC_INIT([tsscds], [2018], [emilio.nunez@usc.es],[tsscds-SOURCE])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/bbfs.f90])
AC_PREFIX_DEFAULT([$HOME/AC_PACKAGE_NAME-AC_PACKAGE_VERSION])

AC_PROG_AWK
if test x"$AWK" != x"gawk" ; then
    AC_MSG_ERROR([Please install GNU gawk package before proceeding.])
fi

AC_PROG_LN_S
AC_PROG_CC
AC_PROG_F77([gfortran])
AC_PROG_FC([gfortran])

if test x"$FC" != x"gfortran" ; then
    AC_MSG_ERROR([Please install gfortran compiler before proceeding.])
fi


AM_PATH_PYTHON([2])

AC_MSG_CHECKING(for numpy python module)
$PYTHON -c "import numpy;" 2>/dev/null
if test $? == 0; then
	AC_MSG_RESULT([found])
else
	AC_MSG_ERROR([not found numpy])
fi

AC_MSG_CHECKING(for scipy python module)
$PYTHON -c "import scipy;" 2>/dev/null
if test $? == 0; then
	AC_MSG_RESULT([found])
else
	AC_MSG_ERROR([not found scipy])
fi

AC_CHECK_PROG(PARALLELS_CHECK,parallel,yes)
if test x"$PARALLELS_CHECK" != x"yes" ; then
    AC_MSG_ERROR([Please install GNU Parallel package before proceeding.])
fi

AC_CHECK_PROG(SQLITE_CHECK,sqlite3,yes)
if test x"$SQLITE_CHECK" != x"yes" ; then
    AC_MSG_ERROR([Please install sqlite3 package before proceeding.])
fi

AC_PATH_PROG(ZENITY, zenity, no)
if test x"$ZENITY" = xno; then
	AC_MSG_ERROR([Please install zenity package before proceeding.])
fi

AC_PATH_PROG(BC_CALCULATOR, bc, no)
if test x"$BC_CALCULATOR" = xno; then
	AC_MSG_ERROR([Please install bc package before proceeding.])
fi


AC_PATH_PROG(svnversioncommand, svnversion)
if test "X$svnversioncommand" = "X" || test "`$svnversioncommand -n '.'`" = "Unversioned directory" || test "`$svnversioncommand -n '.'`" = "exported"; 
then
AC_DEFINE_UNQUOTED(TSSCDS_BUILD, ["`cat share/tsscds_build`"], [repository svn version])
else
AC_DEFINE_UNQUOTED(TSSCDS_BUILD, ["`svnversion -n`"], [repository svn version])
`svnversion -n > share/tsscds_build`
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT([ src/Makefile \
	    src/HL/Makefile\
	     scripts/Makefile\
	     scripts/CATscripts/Makefile\
	     scripts/HLscripts/Makefile\
	     scripts/MOPAC_DEV/Makefile\
      	     docs/Makefile\
	     modules/Makefile\
	     examples/Makefile \
	     share/Makefile \
	     share/CAT/Makefile])

AC_MSG_NOTICE([

+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 1) "Environment Modules" must be installed in this machine. 
    Then, you may add the following to your .bash_profile: 

    module use $prefix/modules		 

 2) "Gaussian09" must be installed in this machine. 
    It must be invoked as:
 
    g09<input>output 
+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

])
