dnl Process this file with autoconf to produce a configure script.
AC_INIT
AC_CONFIG_SRCDIR([src/Abstract/Data.C])
AC_CONFIG_AUX_DIR(config)

AC_CANONICAL_TARGET

if test "$srcdir" = .; then
   ooqp_src=src
else
   ooqp_src="$srcdir/src"
fi
AC_SUBST(ooqp_src)

dnl Checks for programs.
AC_PROG_INSTALL
dnl Use native cc first
if test -z "$CC"; then
   AC_CHECK_PROG(CC,cc,cc)
fi
dnl And optimize, unless told otherwise
if test -z "$CFLAGS"; then CFLAGS=-O; fi
AC_PROG_CC

dnl Use native CC first
if test -z "$CXX"; then
   AC_CHECK_PROG(CXX,CC,CC)
fi
dnl And optimize, unless told otherwise
if test -z "$CXXFLAGS"; then CXXFLAGS=-O; fi
AC_PROG_CXX
AC_CXX_EXPLICIT

AC_ARG_WITH(fortran, [--with-fortran: use a FORTRAN compiler (instead of
	      f2c],,
	      with_fortran=yes )

if test "$with_fortran" != "yes"; then
    FC=false
    F77=false
    dnl Any compilitions of fortran files will fail
    if test -z "$FLIBS"; then
	FLIBS=libf2c.a
    fi
    if test -z "$BLAS"; then
	BLAS=extras/blas/libblas.a
    fi
    if test -z "$MA27LIB"; then
	MA27LIB=extras/MA27/libMA27.a
    fi
    if test -z "$MA57LIB"; then
	MA57LIB=extras/MA57/libMA57.a
    fi
    dnl We cannot test the libraries and so just assume they are Ok
    blas_ok=yes; ma27_ok=yes; ma57_ok=yes;
else
    dnl The user has a fortran compiler

    dnl If FC is set, but F77 is not, assume that FC points to the F77
    dnl compiler.
    if test -z "$F77"; then
	if test -n "$FC"; then
	    F77="$FC"
	    export F77
	fi
    fi

    dnl Use native F77 first
    if test -z "$F77"; then
	AC_CHECK_PROG(F77,f77,f77)
    fi
    dnl And optimize, unless told otherwise
    if test -z "$FFLAGS"; then FFLAGS=-O; fi
    AC_PROG_F77

    dnl Checks for libraries.
    AC_F77_LIBRARY_LDFLAGS

    dnl Look for some Fortran libraries
    AC_LANG([Fortran 77])
    dnl Look for blas - abandon all hope yea who enter here
    blas_ok=no
    if test -z "$BLAS"; then
	dnl The BLAS environment variable was not set
	dnl Look for blas in some typical places
	case "$host" in
	sparc-sun-* )  
	    dnl AC_SEARCH_LIBS(dcopy,sunperf blas,,BLAS="libblas.a")
	    AC_CHECK_LIB(sunperf,dcopy,BLAS=-lsunperf)
	    echo "Checking for dcopy in -xlic_lib=sunperf (will look like I am"
	    echo "testing -lm.)"
	    if test -z "$BLAS"; then
	    oldLIBS="$LIBS"
	    LIBS="$LIBS -xlic_lib=sunperf"
	    AC_CHECK_LIB(m,dcopy,BLAS=-xlic_lib=sunperf)
	    LIBS="$oldLIBS"
	    fi
	    if test -z "$BLAS"; then
	    AC_CHECK_LIB(blas,dcopy,BLAS=-lblas)
	    fi		
	dnl end of the tests for sun
	;;
	*) 
	    AC_CHECK_LIB(blas,dcopy,BLAS=-lblas) 
	;;
	esac
	if test -n "$BLAS"; then
	    dnl     Got it!
	    blas_ok=yes
	fi
    fi
    dnl End the blas environment variable was not set
    if test "$blas_ok" = no; then
	if test -z "$BLAS"; then
	    dnl Still do not have a candidate for the blas library. 
	    dnl  Try the default
	    BLAS=extras/blas/libblas.a
	fi
	AC_CHECK_OBJ([$BLAS],dcopy,blas_ok)
    fi
    dnl Done looking for blas
    
    dnl Look for MA27 
    ma27_ok='no'
    if test -z "$MA27LIB"; then 
	AC_CHECK_LIB(MA27,ma27ad,MA27LIB=-lMA27)
	if test -n "$MA27LIB"; then
	    dnl    Got it!
	    ma27_ok='yes'
	fi
    fi
    if test "$ma27_ok" = "no"; then
	dnl Have not found MA27 yet, try the default
	if test -z "$MA27LIB"; then
	    MA27LIB=extras/MA27/libMA27.a
	fi
	AC_CHECK_OBJ([$MA27LIB],ma27ad,ma27_ok,$MA27LIB)
    fi

    echo "Checking for an old version of MA27"
    dnl Look for an old common block that has been removed!
    AC_CHECK_OBJ([$MA27LIB],ma27td,oldMA27,$MA27LIB)
    if test "$oldMA27" = "yes"; then
	AC_MSG_ERROR( [OOQP will not work with the old version of MA27.\
	    Obtain a new version from HSL.] )
    fi

    dnl Look for MA57 
    ma57_ok='no'
    if test -z "$MA57LIB"; then 
	AC_CHECK_LIB(MA57,ma57ad,MA57LIB=-lMA57,,$BLAS)
	if test -n "$MA57LIB"; then
	    dnl    Got it!
	    ma57_ok='yes'
	fi
    fi
    if test "$ma57_ok" != "yes"; then
	if test -z "$MA57LIB"; then
	    MA57LIB=extras/MA57/libMA57.a
	fi
	AC_CHECK_OBJ([$MA57LIB],ma57id,ma57_ok,"$BLAS")
    fi
fi
dnl End of --with-fortran
    
AC_SUBST(BLAS)
if test "$blas_ok" = no; then
    AC_MSG_WARN( [I could not find an appropriate BLAS library.] )
fi

AC_SUBST(MA27LIB)
if test "$ma27_ok" = "no"; then
    AC_MSG_WARN( [I could not find MA27. Set the MA27LIB environment
	    variable and try again.] )
fi

AC_SUBST(MA57LIB)
if test "$ma57_ok" = "no"; then
   AC_MSG_WARN( [I could not find MA57. To enable MA57, set] )
   AC_MSG_WARN( [the MA57LIB environment variable and try again.] )
fi

AC_LANG([C++])
AC_PROG_RANLIB
AC_LANG([C++])

dnl Checks for header files.
AC_SUBST(INCLUDES)
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.
AC_CHECK_FUNCS(getrusage)

dnl Checks for Petsc version

if test ! -d "$PETSC_DIR"; then
    petsc_enabled='no # Petsc support has been disabled'
    AC_MSG_WARN([Petsc support has been disabled.])
    AC_MSG_WARN([If you intend to use Petsc, set PETSC_DIR and])
    AC_MSG_WARN([PETSC_ARCH and run configure again.])
else
    petsc_enabled='yes'
    if test -f $PETSC_DIR/bmake/common/variables; then
dnl     This is a new version of Petsc!
        petsc_common_variables='${PETSC_DIR}/bmake/common/variables'
    else
    petsc_common_variables='${PETSC_DIR}/bmake/${PETSC_ARCH}/base_variables'
    fi
fi
AC_SUBST(petsc_enabled)
AC_SUBST(petsc_common_variables)
AC_SUBST(PETSC_DIR)

AC_CONFIG_SUBDIRS(src/Mex src/Ampl)

AC_CONFIG_FILES([GNUmakefile PetscMakefile])
AC_OUTPUT



