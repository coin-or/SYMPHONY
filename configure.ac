#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006-2022 Lehigh University. All Rights Reserved.           #
#                                                                           #
# This software is licensed under the Eclipse Public License. Please see    #
# accompanying file for terms.                                              #
#                                                                           #
#===========================================================================#

## $Id: configure.ac 369 2006-06-14 01:19:26Z menal $

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([SYMPHONY],[devel],[https://github.com/coin-or/SYMPHONY/issues/new],[coin-or-symphony],[https://github.com/coin-or/SYMPHONY])

AC_COPYRIGHT([
/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2022 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(include/symphony.h)

# Do some initialization work (version numbers, change prefix default, ...)

AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                        Headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h
AC_COIN_CHECK_MATH_HDRS

# Check for pwd.h
AC_CHECK_HEADERS([pwd.h])

#############################################################################
#                    COIN-OR components                                     #
#############################################################################

m4_foreach_w([myvar],[CoinUtils Osi Cgl Clp],
  [ AC_COIN_CHK_PKG(myvar,[SymphonyLib OsiSymphonyLib])
    if test $coin_has_[]m4_tolower(myvar) != yes ; then
      AC_MSG_ERROR([Required package myvar is not available.])
    fi
  ])

AC_COIN_CHK_PKG(OsiClp,[SymphonyLib],[osi-clp])
AC_COIN_CHK_PKG(OsiDyLP,[SymphonyLib],[osi-dylp])
AC_COIN_CHK_PKG(OsiVol,[SymphonyLib],[osi-vol])
AC_COIN_CHK_PKG(OsiCpx,[SymphonyLib],[osi-cplex])
AC_COIN_CHK_PKG(OsiGlpk,[SymphonyLib],[osi-glpk])
AC_COIN_CHK_PKG(OsiSpx,[SymphonyLib],[osi-soplex])
AC_COIN_CHK_PKG(OsiXpr,[SymphonyLib],[osi-xpress])
AC_COIN_CHK_PKG(OsiGrb,[SymphonyLib],[osi-gurobi])

AC_COIN_CHK_PKG(Glpk,[SymphonyLib],[coinglpk])

# Data Packages

AC_COIN_CHK_PKG(OsiTests,[SymUnitTest],[osi-unittests])
if test $coin_has_ositests != yes ; then
  AC_MSG_ERROR([Required package OsiTests is not available.])
fi
AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],,dataonly)
AC_COIN_CHK_PKG(Miplib3,,[coindatamiplib3],,dataonly)

# Optional package for command line edit support

AC_COIN_CHK_GNU_READLINE([SymphonyLib])

#############################################################################
#                      Check for the LP solver                              #
#############################################################################

AC_MSG_CHECKING([for SYMPHONY default solver])
AC_ARG_WITH([lp-solver],
[AS_HELP_STRING([--with-lp-solver@<:@=lpsolver@:>@],
                [specify the LP solver in small
                letters(default lpsolver=clp)])],
		[sym_lp_solver=$withval],[sym_lp_solver=clp])

found_package=true
case $sym_lp_solver in 
   clp)
    AC_MSG_RESULT(Clp) 
     if [test $coin_has_clp = unavailable || test $coin_has_clp = skipping]; \
     then
         found_package=false
     fi
   ;;
   dylp)
    AC_MSG_RESULT(DyLP) 
     if [test $coin_has_dylp = unavailable]; \
     then
         found_package=false
     fi
   ;;
   vol)
    AC_MSG_RESULT(Vol) 
     if [test $coin_has_vol = unavailable]; \
     then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cpx = false; then 
         found_package=false
     fi
   ;; 
   glpk)
     AC_MSG_RESULT(Glpk)
     if test $coin_has_glpk = false; then 
         found_package=false
     fi
   ;; 
   soplex)
         AC_MSG_RESULT(Soplex)
     if test $coin_has_spx = false; then 
         found_package=false
     fi
   ;; 
   xpress)
         AC_MSG_RESULT(Xpress)
     if test $coin_has_xpr = false; then 
         found_package=false
     fi
   ;; 
   gurobi)
         AC_MSG_RESULT(Gurobi)
     if test $coin_has_grb = false; then 
         found_package=false
     fi
   ;; 
  *)
    AC_MSG_ERROR(Unknown LP solver!)
  ;;
esac

if test $found_package = false; then 
  AC_MSG_ERROR(Couldn't find the $sym_lp_solver package!)
fi

AM_CONDITIONAL(CLP_LP_SOLVER,[test $sym_lp_solver = clp &&
test $coin_has_clp = yes])
AM_CONDITIONAL(DYLP_LP_SOLVER,[test $sym_lp_solver = dylp &&
test $coin_has_dylp = yes])
AM_CONDITIONAL(VOL_LP_SOLVER,[test $sym_lp_solver = vol &&
test $coin_has_vol = yes])
AM_CONDITIONAL(CPLEX_LP_SOLVER,[test $sym_lp_solver = cplex &&
test $coin_has_cpx = yes])
AM_CONDITIONAL(GLPK_LP_SOLVER,[test $sym_lp_solver = glpk &&
test $coin_has_glpk = yes])
AM_CONDITIONAL(SOPLEX_LP_SOLVER,[test $sym_lp_solver = soplex &&
test $coin_has_spx = yes])
AM_CONDITIONAL(XPRESS_LP_SOLVER,[test $sym_lp_solver = xpress &&
test $coin_has_xpr = yes])
AM_CONDITIONAL(GRB_LP_SOLVER,[test $sym_lp_solver = gurobi &&
test $coin_has_grb = yes])

AM_CONDITIONAL(COIN_CXX_IS_CL,[test $enable_msvc = yes])

#############################################################################
#                       Add the options                                     #
#############################################################################

#---------------------------------------------------------------------------#
# USE_CGL_CUTS
#---------------------------------------------------------------------------#
AC_ARG_ENABLE([cgl-cuts],
[AS_HELP_STRING([--disable-cgl-cuts],
                [disable generic cut generation])],
[case "${enableval}" in
   yes) sym_cut_generation=true
   ;;
   no)  sym_cut_generation=no
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-cgl-cuts)
   ;;
esac],
[sym_cut_generation=true])

AM_CONDITIONAL(USE_CGL_CUTS, [test $sym_cut_generation = true])

#---------------------------------------------------------------------------#
#USE_GMPL
#---------------------------------------------------------------------------#
AC_MSG_CHECKING([whether to compile with gmpl])
AC_ARG_WITH([gmpl],
[AS_HELP_STRING([--with-gmpl],
                [compile with GMPL reader])],
		[sym_gmpl=$withval],[sym_gmpl=none])

if test $sym_gmpl = yes; then 
   AC_MSG_RESULT([yes])
   if test $coin_has_glpk = false; then 
      AC_MSG_ERROR(Couldn't find GLPK package to compile in GMPL feature)
   fi     
else
   if test $coin_has_glpk = true && test $sym_gmpl != no; then 
     AC_MSG_RESULT([yes])
     sym_gmpl=yes     
   else 
     AC_MSG_RESULT([no])
   fi
fi

AM_CONDITIONAL(USE_GMPL, [test "$sym_gmpl" = yes])

#---------------------------------------------------------------------------#
# SENSITIVITY_ANALYSIS
#---------------------------------------------------------------------------#

#AC_MSG_CHECKING([whether we want to compile in sensitivity analysis features])
AC_ARG_ENABLE([sensitivity-analysis],
[AS_HELP_STRING([--enable-sensitivity-analysis],
                [compile in the  sensitivity analysis features])],
[case "${enableval}" in
   yes) sym_sensitivity=true
   ;;
   no)  sym_sensitivity=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-sensitivity-analysis)
   ;;
esac],
[sym_sensitivity=true])
#if test $sym_sensitivity = true; then
#  AC_MSG_RESULT([yes])
#else
#  AC_MSG_RESULT([no])
#fi
AM_CONDITIONAL(SENSITIVITY_ANALYSIS, [test $sym_sensitivity = true])

#---------------------------------------------------------------------------#
# ROOT_NODE_ONLY
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([root-only],
[AS_HELP_STRING([--enable-root-only],
                [option to only process the root node])],
[case "${enableval}" in
   yes) sym_root_only=true
   ;;
   no)  sym_root_only=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-root-only)
   ;;
esac],
[sym_root_only=false])

AM_CONDITIONAL(ROOT_NODE_ONLY, [test $sym_root_only = true])

#---------------------------------------------------------------------------#
# FIXME CCMALLOC
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# COMPILE_FRAC_BRANCHING
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([frac-branching],
[AS_HELP_STRING([--enable-frac-branching],
                [whether to compile in the fractional branching option])],
[case "${enableval}" in
   yes) sym_frac_branch=true
   ;;
   no)  sym_frac_branch=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-frac-branching)
   ;;
esac],
[sym_frac_branch=false])

AM_CONDITIONAL(COMPILE_FRAC_BRANCHING, [test $sym_frac_branch = true])


#---------------------------------------------------------------------------#
# DO_TESTS
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([tests],
[AS_HELP_STRING([--enable-tests],
                [whether to perform additional sanity checks (for debugging purposes)])],
[case "${enableval}" in
   yes) sym_tests=true
   ;;
   no)  sym_tests=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-tests)
   ;;
esac],
[sym_tests=false])

AM_CONDITIONAL(DO_TESTS, [test $sym_tests = true])

#---------------------------------------------------------------------------#
# TM_BASIS_TESTS
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([tm-tests],
[AS_HELP_STRING([--enable-tm-tests],
                [more tests...])],
[case "${enableval}" in
   yes) sym_tm_tests=true
   ;;
   no)  sym_tm_tests=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-tm-tests)
   ;;
esac],
[sym_tm_tests=false])

AM_CONDITIONAL(TM_BASIS_TESTS, [test $sym_tm_tests = true])

#---------------------------------------------------------------------------#
# TRACE_PATH 
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([trace-path],
[AS_HELP_STRING([--enable-trace-path],
                [additional debugging options])],
[case "${enableval}" in
   yes) sym_trace_path=true
   ;;
   no)  sym_trace_path=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-trace-path)
   ;;
esac],
[sym_trace_path=false])

AM_CONDITIONAL(TRACE_PATH, [test $sym_trace_path = true])

#---------------------------------------------------------------------------#
# CHECK_CUT_VALIDITY
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([cut-check],
[AS_HELP_STRING([--enable-cut-check],
                [additional debugging options])],
[case "${enableval}" in
   yes) sym_cut_check=true
   ;;
   no)  sym_cut_check=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-cut-check)
   ;;
esac],
[sym_cut_check=false])

AM_CONDITIONAL(CHECK_CUT_VALIDITY, [test $sym_cut_check = true])

#---------------------------------------------------------------------------#
# STATISTICS
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([statistics],
[AS_HELP_STRING([--enable-statistics],
                [additional statistics])],
[case "${enableval}" in
   yes) sym_statistics=true
   ;;
   no)  sym_sttistics=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-statistics)
   ;;
esac],
[sym_statistics=false])

AM_CONDITIONAL(STATISTICS, [test $sym_statistics = true])

#---------------------------------------------------------------------------#
# PSEUDO_COSTS
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([pseudo-costs],
[AS_HELP_STRING([--enable-pseudo-costs],
                [enable some experimental pseudo-cost branching stuff])],
[case "${enableval}" in
   yes) sym_pseudo_costs=true
   ;;
   no)  sym_pseudo_costs=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-pseudo-costs)
   ;;
esac],
[sym_pseudo_costs=false])

AM_CONDITIONAL(PSEUDO_COSTS, [test $sym_pseudo_costs = true])

#---------------------------------------------------------------------------#
# DRAWGRAPH
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([draw-graph],
[AS_HELP_STRING([--enable-draw-graph],
                [enable IGD graph drawing application])],
[case "${enableval}" in
   yes) sym_graph=true
   ;;
   no)  sym_graph=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-graph)
   ;;
esac],
[sym_graph=false])

AM_CONDITIONAL(DRAWGRAPH, [test $sym_graph = true])

#---------------------------------------------------------------------------#
# SIGNAL HANDLER
#---------------------------------------------------------------------------#

AC_ARG_ENABLE([signal-handler],
[AS_HELP_STRING([--disable-signal-handler],
                [disable signal handler])],
[case "${enableval}" in
   yes) sym_sig_catch=true
   ;;
   no)  sym_sig_catch=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --disable-signal-handler)
   ;;
esac],
[sym_sig_catch=true])

AM_CONDITIONAL(SIGHANDLER, [test $sym_sig_catch = true])

#############################################################################
#              Check the OS and add default defs                            #
#############################################################################

SYMDEFS='-DHAS_RANDOM -DHAS_SRANDOM'

ARCH=
case $host in
  *-cygwin*) 
       if test "$enable_doscompile" = mingw; then 
           SYMDEFS="$SYMDEFS -D__MNO_CYGWIN" 
       else 
           SYMDEFS="$SYMDEFS -D__CYGWIN" 
       fi
       ARCH=CYGWIN ;;
  *-mingw*)
       SYMDEFS="$SYMDEFS -D__MNO_CYGWIN"
       ARCH=CYGWIN ;;
  *-darwin*) 
       SYMDEFS="$SYMDEFS -D__DARWIN"
       ARCH=DARWIN ;;
  *-linux-*)
       SYMDEFS="$SYMDEFS -D__LINUX" 
       case $build in 	     
	     ppc*|powerpc*)  
	     ARCH=LINUXPPC ;;
	     *)
             ARCH=LINUX ;;
       esac
       ;;
  *-alpha-*)
       SYMDEFS="$SYMDEFS -D__ALPHA"
       ARCH=ALPHA ;;
  *-aix-*)
       SYMDEFS="$SYMDEFS -D__RS6K" 
       ARCH=RS6K ;;
  *-sun-*|*solaris2*)
       case $build in
            *-i86pc-* | i386-pc*)	
                SYMDEFS="$SYMDEFS -D__X86SOL2"
	        ARCH=X86SOL2 ;;
            *-5.-*)
                np=$(shell mpstat | wc -l -gt -a)
                if test $np = ON; then 
                       SYMDEFS="$SYMDEFS -D__SUNMP"
		       ARCH=SUNMP
                else
                       SYMDEFS="$SYMDEFS -D__SUN4SOL2" 
		       ARCH=SUN4SOL2
                fi 
            ;;
       esac
       ;;
  *)
       SYMDEFS="$SYMDEFS -D__UNKNOWN_OS" ;; 
esac

#############################################################################
#                       Applications                                        #
#############################################################################

AC_MSG_CHECKING(whether to compile application library)
AC_ARG_WITH([application],
[AS_HELP_STRING([--with-application],
                [compile the application library])],
		[sym_appl=$withval],[sym_appl=no])

AM_CONDITIONAL(SYM_WITH_APPL, [test $sym_appl = yes])
AC_MSG_RESULT($sym_appl)
dnl AM_CONDITIONAL(SYM_WITH_APPL, false)
dnl AC_MSG_CHECKING(whether to compile application library)
dnl AC_ARG_WITH([application],
dnl [AS_HELP_STRING([--with-application],
dnl                 [specify the application to be compiled with lower case])],
dnl 		[sym_appl=$withval],[sym_appl=none])

dnl sym_with_appl=true
dnl case $sym_appl in 
dnl      cnrp)
dnl      ;;
dnl      match)
dnl      ;;
dnl      mckp)
dnl      ;;
dnl      mpp)
dnl      ;;
dnl      spp)
dnl      ;;
dnl      spp_with_cuts)
dnl      ;;
dnl      user)
dnl      ;;
dnl      vrp)
dnl      ;;
dnl      none)
dnl      sym_with_appl=false
dnl      ;;
dnl      *)
dnl      AC_MSG_ERROR([Unknown application is requested to be compiled!])
dnl      ;;
dnl esac
dnl AC_MSG_RESULT($sym_appl)
dnl AM_CONDITIONAL(SYM_WITH_APPL, [test $sym_with_appl = true])
dnl AM_CONDITIONAL(SYM_WITH_CNRP, [test $sym_appl = cnrp])
dnl AM_CONDITIONAL(SYM_WITH_MATCH, [test $sym_appl = match])
dnl AM_CONDITIONAL(SYM_WITH_MCKP, [test $sym_appl = mckp])
dnl AM_CONDITIONAL(SYM_WITH_MPP, [test $sym_appl = mpp])
dnl AM_CONDITIONAL(SYM_WITH_SPP, [test $sym_appl = spp])
dnl AM_CONDITIONAL(SYM_WITH_SPP_CUTS, [test $sym_appl = spp_with_cuts])
dnl AM_CONDITIONAL(SYM_WITH_USER, [test $sym_appl = user])
dnl AM_CONDITIONAL(SYM_WITH_VRP, [test $sym_appl = vrp])

#############################################################################
#                       Parallel configuration                              #
#############################################################################

AC_COIN_CHK_LIB(PVM, [SymphonyLib], [-L$PVM_ROOT/lib/$PVM_ARCH -lgpvm3 -lpvm3], [-I$PVM_ROOT/include], [], [pvm_gettid])
#AC_COIN_CHECK_USER_LIBRARY([Pvm], [PVM], [pvm3.h], [pvm_gettid])
AC_ARG_WITH([pvm],
[AS_HELP_STRING([--with-pvm],
                [compile in distributed parallel mode (assuming that 
		pvm is installed and the variable PVM_ROOT is 
		defined.)])],
		[find_pvm=$withval],[find_pvm=no])

if test x"$find_pvm" = xyes; then 
    if test x"$coin_has_pvm" != xyes; then 
        m4_ifvaln([pvm3.h],
                  [AC_CHECK_FILE([$PVM_ROOT/include/pvm3.h],
                  [pvm_header=yes])],
                  [])

	m4_ifvaln([pvm_gettid],
                  [LIBS="-L$PVM_ROOT/lib/$PVM_ARCH -lgpvm3 -lpvm3"
                   pvm_link=no
                   AC_LANG_PUSH(C)
                   AC_MSG_CHECKING([whether symbol pvm_gettid is available with PVM])
                   AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],[[pvm_gettid()]])],
                                  [AC_MSG_RESULT(yes)
                                   pvm_link=yes
                                   break],
                                  [AC_MSG_RESULT(no)])
                   AC_LANG_POP(C)])
        if test x"$pvm_header" != xyes; then
            AC_MSG_ERROR([Cannot find the pvm header files either in system path or in PVM_ROOT path:$PVM_ROOT])
        fi
        if test x"$pvm_link" != xyes ; then
            AC_MSG_ERROR([Cannot find PVM library in system paths or in PVM_ROOT path$PVM_ROOT])
        fi
        if test $pvm_link == yes && test $pvm_header == yes; then
            coin_has_pvm=yes
        fi
    fi
fi

#---------------------------------------------------------------------------#

#AC_MSG_CHECKING([whether to compile in cg])
AC_ARG_WITH([cg],
[AS_HELP_STRING([--without-cg],
                [compile without cut generator module(default=yes)])],
[in_cg=$withval],[in_cg=yes])
AM_CONDITIONAL(SYM_COMPILE_IN_CG, [test "$in_cg" = yes])

#AC_MSG_CHECKING([whether to compile in cp])
AC_ARG_WITH([cp],
[AS_HELP_STRING([--without-cp],
                [compile without cut pool module(defualt=yes)])],
[in_cp=$withval],[in_cp=yes])
AM_CONDITIONAL(SYM_COMPILE_IN_CP, [test "$in_cp" = yes])


#AC_MSG_CHECKING([whether to compile in lp])
AC_ARG_WITH([lp],
[AS_HELP_STRING([--without-lp],
                [compile without LP solver module(default=yes)])],
[in_lp=$withval],[in_lp=yes])
AM_CONDITIONAL(SYM_COMPILE_IN_LP, [test "$in_lp" = yes])

#AC_MSG_CHECKING([whether to compile in tm])
AC_ARG_WITH([tm],
[AS_HELP_STRING([--without-tm],
                [compile without tree manager module (default=yes)])],
[in_tm=$withval],[in_tm=yes])
AM_CONDITIONAL(SYM_COMPILE_IN_TM, [test "$in_tm" = yes])

#sym_def_pconfig=true;
sym_enable_parallel=false;
if [test $in_cg = no || test $in_cp = no ||  
   test $in_lp = no || test $in_tm = no]; then
#   sym_def_pconfig=false;
   if test $coin_has_pvm = no; then 
      AC_MSG_ERROR(Need PVM to enable distributed parallel execution!)
   fi
   sym_enable_parallel=true;
fi   

AM_CONDITIONAL(SYM_PARALLEL, $sym_enable_parallel)

LPEXT=
CPEXT=
TMEXT=
MASTEREXT=
LIBNAME=
SYMPDEFS="$SYMDEFS"

SYMDEFS="$SYMDEFS -DSYM_COMPILE_IN_CG -DSYM_COMPILE_IN_CP -DSYM_COMPILE_IN_LP -DSYM_COMPILE_IN_TM"

if test $in_cg = yes; then 
   LPEXT=_cg
   SYMPDEFS="$SYMPDEFS -DSYM_COMPILE_IN_CG"
fi 
if test $in_cp = yes; then 
   SYMPDEFS="$SYMPDEFS -DSYM_COMPILE_IN_CP"
   CPEXT=_cp
fi 
if test $in_lp = yes; then 
  SYMPDEFS="$SYMPDEFS -DSYM_COMPILE_IN_LP"
  TMEXT="_lp$LPEXT$CPEXT"
else
  TMEXT=$CPEXT
fi 
if test $in_tm = yes; then 
  SYMPDEFS="$SYMPDEFS -DSYM_COMPILE_IN_TM"
  MASTEREXT="_tm$TMEXT"
fi 

WHATTOMAKELIBS="libSym.la"
WHATTOINSTALLLIBS="libSym.la"
WHATTOMAKEBINS="symphony$EXEEXT"
if test $sym_appl = yes ; then 
   WHATTOMAKELIBS="$WHATTOMAKELIBS libSymAppl.la"
   WHATTOINSTALLLIBS="$WHATTOINSTALLLIBS libSymAppl.la"
   WHATTOMAKELIBS="$WHATTOMAKELIBS libSymApplStubs.la"
   WHATTOINSTALLLIBS="$WHATTOINSTALLLIBS libSymApplStubs.la"
fi

if test $sym_enable_parallel = true; then
   SYMPDEFS="$SYMPDEFS -D__PVM__"
   if test $in_cg = no; then 
      WHATTOMAKELIBS=libSym_cg.la
      if test $sym_appl = yes ; then 
      	 WHATTOMAKELIBS="$WHATTOMAKELIBS libSymAppl_cg.la"
      fi
      WHATTOMAKEBINS="symphony_cg$EXEEXT"
      fi 
   if test $in_cp = no; then 
      WHATTOMAKELIBS="$WHATTOMAKELIBS libSym_cp.la"
      if test $sym_appl = yes ; then 
      	 WHATTOMAKELIBS="$WHATTOMAKELIBS libSymAppl_cp.la"
      fi
      WHATTOMAKEBINS="$WHATTOMAKEBINS symphony_cp$EXEEXT"
   fi 
   if test $in_lp = no; then 
      WHATTOMAKELIBS="$WHATTOMAKELIBS libSym_lp$LPEXT.la"
      if test $sym_appl = yes ; then 
      	 WHATTOMAKELIBS="$WHATTOMAKELIBS libSymAppl_lp$LPEXT.la"
      fi
      WHATTOMAKEBINS="$WHATTOMAKEBINS symphony_lp$LPEXT$EXEEXT"
   fi 
   if test $in_tm = no; then 
      WHATTOMAKELIBS="$WHATTOMAKELIBS libSym_tm$TMEXT.la"
      if test $sym_appl = yes ; then 
      	 WHATTOMAKELIBS="$WHATTOMAKELIBS libSymAppl_tm$TMEXT.la"
      fi
      WHATTOMAKEBINS="$WHATTOMAKEBINS symphony_tm$TMEXT$EXEEXT"
    fi 
    MASTEREXT="_m$MASTEREXT"
    WHATTOMAKELIBS="$WHATTOMAKELIBS libSym$MASTEREXT.la"
    WHATTOMAKEBINS="$WHATTOMAKEBINS symphony$MASTEREXT$EXEEXT"
    if test $sym_appl = yes ; then 
       WHATTOMAKELIBS="$WHATTOMAKELIBS libSymAppl$MASTEREXT.la"
       WHATTOINSTALLLIBS="$WHATTOINSTALLLIBS libSymAppl$MASTEREXT.la"
       WHATTOMAKELIBS="$WHATTOMAKELIBS libSymApplStubs$MASTEREXT.la"
       WHATTOINSTALLLIBS="$WHATTOMAKELIBS libSymApplStubs$MASTEREXT.la"
    fi
fi

AC_SUBST(WHATTOMAKELIBS)
AC_SUBST(WHATTOINSTALLLIBS)
AC_SUBST(WHATTOMAKEBINS)
AC_SUBST(LPEXT)
AC_SUBST(CPEXT)
AC_SUBST(TMEXT)
AC_SUBST(MASTEREXT)

##############################################################################
#                           OPENMP defs                                     #
##############################################################################


AC_ARG_ENABLE([openmp],
[AS_HELP_STRING([--enable-openmp],
                [compile in OpenMP features])],
[case "${enableval}" in
   yes) sym_openmp=true
   ;;
   no)  sym_openmp=false
   ;;
   *) AC_MSG_ERROR(Bad value ${enableval} for --enable-openmp)
   ;;
esac],
[sym_openmp=auto])

if test $sym_enable_parallel = true; then
   if test $sym_openmp = true; then
      AC_MSG_ERROR(Cannot use both PVM and OpenMP simultaneously: Disabling OpenMP)
   fi
   sym_openmp=false
fi

if test $sym_openmp = auto; then 
   AC_MSG_CHECKING([whether compiler supports OpenMP])
   case "$CXX" in
     *clang*) 
       AC_MSG_RESULT([no])
       sym_openmp=false
       ;;
     *cl*)
       AC_MSG_RESULT([yes])
       sym_openmp=true
       ;;
     *g++*)
        sym_openmp=true
   	fopenmp_version=4.2
   	gcc_required_version=`echo ${fopenmp_version} | sed -e 's/\./0/g'`
   	get_gcc_version=`$CXX -dumpversion`
   	gcc_version=`echo ${get_gcc_version} | sed -e 's/\./0/g'`
   	major_gcc_version=`echo ${gcc_version:0:3}`
   	#AC_MSG_RESULT([$get_gcc_version])
   	if test ${major_gcc_version} -gt ${gcc_required_version}; then
       	   AC_MSG_RESULT([no])
           sym_openmp=false
        else
       	   AC_MSG_RESULT([yes])
           sym_openmp=true
	fi
        ;;
   esac
fi

AM_CONDITIONAL(OPENMP, [test $sym_openmp = true])
AC_SUBST(SYM_OPENMP_FLAGS)
AC_SUBST(SYM_OPENMP_LIBS)
if test $sym_openmp = true; then
  if test $enable_msvc = yes; then
    SYM_OPENMP_FLAGS=-openmp 
    SYM_OPENMP_LIBS=
  else
    SYM_OPENMP_FLAGS=-fopenmp
    SYM_OPENMP_LIBS=-lgomp
  fi
else
  SYM_OPENMP_FLAGS= 
  SYM_OPENMP_LIBS=
fi

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# To be used in applications: default SYMPHONY defs
AC_SUBST(SYMDEFS)
AC_SUBST(SYMPDEFS)

# Location of Applications for later use
AC_SUBST(SYMAPPLDIR)
SYMAPPLDIR=`pwd`/Applications

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/OsiSym/Makefile
                 test/Makefile
		 Examples/Makefile
		 Applications/Makefile.Applications
                 Applications/CNRP/Makefile
		 Applications/MATCH/Makefile
		 Applications/MCKP/Makefile
		 Applications/MISOCO/Makefile
		 Applications/MPP/Makefile
		 Applications/SPP/Makefile
		 Applications/SPP+CUTS/Makefile
		 Applications/USER/Makefile
		 Applications/VRP/Makefile
		 symphony.pc
		 symphony-app.pc
		 osi-sym.pc:src/OsiSym/osi-sym.pc.in])

AC_COIN_VPATH_LINK([Applications/INSTALL
                    Applications/CNRP/README
                    Applications/CNRP/INSTALL
                    Applications/CNRP/sample.cnrp
                    Applications/MATCH/README
                    Applications/MATCH/INSTALL
                    Applications/MATCH/sample.mat
                    Applications/MCKP/README
                    Applications/MCKP/INSTALL
                    Applications/MCKP/sample.mckp
                    Applications/MISOCO/README
                    Applications/MISOCO/INSTALL
                    Applications/MPP/README
                    Applications/MPP/INSTALL
                    Applications/MPP/sample.mpp
                    Applications/SPP/README
                    Applications/SPP/INSTALL
                    Applications/SPP/sample.spp
                    Applications/SPP+CUTS/README
                    Applications/SPP+CUTS/INSTALL
                    Applications/SPP+CUTS/sample.spp
                    Applications/USER/README
                    Applications/USER/INSTALL
                    Applications/USER/sample.mps
                    Applications/VRP/README
                    Applications/VRP/INSTALL
                    Applications/VRP/sample.vrp])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([include/config.h include/config_sym.h include/config_osisym.h])

AC_COIN_FINALIZE_FLAGS([SymphonyLib OsiSymphonyLib SymUnitTest])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
