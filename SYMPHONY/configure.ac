# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: configure.ac 369 2006-06-14 01:19:26Z menal $

# Authors:  Ted Ralphs, Menal Guzelsoy 

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([SYMPHONY],[5.1],[coin-symphony@list.coin-or.org])

AC_COPYRIGHT([
/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(include/symphony_api.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(SYMPHONY)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# Check if user wants to use readline library
AC_COIN_CHECK_GNU_READLINE
AM_CONDITIONAL(COIN_HAS_READLINE,
               [test $coin_has_readline = yes])

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Cgl)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h])#,[glp_lpx_simplex])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

#############################################################################
#                      Check for the LP solver                              #
#############################################################################

AM_CONDITIONAL(COIN_HAS_CLP,
               false)

AC_MSG_CHECKING(for SYMPHONY default solver)
if test $m4_tolower(coin_has_glpk) = true; then
   AC_MSG_RESULT(glpk)
elif test -r $m4_tolower(coin_has_cplex) = true; then
   AC_MSG_RESULT(cplex)
#FIXME! What is wrong with osl?
dnl elif test -r $coin_has_osl = true; then
dnl   AC_MSG_RESULT(osl)
elif test -r $m4_tolower(coin_has_soplex) = true; then
   AC_MSG_RESULT(soplex)
elif test -r $m4_tolower(coin_has_xpress) = true; then
   AC_MSG_RESULT(xpress)
else
   AC_MSG_RESULT(clp)
   AC_COIN_HAS_PROJECT(Clp)
fi

#############################################################################
#                       Add the options                                     #
#############################################################################

#FIXME! Should we carry this to root config file?
#---------------------------------------------------------------------------#
# USE_CGL_CUTS
#---------------------------------------------------------------------------#
AC_ARG_ENABLE([cgl-cuts],
[AC_HELP_STRING([--disable-cgl-cuts],
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
#FIXME! Check for GLPK! 
#---------------------------------------------------------------------------#
#AC_MSG_CHECKING([whether to compile in gmpl])
AC_ARG_WITH([gmpl],
[AC_HELP_STRING([--with-gmpl],
                [compile with GMPL reader])],
		[sym_gmpl=$withval],[sym_gmpl=no])
AM_CONDITIONAL(USE_GMPL, [test "$sym_gmpl" = yes])

#FIXME! NO WAY TO COMPILE OSISYM! all other packages refer to SYM instead of
# SYMPHONY
#---------------------------------------------------------------------------#
# USE_OSI_INTERFACE
#---------------------------------------------------------------------------#
AC_MSG_CHECKING([whether to compile osi interface])
AC_ARG_WITH([osi-interface],
[AC_HELP_STRING([--with-osi-interface],
                [compile with SymOsi interface])],
		[sym_osi=$withval],[])

if test "$sym_osi" = yes; then
   AC_MSG_RESULT([yes])
 else
   AC_MSG_RESULT([no])
 fi
AM_CONDITIONAL(USE_OSI_INTERFACE, [test "$sym_osi" = yes])

#---------------------------------------------------------------------------#
# SENSITIVITY_ANALYSIS
#---------------------------------------------------------------------------#

#AC_MSG_CHECKING([whether we want to compile in sensitivity analysis features])
AC_ARG_ENABLE([sensitivity-analysis],
[AC_HELP_STRING([--enable-sensitivity-analysis],
                [compile in the  sensitivity analysis features])],
[case "${enableval}" in
   yes) sym_sensitivity=true
   ;;
   no)  sym_sensitivity=false
   ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-sensitivity-analysis)
   ;;
esac],
[sym_sensitivity=false])
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
[AC_HELP_STRING([--enable-root-only],
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
[AC_HELP_STRING([--enable-frac-branching],
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
[AC_HELP_STRING([--enable-tests],
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
[AC_HELP_STRING([--enable-tm-tests],
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
[AC_HELP_STRING([--enable-trace-path],
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
[AC_HELP_STRING([--enable-cut-check],
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
[AC_HELP_STRING([--enable-statistics],
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
[AC_HELP_STRING([--enable-pseudo-costs],
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
[AC_HELP_STRING([--enable-draw-graph],
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


#############################################################################
#              Check the OS and add default defs                            #
#############################################################################

SYMDEFS='-DHAS_RANDOM -DHAS_SRANDOM -D__NONE__'

case $build in
  *-cygwin* | *-mingw*) 
       SYMDEFS="$SYMDEFS -D__CYGWIN" ;;
  *-darwin*) 
       SYMDEFS="$SYMDEFS -D__DARWIN" ;; 
  *-linux-*)
       SYMDEFS="$SYMDEFS -D__LINUX" ;; 
  *-alpha-*)
       SYMDEFS="$SYMDEFS -D__ALPHA" ;; 
  *-aix-*)
       SYMDEFS="$SYMDEFS -D__RS6K" ;; 
  *-sun-*)
       case $build in
            *-i86pc-*)	
                SYMDEFS="$SYMDEFS -D__X86SOL2" ;; 
            *-5.-*)
                SYMDEFS="$SYMDEFS -D__SUN4SOL2" 
                np=$(shell mpstat | wc -l -gt -a)
                case $np in 
                    ON)
                       SYMDEFS="$SYMDEFS -D__SUNMP";;
                esac
            ;;
       esac
       ;;
  *-sgi-*)
       SYMDEFS="$SYMDEFS -D__" ;; 
  *-hp-*)
       SYMDEFS="$SYMDEFS -D__" ;; 
  *-ibm-*)
       SYMDEFS="$SYMDEFS -D__" ;; 
  *)
       AC_MSG_ERROR([Unknown OS]) ;;
esac

#############################################################################
#                       Parallel configuration                              #
#############################################################################
AC_MSG_CHECKING([whether to compile in the parallel configuration])
AC_ARG_WITH([pvm],
[AC_HELP_STRING([--with-pvm],
                [set the root directory for pvm])],
		[PVM_ROOT=$withval],[PVM_ROOT=none])
sym_parallel=yes
case $PVM_ROOT in
     none)
    sym_parallel=false
    AC_MSG_RESULT(no)
     ;;
     *)
     m4_ifvaln([pvmalloc.h],[AC_CHECK_FILE([$PVM_ROOT/src/pvmalloc.h],[],
                 [AC_MSG_ERROR([Cannot find the pvm header files 
in $PVM_ROOT/src])])])
    AC_MSG_RESULT(yes)
    ;;
esac

AM_CONDITIONAL(SYM_PARALLEL, [test "$sym_parallel" = yes])
AC_SUBST(PVM_ROOT)

#---------------------------------------------------------------------------#

#AC_MSG_CHECKING([whether to compile in cg])
AC_ARG_WITH([cg],
[AC_HELP_STRING([--without-cg],
                [compile without cut generator module])],
[in_cg=$withval],[in_cg=yes])
AM_CONDITIONAL(COMPILE_IN_CG, [test "$in_cg" = yes])

#AC_MSG_CHECKING([whether to compile in cp])
AC_ARG_WITH([cg],
[AC_HELP_STRING([--without-cp],
                [compile without cut pool module])],
[in_cp=$withval],[in_cp=yes])
AM_CONDITIONAL(COMPILE_IN_CP, [test "$in_cp" = yes])


#AC_MSG_CHECKING([whether to compile in lp])
AC_ARG_WITH([lp],
[AC_HELP_STRING([--without-lp],
                [compile without LP solver module])],
[in_lp=$withval],[in_lp=yes])
AM_CONDITIONAL(COMPILE_IN_LP, [test "$in_lp" = yes])

#AC_MSG_CHECKING([whether to compile in tm])
AC_ARG_WITH([tm],
[AC_HELP_STRING([--without-tm],
                [compile without tree manager module])],
[in_tm=$withval],[in_tm=yes])
AM_CONDITIONAL(COMPILE_IN_TM, [test "$in_tm" = yes])


#FIXME-For now, add defaults to SYMDEFS

SYMDEFS="$SYMDEFS -DCOMPILE_IN_CG -DCOMPILE_IN_CP -DCOMPILE_IN_LP -DCOMPILE_IN_TM "

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
#AC_COIN_HAS_ASL

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

#AC_COIN_VPATH_LINK(examples/sudoku_sample.csv)
#AC_COIN_VPATH_LINK(examples/quad.mps)
#AC_COIN_VPATH_LINK(examples/quad2.mps)


#############################################################################
#                       Applications                                        #
#############################################################################
AM_CONDITIONAL(SYM_WITH_APPL, false)
AC_MSG_CHECKING(for the application to be compiled)
AC_ARG_WITH([application],
[AC_HELP_STRING([--with-application],
                [specify the application to be compiled with lower case])],
		[sym_appl=$withval],[sym_appl=none])

sym_with_appl=true
case $sym_appl in 
     cnrp)
     ;;
     match)
     ;;
     mckp)
     ;;
     mpp)
     ;;
     spp)
     ;;
     spp_with_cuts)
     ;;
     user)
     ;;
     vrp)
     ;;
     none)
     sym_with_appl=false
     ;;
     *)
     AC_MSG_ERROR([Unknown application is requested to be compiled!])
     ;;
esac
AC_MSG_RESULT($sym_appl)
AM_CONDITIONAL(SYM_WITH_APPL, [test $sym_with_appl = true])
AM_CONDITIONAL(SYM_WITH_CNRP, [test $sym_appl = cnrp])
AM_CONDITIONAL(SYM_WITH_MATCH, [test $sym_appl = match])
AM_CONDITIONAL(SYM_WITH_MCKP, [test $sym_appl = mckp])
AM_CONDITIONAL(SYM_WITH_MPP, [test $sym_appl = mpp])
AM_CONDITIONAL(SYM_WITH_SPP, [test $sym_appl = spp])
AM_CONDITIONAL(SYM_WITH_SPP_CUTS, [test $sym_appl = spp_with_cuts])
AM_CONDITIONAL(SYM_WITH_USER, [test $sym_appl = user])
AM_CONDITIONAL(SYM_WITH_VRP, [test $sym_appl = vrp])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################
# Location of SYMPHONY for later use
AC_SUBST(SYMINCDIR)
AC_SUBST(SYMSRCDIR)
SYMINCDIR=`pwd`/include
SYMSRCDIR=`pwd`/src 

# To be used in applications: default SYMPHONY defs
AC_SUBST(SYMDEFS)
AC_SUBST(COIN_HAS_CLP)

# Location of Applications for later use
AC_SUBST(SYMAPPLDIR)
SYMAPPLDIR=`pwd`/Applications

#FIXME-Debugging conditional, for now

AM_CONDITIONAL(SYM_DEBUG, [ test $coin_debug_compile = true])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
		 Applications/CNRP/Makefile
		 Applications/USER/Makefile
                 sym_addlibs.txt])

#src/GMPL/Makefile 

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([include/config_sym.h])



# Finally, we let configure write all the output...
AC_COIN_FINALIZE
