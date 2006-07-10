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
AC_COIN_CHECK_GNU_READLINE:

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Cgl)

# ToDo: find out when CBC_HAS_CLP should be defined!  For now, I
# replaced all CBC_USE by COIN_HAS #ifdef's

# Decide which one is the default solver
AC_MSG_CHECKING(for SYMPHONY default solver)
AC_ARG_WITH([sym-default-solver],
            AC_HELP_STRING([--with-sym-default-solver],
                           [specify default solver for SYMPHONY in lower case letters (if not given, clp is assumed)]),
            [sym_default_solver=$with_val],
            [sym_default_solver=clp])
AC_MSG_RESULT($sym_default_solver)

case $sym_default_solver in
  clp)
    AC_DEFINE([SYM_DEFAULT_SOLVER],"clp",
              [Default solver is chosen...])
    ;;
#  glpk)
#    AC_DEFINE([SYM_DEFAULT_SOLVER],"glpk",[GLPK is chosen...])
#    AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])
#    ;;
  *)
    AC_MSG_ERROR([Unknown solver $sym_default_solver chosen as SYMPHONY's default solver])
    ;;
esac

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

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

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/GMPL/Makefile 
                 sym_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([include/config_sym.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
