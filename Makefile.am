#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006-2007 Lehigh University. All Rights Reserved.           #
#                                                                           #
# This software is licensed under the Common Public License. Please see     #
# accompanying file for terms.                                              #
#                                                                           #
#===========================================================================#

## $Id$

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = 

########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = CoinUtils Data/Sample Data/miplib3 Data/Netlib Osi Cgl 

if ADD_CLP
SUBDIRS += Clp
endif 

SUBDIRS += SYMPHONY FlopC++

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd SYMPHONY; $(MAKE) test

ptest: all
	cd SYMPHONY; $(MAKE) ptest

unitTest: test

punitTest: ptest

tests: all
	for dir in $(SUBDIRS); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; $(MAKE) test) \
	  fi; \
	done

unitTests: tests

.PHONY: test unitTest ptest punitTest tests unitTests

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
