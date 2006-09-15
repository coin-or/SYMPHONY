# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

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

SUBDIRS += SYMPHONY

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd SYMPHONY; make test

unitTest: test

tests: all
	for dir in $(SUBDIRS); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; make test) \
	  fi; \
	done

unitTests: tests

.PHONY: test unitTest tests unitTests

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
