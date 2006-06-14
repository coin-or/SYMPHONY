# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 344 2006-06-01 19:17:05Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = $(subdirs)

########################################################################
#                           Extra Targets                              #
########################################################################

#test: all
#	cd Cbc; make test

#unitTest: test

tests: all
	for dir in $(subdirs); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; make test) \
	  fi; \
	done

unitTests: tests

.PHONY: tests unitTests

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
