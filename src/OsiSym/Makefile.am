#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006 Lehigh University. All Rights Reserved.                #
#                                                                           #
# This software is licensed under the Eclipse Public License. Please see    #
# accompanying file for terms.                                              #
#                                                                           #
#===========================================================================#

## $Id: Makefile.am 1152 2007-12-28 03:59:07Z andreasw $

########################################################################
#                             libOsiSym                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiSym.la

# List all source files for this library, including headers
libOsiSym_la_SOURCES = \
	OsiSymSolverInterface.cpp OsiSymSolverInterface.hpp \
	OsiSymSolverParameters.hpp \
	SymWarmStart.cpp SymWarmStart.hpp

# List all additionally required libraries
libOsiSym_la_LIBADD = $(OSISYMPHONYLIB_LFAGLS) ../libSym.la

# This is for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS =  -I`$(CYGPATH_W) $(top_srcdir)/include` \
               -I`$(CYGPATH_W) $(top_builddir)/include` 
AM_CPPFLAGS += $(OSISYMPHONYLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = \
	OsiSymSolverInterface.hpp \
	OsiSymSolverParameters.hpp \
	SymWarmStart.hpp

install-exec-local:
	$(install_sh_DATA) ../../include/config_osisym.h $(DESTDIR)$(includecoindir)/OsiSymConfig.hpp

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiSymConfig.h
