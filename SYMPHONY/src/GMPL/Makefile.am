# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 387 2006-06-22 23:15:50Z menal $

# Authors: Ted Ralphs, Menal Guzelsoy           

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libGmpl                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
noinst_LTLIBRARIES = libGmpl.la

# List all source files for this library, including headers
libGmpl_la_SOURCES = \
	glpavl.h \
	glpdmp.h \
	glplib.h \
	glpmpl.h \
	glprng.h \
	glpmpl1.c \
	glpmpl2.c \
	glpmpl3.c \
	glpmpl4.c \
	glpdmp.c \
	glpavl.c \
	glprng.c \
	glplib1a.c \
	glplib2.c \
	glplib3.c

# This is for libtool (on Windows)
libGmpl_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) .`
CFLAGS =
ConfigHeader = 

AM_CPPFLAGS += -DHAVE_LIBM=1 -DSTDC_HEADERS=1
# install-exec-local:
# 	echo "#ifndef HAVE_CONFIG_H" >bla
# 	echo "#define HAVE_CONFIG_H" >>bla
# 	echo "#endif" >> bla
# 	cat $(srcdir)/$(ConfigHeader) >> bla
# 	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
# 	rm -f bla

# uninstall-local:
# 	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)
