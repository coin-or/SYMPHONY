###############################################################################
#                                                                             #
# This file is part of the SYMPHONY Branch, Cut, and Price Callable           #
# Library.                                                                    #
#                                                                             #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and      #
# Laci Ladanyi (ladanyi@us.ibm.com).                                          #
#                                                                             #
# (c) Copyright 2004 Ted Ralphs and Lehigh University.                        #
# All Rights Reserved.                                                        #
#                                                                             #
# The authors of this file are Menal Guzelsoy and Ted Ralphs.                 #
#                                                                             #
# This software is licensed under the Common Public License. Please see       #
# accompanying file for terms.                                                #
#                                                                             #
###############################################################################

# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
#OptLevel := -O1
OptLevel := -g -w

LIBNAME := OsiSym

LIBSRC :=
LIBSRC += OsiSymSolverInterface.cpp SymWarmStart.cpp

SYMROOT := $(HOME)/SYMPHONY

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

export ExtraIncDir  := \
    ${SymIncDir}  ${SymUserIncDir} ${OsiIncDir}  ${CoinIncDir}
export ExtraLibDir  := \
    ${SymLibDir}  ${SymUserLibDir} ${OsiLibDir}  ${CoinLibDir}
export ExtraLibName := \
    ${SymLibName} ${SymUserLibName} ${OsiLibName} ${CoinLibName}
export ExtraDefine  := \
    ${SymDefine}  ${OsiDefine}  ${CoinDefine}

export LibType OptLevel LIBNAME LIBSRC

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install libosi library clean doc

default: install

install library : % :
	(cd $(CoinDir)/Osi && $(MAKE) -f Makefile.Osi install)
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*
#	(cd $(SYMROOT) && $(MAKE))

clean doc : % :
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*


