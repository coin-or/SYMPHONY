##############################################################################
##############################################################################
#                                                                            #
# This file is part of the SYMPHONY Branch, Cut, and Price Library.          #
#                                                                            #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and     #
# Laci Ladanyi (ladanyi@us.ibm.com).                                         #
#                                                                            #
# (c) Copyright 2004 Ted Ralphs. All Rights Reserved.                        #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
##############################################################################
##############################################################################

##############################################################################
##############################################################################
#
# This makefile is for Microsoft Visual C++ usage only! In order to compile 
# this application, simply type the following command:
#
# nmake /f mckp.mak 
#
# The executable "symphony.exe" for this application will be created in 
# .\Debug directory. By default, SYMPHONY is set up to use the CLP
# optimization solver via COIN_OSI's CLP interface and to use the CGL cuts. 
# However, you are free to  specify your own settings for the executable via 
# the following variables.
# (you can download nmake.exe from  "http://download.microsoft.com/download/
# vc15/Patch/1.52/W95/EN-US/Nmake15.exe" if you need that.)
##############################################################################

##############################################################################
# The SYMPHONYROOT environment variable specifies the root directory for the 
# SYMPHONY callable library. If this file is not in the SYMPHONY root directory, 
# change this variable to the correct path. In addition, make sure that you have 
# created the SYMPHONY library (which should be in $(SYMPHONYROOT)\WIN32\Debug 
# subdirectory) before. See the README file of SYMPHONY for instructions.
##############################################################################

SYMPHONYROOT=..\..

##############################################################################
# OUTDIR variable specifies where to create the executable file, 
# "symphony.exe", the corresponding objects and the dependencies.  
##############################################################################

OUTDIR=.\Debug

##############################################################################
##############################################################################
#
# Compiling and Linking...
#
##############################################################################
##############################################################################

DEFINITIONS = "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB"

ALL_INCDIR = "$(SYMPHONYROOT)\include" /I "..\include"
ALL_LIB = $(SYMPHONYROOT)\WIN32\Debug\symphonyLib.lib

.SILENT:

CPP=cl.exe
CPPFLAGS= /nologo /MLd /W2 /GR /Gm /YX /GX /ZI /Od \
	/I $(ALL_INCDIR) /D $(DEFINITIONS) \
	/Fp"$(OUTDIR)\mckp.pch" /Fo"$(OUTDIR)\\" /Fd"$(OUTDIR)\mckp.pdb" \
	/FD /GZ /c /Tp

.c.obj: 
	$(CPP) $(CPPFLAGS) "$*.c"
	 
ALL : "$(OUTDIR)" "APPL_MESSAGE" "APPL_OBJECTS" mckp_exe 

CLEAN:
	del /Q $(OUTDIR)\*.obj
        del /Q $(OUTDIR)\mckp.exe 
        del /Q $(OUTDIR)\mckp.idb
        del /Q $(OUTDIR)\mckp.pdb
	del /Q $(OUTDIR)\mckp.pch

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

APPL_MESSAGE: 
	echo Compiling application files...

APPL_OBJECTS : \
	..\mckp_main.obj
	echo Application files compiled successfully...
	echo ...

LINK_OBJECTS= \
	$(OUTDIR)\mckp_main.obj
               	          
mckp_exe : $(LINK_OBJECTS)
	echo Linking...
	$(CPP) /nologo /Fe"$(OUTDIR)\mckp.exe" $(ALL_LIB) $**
	echo "mckp.exe" created successfully...
