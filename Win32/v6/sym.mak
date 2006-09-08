##############################################################################
##############################################################################
#                                                                            #
# This file is part of the SYMPHONY Branch, Cut, and Price Library.          #
#                                                                            #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and     #
# Laci Ladanyi (ladanyi@us.ibm.com).                                         #
#                                                                            #
# (c) Copyright 2005 Ted Ralphs. All Rights Reserved.                        #
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
# nmake /f sym.mak 
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
# source code. If this file is not in the SYMPHONY root directory, change this
# variable to the correct path.
##############################################################################

SYMPHONYROOT=..\..

##############################################################################
# COINROOT is the path to the root directory of the COIN libraries. Many of
# the new features of COIN require the COIN libraries to be installed.
##############################################################################

COINROOT = ..\$(SYMPHONYROOT)

##############################################################################
# OUTDIR variable specifies where to create the executable file, 
# "symphony.exe", the corresponding objects and the dependencies.  
##############################################################################

OUTDIR=.\Debug

##############################################################################
##############################################################################
# LP solver dependent definitions
##############################################################################
##############################################################################

##############################################################################
##############################################################################
#You must define an LP solver in order to use the software. By default, this 
# option is set to OsI_CLP. See the corresponding "LPINCDIR" and "LPLIB" 
# variables used to put the lp solver include files and the libraries on path
# and make the necessary changes if you require.
##############################################################################
##############################################################################

##############################################################################
# CPLEX definitions
##############################################################################

# Uncomment the line below if you want to use CPLEX and specify the 
# corresponding paths to the solver files and libraries. 

#LP_SOLVER = CPLEX
!IF "$(LP_SOLVER)" == "CPLEX"
LPINCDIR = "C:\ILOG\cplex81\include\ilcplex"
LPLIB = "C:\ILOG\cplex81\lib\msvc6\stat_sta\cplex81.lib"
!ENDIF

##############################################################################
# OSL definitions
##############################################################################

# Uncomment the line below if you want to use OSL and specify the 
# corresponding paths to the solver files and libraries. 

#LP_SOLVER = OSL
!IF "$(LP_SOLVER)" == "OSL"
LPINCDIR = "C:\Program Files\IbmOslV3Lib\osllib\include"
LPLIB = "C:\Program Files\IbmOslV3Lib\osllib\lib\oslmd6030.lib"
!ENDIF

##############################################################################
# OSI definitions
##############################################################################

# Uncomment the line below if you want to use OSI interface and specify the 
# corresponding paths to the solver files and libraries. 

LP_SOLVER = OSI
OSI_INTERFACE = CLP

!IF "$(LP_SOLVER)" == "OSI"
LPINCDIR = \
	"$(COINROOT)\Coin\include" /I\
	"$(COINROOT)\Osi\include"
LPLIB = \
	"$(COINROOT)\Win\coinLib\Debug\coinLib.lib" \
	"$(COINROOT)\Win\osiLib\Debug\osiLib.lib"
!ENDIF


!IF "$(OSI_INTERFACE)" == "CPLEX"
LPINCDIR = $(LPINCDIR) /I\
	"C:\ILOG\cplex81\include\ilcplex" /I\
	"$(COINROOT)\Osi\OsiCpx\include"
LPLIB = $(LPLIB) \
	"C:\ILOG\cplex81\lib\msvc6\stat_sta\cplex81.lib" \
	"$(COINROOT)\Win\osiCpxLib\Debug\osiCpxLib.lib"
!ENDIF


!IF "$(OSI_INTERFACE)" == "OSL"
LPINCDIR = $(LPINCDIR) /I\
	"C:\Program Files\IbmOslV3Lib\osllib\include" /I\
	"$(COINROOT)\Osi\OsiOsl\include"
LPLIB = $(LPLIB) \
	"C:\Program Files\IbmOslV3Lib\osllib\lib\oslmd6030.lib" \
        "$(COINROOT)\Win\osiOslLib\Debug\osiOslLib.lib"
!ENDIF


!IF "$(OSI_INTERFACE)" == "CLP"
LPINCDIR = $(LPINCDIR) /I\
	"$(COINROOT)\Clp\include" /I\
	"$(COINROOT)\Osi\OsiClp\include"
LPLIB = $(LPLIB) \
	"$(COINROOT)\Win\clpLib\Debug\clpLib.lib" \
	"$(COINROOT)\Win\osiClpLib\Debug\osiClpLib.lib"
!ENDIF

!IF "$(OSI_INTERFACE)" == "XPRESS"
LPINCDIR = $(LPINCDIR) /I\
	"C:\" /I\
	"$(COINROOT)\Osi\OsiXpr\include"
LPLIB = $(LPLIB) \
	"C:\" \
	"$(COINROOT)\Win\osiXprLib\Debug\osiXprLib.lib"
!ENDIF

!IF "$(OSI_INTERFACE)" == "SOPLEX"
LPINCDIR = $(LPINCDIR) /I\
	"C:\" /I\
	"$(COINROOT)\Osi\OsiSpx\include"
LPLIB = $(LPLIB) \
	"C:\" \
	"$(COINROOT)\Win\osiSpxLib\Debug\osiSpxLib.lib"
!ENDIF

!IF "$(OSI_INTERFACE)" == "DYLP"
LPINCDIR = $(LPINCDIR) /I\
	"C:\" /I\
	"$(COINROOT)\Osi\OsiDylp\include"
LPLIB = $(LPLIB) \
	"C:\" \
	"$(COINROOT)\Win\osiDylpLib\Debug\osiDylpLib.lib"
!ENDIF


!IF "$(OSI_INTERFACE)" == "GLPK"
LPINCDIR = $(LPINCDIR) /I\
	"C:\GLPK\glpk-4.0\include" /I\
	"$(COINROOT)\Osi\OsiGlpk\include"
LPLIB = $(LPLIB) \
	"C:\GLPK\glpk-4.0\glpk.lib" \
	"$(COINROOT)\Win\osiGlpkLib\Debug\osiGlpkLib.lib"
!ENDIF


##############################################################################
# Besides the above variables, you have to set your environment path to solver
# specific dynamic libraries if there exists any. For instance, you have to 
# set your path to where "cplex81.dll" is, something like: 
#
#              "set path = %path%;C:\ILOG\cplex81\bin\msvc6 " 
#
# if you are using CPLEX 8.1 and Visual C++ 6. 
##############################################################################

##############################################################################
# SOLVER definition for SYMPHONY
##############################################################################

!IF "$(LP_SOLVER)" == "OSI"
DEFINITIONS ="__$(LP_SOLVER)_$(OSI_INTERFACE)__"
!ELSE
DEFINITIONS ="__$(LP_SOLVER)__"
!ENDIF

##############################################################################
# GLPMPL definitions. The user should set "USE_GLPMPL" variable to "TRUE" and 
# specify the paths for "glpk" files if she wants to read in glpmpl files.  
##############################################################################

USE_GLPMPL = TRUE

!IF "$(USE_GLPMPL)" == "TRUE"
DEFINITIONS = $(DEFINITIONS) /D "USE_GLPMPL"
!ENDIF

##############################################################################
##############################################################################
# Generate generic cutting planes. If you are using the OSI interface, you 
# can now add generic cutting planes from the CGL by setting the flag below.
# Which cutting planes are added can be controlled by SYMPHONY parameters (see
# the user's manual
##############################################################################
##############################################################################

USE_CGL_CUTS = TRUE

!IF "$(USE_CGL_CUTS)" == "TRUE"
LPINCDIR = $(LPINCDIR) /I "$(COINROOT)\Cgl\include"
LPLIB = $(LPLIB) "$(COINROOT)\Win\cglLib\Debug\cglLib.lib"
DEFINITIONS= $(DEFINITIONS) /D "USE_CGL_CUTS"
!ENDIF

##############################################################################
# If you wish to compile and use the SYMPHONY callable library through the 
# SYMPHONY OSI interface, set USE_OSI_INTERFACE to TRUE below. Note that
# you must have COIN installed to use this capability. See below to set the 
# path to the COIN directories. 
##############################################################################

USE_OSI_INTERFACE = FALSE

!IF "$(USE_OSI_INTERFACE)" == "TRUE"
ALL_INCDIR = $(LPINCDIR) /I "$(COINROOT)\Osi\OsiSym\include"
ALL_LIB = $(LPLIB) "$(COINROOT)\Win\osiSymLib\Debug\osiSymLib.lib"
!ELSE
ALL_INCDIR = $(LPINCDIR)
ALL_LIB = $(LPLIB)
!ENDIF


##############################################################################
##############################################################################
#
# Compiling and Linking...
#
##############################################################################
##############################################################################

DEFINITIONS = $(DEFINITIONS) /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" \
	/D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" \
	/D "COMPILE_IN_TM"

GMPL_INCDIR = "$(SYMPHONYROOT)\src\GMPL"
ALL_INCDIR =$(ALL_INCDIR) /I $(GMPL_INCDIR) /I "$(SYMPHONYROOT)\include"

.SILENT:

CPP=cl.exe
CPPFLAGS= /nologo /MLd /W2 /GR /Gm /YX /GX /ZI /Od \
	/I $(ALL_INCDIR) /D $(DEFINITIONS) \
	/Fp"$(OUTDIR)\sym.pch" /Fo"$(OUTDIR)\\" /Fd"$(OUTDIR)\sym.pdb" \
	/FD /GZ /c /Tp

CFLAGS= /nologo /MLd /W2 /GR /Gm /YX /GX /ZI /Od \
	/I $(GMPL_INCDIR) /D $(DEFINITIONS) \
	/Fp"$(OUTDIR)\sym.pch" /Fo"$(OUTDIR)\\" /Fd"$(OUTDIR)\sym.pdb" \
	/FD /GZ /c
.c.obj: 
	$(CPP) $(CPPFLAGS) "$*.c"
	 
.c.cobj: 
	$(CPP) $(CFLAGS) "$*.c"

ALL : "$(OUTDIR)" "LIB_MESSAGE" sym_lib "SYMPHONY_MESSAGE" "OBJECTS" sym_exe 

CLEAN:
	del /Q $(OUTDIR)\*.obj
	del /Q $(OUTDIR)\symphony.exe 
	del /Q $(OUTDIR)\symphonyLib.lib
	del /Q $(OUTDIR)\sym.idb
	del /Q $(OUTDIR)\sym.pdb
	del /Q $(OUTDIR)\sym.pch

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

LIB_MESSAGE:
	echo Creating SYMPHONY library...	

SYMPHONY_MESSAGE:
	echo Compiling SYMPHONY main function...

sym_lib : \
	$(SYMPHONYROOT)\src\Common\pack_array.obj \
	$(SYMPHONYROOT)\src\Common\pack_cut.obj \
	$(SYMPHONYROOT)\src\Common\proccomm.obj \
	$(SYMPHONYROOT)\src\Common\qsortucb.obj \
	$(SYMPHONYROOT)\src\Common\qsortucb_di.obj \
	$(SYMPHONYROOT)\src\Common\qsortucb_i.obj \
	$(SYMPHONYROOT)\src\Common\qsortucb_ic.obj \
	$(SYMPHONYROOT)\src\Common\qsortucb_id.obj \
	$(SYMPHONYROOT)\src\Common\qsortucb_ii.obj \
	$(SYMPHONYROOT)\src\Common\timemeas.obj \
	$(SYMPHONYROOT)\src\CutGen\cg_func.obj \
	$(SYMPHONYROOT)\src\CutGen\cg_proccomm.obj \
	$(SYMPHONYROOT)\src\CutGen\cg_wrapper.obj \
	$(SYMPHONYROOT)\src\CutGen\cut_gen.obj \
	$(SYMPHONYROOT)\src\CutPool\cp_func.obj \
	$(SYMPHONYROOT)\src\CutPool\cp_proccomm.obj \
	$(SYMPHONYROOT)\src\CutPool\cp_wrapper.obj \
	$(SYMPHONYROOT)\src\CutPool\cut_pool.obj \
	$(SYMPHONYROOT)\src\LP\lp.obj \
	$(SYMPHONYROOT)\src\LP\lp_branch.obj \
	$(SYMPHONYROOT)\src\LP\lp_free.obj \
	$(SYMPHONYROOT)\src\LP\lp_genfunc.obj \
	$(SYMPHONYROOT)\src\LP\lp_proccomm.obj \
	$(SYMPHONYROOT)\src\LP\lp_rowfunc.obj \
	$(SYMPHONYROOT)\src\LP\lp_solver.obj \
	$(SYMPHONYROOT)\src\LP\lp_varfunc.obj \
	$(SYMPHONYROOT)\src\LP\lp_wrapper.obj \
	$(SYMPHONYROOT)\src\Master\master.obj \
	$(SYMPHONYROOT)\src\Master\master_func.obj \
	$(SYMPHONYROOT)\src\Master\master_io.obj \
	$(SYMPHONYROOT)\src\Master\master_wrapper.obj \
	$(SYMPHONYROOT)\src\TreeManager\tm_func.obj \
	$(SYMPHONYROOT)\src\TreeManager\tm_proccomm.obj \
	$(SYMPHONYROOT)\src\TreeManager\treemanager.obj \
	$(SYMPHONYROOT)\src\GMPL\glpavl.cobj \
	$(SYMPHONYROOT)\src\GMPL\glpdmp.cobj \
	$(SYMPHONYROOT)\src\GMPL\glplib1a.cobj \
	$(SYMPHONYROOT)\src\GMPL\glplib2.cobj \
	$(SYMPHONYROOT)\src\GMPL\glplib3.cobj \
	$(SYMPHONYROOT)\src\GMPL\glpmpl1.cobj \
	$(SYMPHONYROOT)\src\GMPL\glpmpl2.cobj \
	$(SYMPHONYROOT)\src\GMPL\glpmpl3.cobj \
	$(SYMPHONYROOT)\src\GMPL\glpmpl4.cobj \
	$(SYMPHONYROOT)\src\GMPL\glprng.cobj
	lib.exe /nologo /out:$(OUTDIR)\symphonyLib.lib $(OUTDIR)\*.obj
	echo "symphonyLib.lib" created successfully...
	echo ...

LINK_OBJECTS= \
	$(OUTDIR)\main.obj

OBJECTS : \
	$(SYMPHONYROOT)\src\Master\main.obj
	echo main compiled successfully...
	echo ...	
               	          
sym_exe : $(LINK_OBJECTS) $(OUTDIR)\symphonyLib.lib
	echo Linking...
	$(CPP) /nologo /W3 /Fe"$(OUTDIR)\symphony.exe" $(ALL_LIB) $**
	echo "symphony.exe" created successfully...
