##############################################################################
##############################################################################
#                                                                            #
# This file is part of the SYMPHONY Branch, Cut, and Price Library.          #
#                                                                            #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and     #
# Laci Ladanyi (ladanyi@us.ibm.com).                                         #
#                                                                            #
# (c) Copyright 2000, 2001, 2002  Ted Ralphs. All Rights Reserved.           #
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
# nmake /f mpp.mak 
#
# The executable "symphony.exe" for this application will be created in 
# .\Debug directory. By default, SYMPHONY is set up to use the CPLEX 8.1
# optimization solver via COIN-OSI's CPLEX interface and the COIN-CGL cuts. 
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

SYMPHONYROOT=..\..\

##############################################################################
# COINROOT is the path to the root directory of the COIN libraries. Many of
# the new features of COIN require the COIN libraries to be installed.
##############################################################################

COINROOT = C:\COIN

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
# option is set to OsI_CPLEX. See the corresponding "LPINCDIR" and "LPLIB" 
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
OSI_INTERFACE = CPLEX

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

USE_GLPMPL = FALSE

!IF "$(USE_GLPMPL)" == "TRUE"
LPINCDIR = $LPINCDIR) C:\GLPK\glpk-4.0\include
LPLIB = $(LPLIB) C:\GLPK\glpk-4.0\glpk.lib
DEFINITIONs = $(DEFINITIONS) /D "USE_GLPMPL"
!ENDIF

##############################################################################
##############################################################################
# Generate generic cutting planes. If you are using the OSI interface, you 
# can now add generic cutting planes from the CGL by setting the flag below.
# Which cutting planes are added can be controlled by SYMPHONY parameters (see
# the user's manual
##############################################################################
##############################################################################

USE_CGL_CUTS = FALSE

!IF "$(USE_CGL_CUTS)" == "TRUE"
LPINCDIR = $(LPINCDIR) /I "$(COINROOT)\Cgl\include"
LPLIB = $(LPLIB) "$(COINROOT)\Win\cglLib\Debug\cglLib.lib"
DEFINITIONS= $(DEFINITIONS) /D "USE_CGL_CUTS"
!ENDIF

##############################################################################
##############################################################################
#
# Compiling and Linking...
#
##############################################################################
##############################################################################

.SILENT:

CPP=cl.exe
CPPFLAGS= /nologo /MLd /W3 /GR /Gm /YX /GX /ZI /Od \
	/I $(LPINCDIR) /I "$(SYMPHONYROOT)\include" /I "..\include" \
	/D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" \
	/D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" \
	/D "COMPILE_IN_TM" /D $(DEFINITIONS) \
	/Fp"$(OUTDIR)\mpp.pch" /Fo"$(OUTDIR)\\" /Fd"$(OUTDIR)\mpp.pdb" \
	/FD /GZ /c /Tp

.c.obj: 
	$(CPP) $(CPPFLAGS) "$*.c"
	 
ALL : "APPL_MESSAGE" mpp.lib "SYMPHONY_MESSAGE" "OBJECTS" symphony.exe 

CLEAN:
	del /Q $(OUTDIR)\*.obj
        del /Q $(OUTDIR)\symphony.exe 
        del /Q $(OUTDIR)\mpp.lib
        del /Q $(OUTDIR)\mpp.idb
        del /Q $(OUTDIR)\mpp.pdb
	del /Q $(OUTDIR)\mpp.pch

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"


APPL_MESSAGE: 
	echo Compiling application files...

SYMPHONY_MESSAGE:
	echo Compiling SYMPHONY files...	

mpp.lib : \
	..\CutGen\mpp_cg.obj \
	..\CutPool\mpp_cp.obj \
	..\DrawGraph\mpp_dg.obj \
	..\LP\mpp_lp.obj \
	..\LP\mpp_lp_branch.obj \
	..\Master\mpp_master.obj
	lib.exe /nologo /out:$(OUTDIR)\mpp.lib $(OUTDIR)\*.obj
	echo Application files compiled successfully...
	echo ...

LINK_OBJS= \
	$(OUTDIR)\pack_array.obj \
	$(OUTDIR)\pack_cut.obj \
	$(OUTDIR)\proccomm.obj \
	$(OUTDIR)\qsortucb.obj \
	$(OUTDIR)\qsortucb_di.obj \
	$(OUTDIR)\qsortucb_i.obj \
	$(OUTDIR)\qsortucb_ic.obj \
	$(OUTDIR)\qsortucb_id.obj \
	$(OUTDIR)\qsortucb_ii.obj \
	$(OUTDIR)\timemeas.obj \
	$(OUTDIR)\cg_func.obj \
	$(OUTDIR)\cg_proccomm.obj \
	$(OUTDIR)\cg_wrapper.obj \
	$(OUTDIR)\cut_gen.obj \
	$(OUTDIR)\cp_func.obj \
	$(OUTDIR)\cp_proccomm.obj \
	$(OUTDIR)\cp_wrapper.obj \
	$(OUTDIR)\cut_pool.obj \
	$(OUTDIR)\lp.obj \
	$(OUTDIR)\lp_branch.obj \
	$(OUTDIR)\lp_free.obj \
	$(OUTDIR)\lp_genfunc.obj \
	$(OUTDIR)\lp_proccomm.obj \
	$(OUTDIR)\lp_rowfunc.obj \
	$(OUTDIR)\lp_solver.obj \
	$(OUTDIR)\lp_varfunc.obj \
	$(OUTDIR)\lp_wrapper.obj \
	$(OUTDIR)\master.obj \
	$(OUTDIR)\master_io.obj \
	$(OUTDIR)\master_wrapper.obj \
	$(OUTDIR)\tm_func.obj \
	$(OUTDIR)\tm_proccomm.obj \
	$(OUTDIR)\treemanager.obj 


OBJECTS : \
	$(SYMPHONYROOT)\Common\pack_array.obj \
	$(SYMPHONYROOT)\Common\pack_cut.obj \
	$(SYMPHONYROOT)\Common\proccomm.obj \
	$(SYMPHONYROOT)\Common\qsortucb.obj \
	$(SYMPHONYROOT)\Common\qsortucb_di.obj \
	$(SYMPHONYROOT)\Common\qsortucb_i.obj \
	$(SYMPHONYROOT)\Common\qsortucb_ic.obj \
	$(SYMPHONYROOT)\Common\qsortucb_id.obj \
	$(SYMPHONYROOT)\Common\qsortucb_ii.obj \
	$(SYMPHONYROOT)\Common\timemeas.obj \
	$(SYMPHONYROOT)\CutGen\cg_func.obj \
	$(SYMPHONYROOT)\CutGen\cg_proccomm.obj \
	$(SYMPHONYROOT)\CutGen\cg_wrapper.obj \
	$(SYMPHONYROOT)\CutGen\cut_gen.obj \
	$(SYMPHONYROOT)\CutPool\cp_func.obj \
	$(SYMPHONYROOT)\CutPool\cp_proccomm.obj \
	$(SYMPHONYROOT)\CutPool\cp_wrapper.obj \
	$(SYMPHONYROOT)\CutPool\cut_pool.obj \
	$(SYMPHONYROOT)\LP\lp.obj \
	$(SYMPHONYROOT)\LP\lp_branch.obj \
	$(SYMPHONYROOT)\LP\lp_free.obj \
	$(SYMPHONYROOT)\LP\lp_genfunc.obj \
	$(SYMPHONYROOT)\LP\lp_proccomm.obj \
	$(SYMPHONYROOT)\LP\lp_rowfunc.obj \
	$(SYMPHONYROOT)\LP\lp_solver.obj \
	$(SYMPHONYROOT)\LP\lp_varfunc.obj \
	$(SYMPHONYROOT)\LP\lp_wrapper.obj \
	$(SYMPHONYROOT)\Master\master.obj \
	$(SYMPHONYROOT)\Master\master_io.obj \
	$(SYMPHONYROOT)\Master\master_wrapper.obj \
	$(SYMPHONYROOT)\TreeManager\tm_func.obj \
	$(SYMPHONYROOT)\TreeManager\tm_proccomm.obj \
	$(SYMPHONYROOT)\TreeManager\treemanager.obj
	echo SYMPHONY files compiled successfully...	
               	          
symphony.exe : $(LINK_OBJS) $(OUTDIR)\mpp.lib
	echo Linking...
	$(CPP) /nologo /W3 /Fe"$(OUTDIR)\symphony.exe" $(LPLIB) $**
	echo "symphony.exe" created successfully...
