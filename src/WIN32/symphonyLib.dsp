# Microsoft Developer Studio Project File - Name="symphonyLib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=symphonyLib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "symphonyLib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "symphonyLib.mak" CFG="symphonyLib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "symphonyLib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "symphonyLib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "symphonyLib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "symphonyLib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "symphonyLib___Win32_Debug"
# PROP BASE Intermediate_Dir "symphonyLib___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W2 /Gm /GR /GX /ZI /Od /I "C:\COIN\Osi\include" /I "C:\COIN\Osi\OsiClp\include" /I "C:\COIN\Clp\include" /I "C:\COIN\Coin\include" /I "C:\COIN\Cgl\include" /I "..\include" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "INTEL" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "__OSI_CLP__" /D "USE_CGL_CUTS" /FD /GZ /c /Tp
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "symphonyLib - Win32 Release"
# Name "symphonyLib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "Common"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\Common\pack_array.c
# End Source File
# Begin Source File

SOURCE=..\Common\pack_cut.c
# End Source File
# Begin Source File

SOURCE=..\Common\proccomm.c
# End Source File
# Begin Source File

SOURCE=..\Common\qsortucb.c
# End Source File
# Begin Source File

SOURCE=..\Common\qsortucb_di.c
# End Source File
# Begin Source File

SOURCE=..\Common\qsortucb_i.c
# End Source File
# Begin Source File

SOURCE=..\Common\qsortucb_ic.c
# End Source File
# Begin Source File

SOURCE=..\Common\qsortucb_id.c
# End Source File
# Begin Source File

SOURCE=..\Common\qsortucb_ii.c
# End Source File
# Begin Source File

SOURCE=..\Common\timemeas.c
# End Source File
# End Group
# Begin Group "CutGen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\CutGen\cg_func.c
# End Source File
# Begin Source File

SOURCE=..\CutGen\cg_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\CutGen\cg_wrapper.c
# End Source File
# Begin Source File

SOURCE=..\CutGen\cut_gen.c
# End Source File
# End Group
# Begin Group "CutPool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\CutPool\cp_func.c
# End Source File
# Begin Source File

SOURCE=..\CutPool\cp_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\CutPool\cp_wrapper.c
# End Source File
# Begin Source File

SOURCE=..\CutPool\cut_pool.c
# End Source File
# End Group
# Begin Group "LP"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\LP\lp.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_branch.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_free.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_genfunc.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_rowfunc.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_solver.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_varfunc.c
# End Source File
# Begin Source File

SOURCE=..\LP\lp_wrapper.c
# End Source File
# End Group
# Begin Group "Master"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\Master\master.c
# End Source File
# Begin Source File

SOURCE=..\Master\master_func.c
# End Source File
# Begin Source File

SOURCE=..\Master\master_io.c
# End Source File
# Begin Source File

SOURCE=..\Master\master_wrapper.c
# End Source File
# End Group
# Begin Group "TreeManager"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\TreeManager\tm_func.c
# End Source File
# Begin Source File

SOURCE=..\TreeManager\tm_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\TreeManager\treemanager.c
# End Source File
# End Group
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
