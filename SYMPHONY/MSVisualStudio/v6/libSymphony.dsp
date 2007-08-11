# Microsoft Developer Studio Project File - Name="libSymphony" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libSymphony - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libSymphony.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libSymphony.mak" CFG="libSymphony - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libSymphony - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "libSymphony - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libSymphony - Win32 Release"

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
# ADD CPP /nologo /MLd /W2 /Gm /GR /GX /ZI /Od /I "..\..\..\Osi\src" /I "..\..\..\Osi\src\OsiClp" /I "..\..\..\Clp\src" /I "..\..\..\CoinUtils\src" /I "..\..\..\Cgl\src" /I "..\..\..\Cgl\src\CglLiftAndProject" /I "..\..\..\Cgl\src\CglLandP" /I "..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\Cgl\src\CglClique" /I "..\..\..\Cgl\src\CglOddHole" /I "..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\Cgl\src\CglGomory" /I "..\..\..\Cgl\src\CglSimpleRounding" /I "..\..\..\Cgl\src\CglProbing" /I "..\..\..\BuildTools\headers" /I "..\..\include" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "INTEL" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "__OSI_CLP__" /D "USE_CGL_CUTS" /Fo"Debug/" /Fd"Debug/" /FD /GZ /c /Tp
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libSymphony - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W2 /Gm /GR /GX /ZI /Od /I "..\..\..\Osi\src" /I "..\..\..\Osi\src\OsiClp" /I "..\..\..\Clp\src" /I "..\..\..\CoinUtils\src" /I "..\..\..\Cgl\src" /I "..\..\..\Cgl\src\CglLiftAndProject" /I "..\..\..\Cgl\src\CglLandP" /I "..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\Cgl\src\CglClique" /I "..\..\..\Cgl\src\CglOddHole" /I "..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\Cgl\src\CglGomory" /I "..\..\..\Cgl\src\CglSimpleRounding" /I "..\..\..\Cgl\src\CglTwomir" /I "..\..\..\Cgl\src\CglProbing" /I "..\..\..\Cgl\src\CglRedSplit" /I "..\..\..\BuildTools\headers" /I "..\..\include" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "INTEL" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "__OSI_CLP__" /D "USE_CGL_CUTS" /FD /GZ /c /Tp
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

# Name "libSymphony - Win32 Release"
# Name "libSymphony - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "Common"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\Common\pack_array.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\pack_cut.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\qsortucb.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\qsortucb_di.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\qsortucb_i.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\qsortucb_ic.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\qsortucb_id.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\qsortucb_ii.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\timemeas.c
# End Source File
# End Group
# Begin Group "CutGen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\CutGen\cg_func.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutGen\cg_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutGen\cg_wrapper.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutGen\cut_gen.c
# End Source File
# End Group
# Begin Group "CutPool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\CutPool\cp_func.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutPool\cp_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutPool\cp_wrapper.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutPool\cut_pool.c
# End Source File
# End Group
# Begin Group "LP"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\LP\lp.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_branch.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_free.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_genfunc.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_rowfunc.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_solver.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_varfunc.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\lp_wrapper.c
# End Source File
# End Group
# Begin Group "Master"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\Master\master.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\master_func.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\master_io.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\master_wrapper.c
# End Source File
# End Group
# Begin Group "TreeManager"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\TreeManager\tm_func.c
# End Source File
# Begin Source File

SOURCE=..\..\src\TreeManager\tm_proccomm.c
# End Source File
# End Group
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
