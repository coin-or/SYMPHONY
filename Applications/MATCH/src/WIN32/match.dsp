# Microsoft Developer Studio Project File - Name="match" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=match - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "match.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "match.mak" CFG="match - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "match - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "match - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "match - Win32 Release"

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

!ELSEIF  "$(CFG)" == "match - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "match___Win32_Debug"
# PROP BASE Intermediate_Dir "match___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W2 /Gm /GR /GX /ZI /Od /I "..\..\include" /I "..\include" /I "C:\COIN\Osi\include" /I "C:\Coin\Coin\include" /I "C:\COIN\Cgl\include" /I "C:\Coin\Clp\include" /I "C:\Coin\Osi\OsiClp\include" /I "C:\Coin\Osi\OsiSym\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "__OSI_CLP__" /D "USE_CGL_CUTS" /Fp"Debug/vrp.pch" /YX /FD /GZ /c /Tp
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

# Name "match - Win32 Release"
# Name "match - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "CutGen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\CutGen\user_cg.c
# End Source File
# End Group
# Begin Group "CutPool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\CutPool\user_cp.c
# End Source File
# End Group
# Begin Group "DrawGraph"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\DrawGraph\user_dg.c
# End Source File
# End Group
# Begin Group "LP"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\LP\user_lp.c
# End Source File
# Begin Source File

SOURCE=..\LP\user_lp_branch.c
# End Source File
# End Group
# Begin Group "Master"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\Master\user_main.c
# End Source File
# Begin Source File

SOURCE=..\Master\user_master.c
# End Source File
# End Group
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# Begin Group "Common"

# PROP Default_Filter ""
# End Group
# End Target
# End Project
