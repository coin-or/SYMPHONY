# Microsoft Developer Studio Project File - Name="cnrp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=cnrp - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cnrp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cnrp.mak" CFG="cnrp - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cnrp - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "cnrp - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "cnrp - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "cnrp - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W2 /Gm /GR /GX /ZI /Od /I "..\..\include" /I "..\include" /I "C:\COIN\Osi\include" /I "C:\Coin\Coin\include" /I "C:\Coin\Clp\include" /I "C:\Coin\Osi\OsiClp\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "DIRECTED_X_VARS" /D "ADD_FLOW_VARS" /D "SAVE_CUT_POOL" /D "MULTI_CRITERIA" /D "__OSI_CLP__" /Fp"Debug/cnrp2.pch" /YX /FD /GZ /c /Tp
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "cnrp - Win32 Release"
# Name "cnrp - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter ""
# Begin Group "Common"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\Common\cnrp_macros.c
# End Source File
# Begin Source File

SOURCE=..\Common\compute_cost.c
# End Source File
# Begin Source File

SOURCE=..\Common\network.c
# End Source File
# End Group
# Begin Group "TreeManager"

# PROP Default_Filter ""
# End Group
# Begin Group "LP"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\LP\cnrp_lp.c
# End Source File
# Begin Source File

SOURCE=..\LP\cnrp_lp_branch.c
# End Source File
# End Group
# Begin Group "CutPool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\CutPool\cnrp_cp.c
# End Source File
# End Group
# Begin Group "CutGen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\CutGen\biconnected.c
# End Source File
# Begin Source File

SOURCE=..\CutGen\cnrp_cg.c
# End Source File
# Begin Source File

SOURCE=..\CutGen\shrink.c
# End Source File
# End Group
# Begin Group "Master"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\Master\cnrp_io.c
# End Source File
# Begin Source File

SOURCE=..\Master\cnrp_main.c
# End Source File
# Begin Source File

SOURCE=..\Master\cnrp_master.c
# End Source File
# Begin Source File

SOURCE=..\Master\cnrp_master_functions.c
# End Source File
# Begin Source File

SOURCE=..\Master\small_graph.c
# End Source File
# End Group
# Begin Group "DrawGraph"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\DrawGraph\cnrp_dg_functions.c
# End Source File
# End Group
# End Group
# Begin Source File

SOURCE=..\..\..\COIN\Win\clpLib\Debug\clpLib.lib
# End Source File
# Begin Source File

SOURCE=..\..\..\COIN\Win\coinLib\Debug\coinLib.lib
# End Source File
# Begin Source File

SOURCE=..\..\..\COIN\Win\cglLib\Debug\cglLib.lib
# End Source File
# Begin Source File

SOURCE=..\..\..\COIN\Win\osiLib\Debug\osiLib.lib
# End Source File
# Begin Source File

SOURCE=..\..\..\COIN\Win\osiClpLib\Debug\osiClpLib.lib
# End Source File
# Begin Source File

SOURCE=.\Debug\symphonyLib.lib
# End Source File
# End Target
# End Project
