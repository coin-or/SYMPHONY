# Microsoft Developer Studio Project File - Name="vrp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=vrp - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "vrp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "vrp.mak" CFG="vrp - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "vrp - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "vrp - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "vrp - Win32 Release"

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

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

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
# ADD CPP /nologo /W2 /Gm /GR /GX /ZI /Od /I "..\..\..\..\include" /I "..\..\include" /I "..\..\..\..\..\Osi\src" /I "..\..\..\..\..\CoinUtils\src" /I "..\..\..\..\..\Clp\src" /I "..\..\..\..\..\Cgl\src" /I "..\..\..\..\..\Cgl\src\CglLiftAndProject" /I "..\..\..\..\..\Cgl\src\CglGomory" /I "..\..\..\..\..\Cgl\src\CglClique" /I "..\..\..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\..\..\Cgl\src\CglProbing" /I "..\..\..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\..\..\Cgl\src\CglOddHole" /I "..\..\..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\..\..\Cgl\src\CglSimpleRounding" /I "..\..\..\..\..\Osi\src\OsiClp" /I "..\..\..\..\..\BuildTools\headers" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "__OSI_CLP__" /Fp"Debug/vrp2.pch" /YX /FD /GZ /c /Tp
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

# Name "vrp - Win32 Release"
# Name "vrp - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter ""
# Begin Group "Common"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\Common\compute_cost.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\network.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Common\vrp_macros.c
# End Source File
# End Group
# Begin Group "CutGen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\CutGen\biconnected.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutGen\shrink.c
# End Source File
# Begin Source File

SOURCE=..\..\src\CutGen\tsp.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\CutGen\vrp_cg.c
# End Source File
# End Group
# Begin Group "CutPool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\CutPool\vrp_cp.c
# End Source File
# End Group
# Begin Group "DrawGraph"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\DrawGraph\vrp_dg_functions.c
# End Source File
# End Group
# Begin Group "LP"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\LP\vrp_lp.c
# End Source File
# Begin Source File

SOURCE=..\..\src\LP\vrp_lp_branch.c
# End Source File
# End Group
# Begin Group "Master"

# PROP Default_Filter ""
# Begin Group "Heuristics"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\src\Master\Heuristics\cluster_heur.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\Master\Heuristics\exchange_heur.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\Master\Heuristics\lower_bound.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\Master\Heuristics\receive_rout.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\Master\Heuristics\route_heur.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\Master\Heuristics\start_heurs.c

!IF  "$(CFG)" == "vrp - Win32 Release"

!ELSEIF  "$(CFG)" == "vrp - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# End Group
# Begin Source File

SOURCE=..\..\src\Master\small_graph.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\vrp_io.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\vrp_main.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\vrp_master.c
# End Source File
# Begin Source File

SOURCE=..\..\src\Master\vrp_master_functions.c
# End Source File
# End Group
# End Group
# End Target
# End Project
