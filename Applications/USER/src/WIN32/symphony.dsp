# Microsoft Developer Studio Project File - Name="symphony" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=symphony - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "symphony.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "symphony.mak" CFG="symphony - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "symphony - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "symphony - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "symphony - Win32 Release"

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

!ELSEIF  "$(CFG)" == "symphony - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "symphony___Win32_Debug"
# PROP BASE Intermediate_Dir "symphony___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "C:\COIN\Osi\include" /I "C:\COIN\Osi\OsiCpx\include" /I "C:\COIN\Coin\include" /I "C:\COIN\Cgl\include" /I "..\..\include" /I "C:\ILOG\cplex81\include\ilcplex" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "INTEL" /D "COMPILE_IN_CG" /D "COMPILE_IN_CP" /D "COMPILE_IN_LP" /D "COMPILE_IN_TM" /D "__OSI_CPLEX__" /D "USE_CGL_CUTS" /FD /GZ /c /Tp
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

# Name "symphony - Win32 Release"
# Name "symphony - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "Common"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Common\pack_array.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\pack_cut.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\qsortucb.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\qsortucb_di.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\qsortucb_i.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\qsortucb_ic.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\qsortucb_id.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\qsortucb_ii.c
# End Source File
# Begin Source File

SOURCE=..\..\Common\timemeas.c
# End Source File
# End Group
# Begin Group "CutGen"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\CutGen\cg_func.c
# End Source File
# Begin Source File

SOURCE=..\..\CutGen\cg_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\CutGen\cg_wrapper.c
# End Source File
# Begin Source File

SOURCE=..\..\CutGen\cut_gen.c
# End Source File
# End Group
# Begin Group "CutPool"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\CutPool\cp_func.c
# End Source File
# Begin Source File

SOURCE=..\..\CutPool\cp_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\CutPool\cp_wrapper.c
# End Source File
# Begin Source File

SOURCE=..\..\CutPool\cut_pool.c
# End Source File
# End Group
# Begin Group "LP"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\LP\lp.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_branch.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_free.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_genfunc.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_pseudo_branch.c

!IF  "$(CFG)" == "symphony - Win32 Release"

!ELSEIF  "$(CFG)" == "symphony - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_rowfunc.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_solver.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_varfunc.c
# End Source File
# Begin Source File

SOURCE=..\..\LP\lp_wrapper.c
# End Source File
# End Group
# Begin Group "Master"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Master\master.c
# End Source File
# Begin Source File

SOURCE=..\..\Master\master_io.c
# End Source File
# Begin Source File

SOURCE=..\..\Master\master_wrapper.c
# End Source File
# End Group
# Begin Group "TreeManager"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\TreeManager\tm_func.c
# End Source File
# Begin Source File

SOURCE=..\..\TreeManager\tm_proccomm.c
# End Source File
# Begin Source File

SOURCE=..\..\TreeManager\treemanager.c
# End Source File
# End Group
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\Debug\user.lib
# End Source File
# Begin Source File

SOURCE=C:\COIN\Win\cglLib\Debug\cglLib.lib
# End Source File
# Begin Source File

SOURCE=C:\COIN\Win\coinLib\Debug\coinLib.lib
# End Source File
# Begin Source File

SOURCE=C:\COIN\Win\osiLib\Debug\osiLib.lib
# End Source File
# Begin Source File

SOURCE=C:\COIN\Win\osiCpxLib\Debug\osiCpxLib.lib
# End Source File
# Begin Source File

SOURCE=C:\ILOG\cplex81\lib\msvc6\stat_sta\cplex81.lib
# End Source File
# End Target
# End Project
