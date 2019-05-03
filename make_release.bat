@echo off

REM
REM This script loads the required modules, executes cmake, make and make install.
REM ------------------------------------------------------------------------------

set CMAKE_FIND_SCRIPTS=T:/modulesystem/tools/cmake/cmake-find-scripts
set VISTA_CMAKE_COMMON=T:\modulesystem\tools\vista\vista-dlr-develop\src\VistaCMakeCommon

SET ParaView_DIR=D:\SoftwareEntwicklung\Catalyst\install\win7.x86_64.msvc15.release\lib\cmake\paraview-5.6

REM create some required variables -----------------------------------------------

REM this directory should contain the top-level CMakeLists.txt - it is assumed to
REM reside in the same directory as this script. If the script is called from
REM elsewhere it does not matter

set CMAKE_DIR=%~dp0

REM get the current directory - this is the default location for the build and
REM install directory

set CURRENT_DIR=%cd%

REM the build directory can be passed as first parameter

set BUILD_DIR=%CURRENT_DIR%\build\win7.x86_64.msvc15.release

REM the install directory can be passed as second parameter

set INSTALL_DIR=%CURRENT_DIR%\install\win7.x86_64.msvc15.release

REM create build directory if neccessary -----------------------------------------

IF EXIST %BUILD_DIR% GOTO BUILD_DIR_CREATED
    MKDIR %BUILD_DIR%
:BUILD_DIR_CREATED

REM configure cmake --------------------------------------------------------------

cd %BUILD_DIR%

cmake -G "Visual Studio 15 Win64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR% %CMAKE_DIR%^
      -DParaView_DIR=%ParaView_DIR%^
	  -DUSE_CATALYST=On
	  
msbuild INSTALL.vcxproj /p:Configuration=Release /m /v:m

cd %CURRENT_DIR%

@echo on
