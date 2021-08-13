@echo off

rem ---------------------------------------------------------------------------------------------- #
rem                              This file is part of CosmoScout VR                                #
rem     and may be used under the terms of the MIT license. See the LICENSE file for details.      #
rem                       Copyright: (c) 2019 German Aerospace Center (DLR)                        #
rem ---------------------------------------------------------------------------------------------- #

rem ---------------------------------------------------------------------------------------------- #
rem Make sure to run "git submodule update --init" before executing this script!                   #
rem Default build mode is release, if "set COSMOSCOUT_DEBUG_BUILD=true" is executed before, all    #
rem dependecies will be built in debug mode.                                                       #
rem Usage:                                                                                         #
rem    make_externals.bat [additional CMake flags, defaults to -G "Visual Studio 15 Win64"]        #
rem Examples:                                                                                      #
rem    make_externals.bat                                                                          #
rem    make_externals.bat -G "Visual Studio 15 Win64"                                              #
rem    make_externals.bat -G "Visual Studio 16 2019" -A x64                                        #
rem ---------------------------------------------------------------------------------------------- #

rem The CMake generator and other flags can be passed as parameters.
set CMAKE_FLAGS=-G "Visual Studio 15 Win64"
IF NOT "%~1"=="" (
  SET CMAKE_FLAGS=%*
)

rem Check if ComoScout VR debug build is enabled with "set COSMOSCOUT_DEBUG_BUILD=true".
IF "%COSMOSCOUT_DEBUG_BUILD%"=="true" (
  echo CosmoScout VR debug build is enabled!
  set BUILD_TYPE=debug
) else (
  set BUILD_TYPE=release
)

rem Create some required variables. ----------------------------------------------------------------

rem This directory should contain all submodules - they are assumed to reside in the subdirectory 
rem "externals" next to this script.
set EXTERNALS_DIR=%~dp0\externals

rem Get the current directory - this is the default location for the build and install directory.
set CURRENT_DIR=%cd%

rem The build directory.
set BUILD_DIR=%CURRENT_DIR%\build\windows-externals-%BUILD_TYPE%

rem The install directory.
set INSTALL_DIR=%CURRENT_DIR%\install\windows-externals-%BUILD_TYPE%

rem Create some default installation directories.
cmake -E make_directory "%INSTALL_DIR%/lib"
cmake -E make_directory "%INSTALL_DIR%/share"
cmake -E make_directory "%INSTALL_DIR%/bin"
cmake -E make_directory "%INSTALL_DIR%/include"


rem Prepare windows build
set Qt5_DIR="C:/Qt/Qt5.14.2/5.14.2/msvc2017_64/lib/cmake/Qt5"

rem MPIR --------------------------------------------------------------------------------------------
:mpir

echo .
echo Building and installing MPIR Library ...
echo .

rem mkdir "%BUILD_DIR%/mpir" && cd "%BUILD_DIR%/mpir"
msbuild %EXTERNALS_DIR%\mpir\msvc\vs17\dll_mpir_gc\dll_mpir_gc.vcxproj /p:Configuration=Release /p:Platform=x64 /p:WindowsTargetPlatformVersion=10.0.17763.0
robocopy %EXTERNALS_DIR%/mpir/dll/x64/Release %INSTALL_DIR%/mpir /s /e 

pause

cd "%CURRENT_DIR%"
echo Finished successfully.

@echo on
