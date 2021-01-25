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

rem Paraview --------------------------------------------------------------------------------------------
:paraview

rem cd "%EXTERNALS_DIR%/ttk/paraview/patch/"
rem patch-paraview-msvc.cmd "%EXTERNALS_DIR%/paraview-5.6"

echo .
echo Building and installing Paraview 5.8 ...
echo .

cmake -E make_directory "%BUILD_DIR%/paraview" && cd "%BUILD_DIR%/paraview"
cmake %CMAKE_FLAGS% -DCMAKE_INSTALL_PREFIX="%INSTALL_DIR%"^
      -DCMAKE_INSTALL_LIBDIR=lib^
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON^
      -DPARAVIEW_USE_PYTHON=ON^
      -DVTK_PYTHON_VERSION=3^
      -DPARAVIEW_USE_QT=ON^
      -DPARAVIEW_USE_MPI=ON^
      -DPARAVIEW_USE_VTKM=OFF^
      -DQt5_DIR="C:/Qt/Qt5.14.2/5.14.2/msvc2017_64/lib/cmake/Qt5"^
      -DCMAKE_BUILD_TYPE="%BUILD_TYPE%" "%EXTERNALS_DIR%/paraview-5.6"
cmake --build . --config "%BUILD_TYPE%" --target install --parallel 12

rem # EIGEN -----------------------------------------------------------------------------------------
:eigen

echo .
echo Building and installing Eigen

cmake -E make_directory "%BUILD_DIR%/eigen" && cd "%BUILD_DIR%/eigen"
cmake %CMAKE_FLAGS% -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR%^
       -DQt5_DIR=C:/Qt/Qt5.14.2/5.14.2/msvc2017_64/lib/cmake/Qt5^
       -DCMAKE_BUILD_TYPE=%BUILD_TYPE% "%EXTERNALS_DIR%/eigen"
cmake --build . --config %BUILD_TYPE% --target install --parallel 12

echo .

rem # TTK -----------------------------------------------------------------------------------------
:ttk

echo .
echo Building and installing TTK
echo .

cmake -E remove_directory "%EXTERNALS_DIR%/ttk/paraview/WRLExporter"
cmake -E remove_directory "%EXTERNALS_DIR%/ttk/core/vtk/ttkWRLExporter"

cmake -E make_directory "%BUILD_DIR%/ttk" && cd "%BUILD_DIR%/ttk"
cmake %CMAKE_FLAGS% -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR%^
      -DParaView_DIR="%BUILD_DIR%/paraview"^
      -DEigen3_DIR="%INSTALL_DIR%/share/eigen3/cmake"^
      -DTTK_INSTALL_PLUGIN_DIR="%INSTALL_DIR%/lib/plugins"^
      -DTTK_ENABLE_ZLIB=OFF^
      -DTTK_ENABLE_KAMIKAZE=On^
      -DTTK_ENABLE_MPI=ON^
      -DTTK_BUILD_STANDALONE_APPS=OFF^
      -DVTK_MODULE_ENABLE_ttkCinemaImaging=DONT_WANT^
      -DVTK_MODULE_ENABLE_ttkUserInterfaceBase=DONT_WANT^
	    -DQt5_DIR=C:/Qt/Qt5.14.2/5.14.2/msvc2017_64/lib/cmake/Qt5^
      -DCMAKE_CXX_FLAGS="/bigobj /EHsc /UBOOST_NO_EXCEPTIONS"^
      -DCMAKE_BUILD_TYPE=%BUILD_TYPE% "%EXTERNALS_DIR%/ttk"
cmake --build . --config %BUILD_TYPE% --target install --parallel 12

pause

cd "%CURRENT_DIR%"
echo Finished successfully.

@echo on
