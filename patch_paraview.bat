@echo off

rem ---------------------------------------------------------------------------------------------- #
rem                              This file is part of CosmoScout VR                                #
rem     and may be used under the terms of the MIT license. See the LICENSE file for details.      #
rem                       Copyright: (c) 2019 German Aerospace Center (DLR)                        #
rem ---------------------------------------------------------------------------------------------- #

rem ---------------------------------------------------------------------------------------------- #
rem Make sure to run "git submodule update --init" before executing this script!                   #
rem ---------------------------------------------------------------------------------------------- #

rem Create some required variables. ----------------------------------------------------------------

rem This directory should contain all submodules - they are assumed to reside in the subdirectory 
rem "externals" next to this script.
set EXTERNALS_DIR=%~dp0\externals

rem Paraview --------------------------------------------------------------------------------------------

cd "%EXTERNALS_DIR%/ttk/paraview/patch/"
patch-paraview-msvc.cmd "%EXTERNALS_DIR%/paraview-5.6"

@echo on
