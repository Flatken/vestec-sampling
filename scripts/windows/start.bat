@echo off
set PATH=%PATH%;D:\SoftwareEntwicklung\Catalyst\install\win7.x86_64.msvc15.release\bin

mpiexec /genv PATH D:\SoftwareEntwicklung\Catalyst\install\win7.x86_64.msvc15.release\bin -hosts 1 localhost -cores 1 VestecSampling.exe "D:\\vr_data\\VESTEC\\diseases\\Output"
@echo on
