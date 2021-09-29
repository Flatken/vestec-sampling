cmake_minimum_required(VERSION 2.8)

project(VestecSampling)

# ensure local modules (for dependencies etc.) are found
list(
    APPEND CMAKE_MODULE_PATH
    "$ENV{CMAKE_FIND_SCRIPTS}"
)

# CMAKE_BUILD_TYPE must be set (except for Visual Studio)
# ==============================================================================
if(NOT MSVC)
    if(NOT CMAKE_BUILD_TYPE OR
       (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release" AND
        NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug"      ))

        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release or Debug" FORCE)
    endif()
endif()

# Compiler settings
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/lib"
)

# Subdirectories
# ==============================================================================
add_subdirectory(src)
add_subdirectory(scripts)
 
# build doxygen documentation
# ==============================================================================
option(ENABLE_DOCBUILD "Create doxygen documentation" OFF)

if(ENABLE_DOCBUILD)
    find_package(Doxygen REQUIRED)
    add_subdirectory(doc)
endif()
