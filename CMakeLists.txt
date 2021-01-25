cmake_minimum_required(VERSION 2.8)

project(VestecSampling)
find_package (Cnl REQUIRED)

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

# subdirectories
# test support (requires that GTest is available)
# ==============================================================================
 option(ENABLE_TESTING "Perform tests" OFF)

 if(ENABLE_TESTING)
    enable_testing()
    find_package(GTest REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
    add_subdirectory(test)
 endif()


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
