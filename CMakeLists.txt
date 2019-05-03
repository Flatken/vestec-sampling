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
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    # compiler settings for Visual Studio (only tested 10, aka MSVC 1600)
    if(${MSVC_VERSION} EQUAL   1600 OR
       ${MSVC_VERSION} GREATER 1600)
        # make use of SSE2 instructions
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
        # disable warning C4251: "'T' needs to have dll-interface to be used ..."
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4251")
        # disable warning C4275: "non dll-interface class 'T' usedas base ..."
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4275")
		# Enable Multicore Compiling
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    else()
        message(WARNING "MSVC_VERSION '${MSVC_VERSION}' not recognized - "
                "can not configure compiler flags")
    endif()
else()
    add_definitions("-std=c++11")
    set(BUILD_SHARED_LIBS "Off")
endif()

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
