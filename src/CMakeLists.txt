
# Dependencies ==============================================================
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(BUILD_SHARED_LIBS ON)

find_package(ParaView REQUIRED)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP)

# Target for the emulator====================================================
set(TARGET_NAME VestecCatalystEmulator)

# ============== ParaView Catalyst dependecy =================================
set(USE_CATALYST ON CACHE BOOL "Link the emulator against paraview catalyst")

# Headers / Sources for the emulator
# ===========================================================================
file(GLOB CPP_SOURCE_FILES *.cpp)
file(GLOB CPP_HEADER_FILES *.hpp)

# Headers / Sources for the catalyst adaptors
# ===========================================================================
file(GLOB CPP_SOURCE_FILES_ADAPTORS Adaptors/*.cpp)
file(GLOB CPP_HEADER_FILES_ADAPTORS Adaptors/*.hpp)

# Build VESTEC ParaView plugins
# ==============================================================================	
set(BUILD_VESTEC_PLUGINS ON CACHE BOOL "Build VESTEC plugins for ParaView and Catalyst")
if(BUILD_VESTEC_PLUGINS)
	add_subdirectory(plugins)
endif()	

# Build Adaptor library
# ==============================================================================
if(USE_CATALYST)
	add_library(Adaptor ${CPP_SOURCE_FILES_ADAPTORS} ${CPP_HEADER_FILES_ADAPTORS})
	target_link_libraries(Adaptor  
		ParaView::Catalyst 
		ParaView::PythonCatalyst 
		${VTK_LIBRARIES}
		)
	add_definitions("-DUSE_CATALYST")
	if(NOT PARAVIEW_USE_MPI)
		message(SEND_ERROR "ParaView must be built with MPI enabled")
	endif()
	
	install(TARGETS Adaptor
		RUNTIME DESTINATION "bin"
		LIBRARY DESTINATION "lib"
		ARCHIVE DESTINATION "lib"
	)

else()
	find_package(MPI REQUIRED)
	include_directories(${MPI_CXX_INCLUDE_PATH})
endif()



# Build VESTEC evaluation tool
# ==============================================================================	
add_executable(${TARGET_NAME} ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES})
if(USE_CATALYST)
	IF (WIN32)
		target_link_libraries(${TARGET_NAME} LINK_PRIVATE Adaptor OpenMP::OpenMP_CXX ParaView::Catalyst ParaView::PythonCatalyst ${VTK_LIBRARIES})
	ELSE()
		target_link_libraries(${TARGET_NAME} LINK_PRIVATE Adaptor stdc++fs OpenMP::OpenMP_CXX ParaView::Catalyst ParaView::PythonCatalyst ${VTK_LIBRARIES})
	ENDIF()
else()
	IF (WIN32)
		target_link_libraries(${TARGET_NAME} LINK_PRIVATE ${MPI_LIBRARIES} OpenMP::OpenMP_CXX)
	ELSE()
		target_link_libraries(${TARGET_NAME} LINK_PRIVATE ${MPI_LIBRARIES} stdc++fs OpenMP::OpenMP_CXX)
    ENDIF()
endif()

target_include_directories(${TARGET_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/..
)

# Install
# ==============================================================================
install(TARGETS ${TARGET_NAME}
		RUNTIME DESTINATION "bin"
		LIBRARY DESTINATION "lib"
		ARCHIVE DESTINATION "lib"
)

install(
  DIRECTORY ${EXTERNALS_DIR}/lib/
  DESTINATION "lib"
  USE_SOURCE_PERMISSIONS
)

install(
  DIRECTORY ${EXTERNALS_DIR}/bin/
  DESTINATION "bin"
  USE_SOURCE_PERMISSIONS
)

install(
    FILES 
        ${CPP_HEADER_FILES}
    DESTINATION 
        "include/VestecSampling"
)

