set(VESTEC_XMLS             "" CACHE INTERNAL "")
set(VESTEC_MODULES          "" CACHE INTERNAL "")
set(VESTEC_VTK_MODULE_FILES "" CACHE INTERNAL "")
# set(VESTEC_SOURCES             "" CACHE INTERNAL "")
# set(VESTEC_HEADERS             "" CACHE INTERNAL "")
# set(VESTEC_LINKS             "" CACHE INTERNAL "")

function(vestec_add_paraview_plugin library)
	message("argn: " ${ARGN})
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION;"
		"SOURCES;HEADERS;LINKS" ${ARGN})
	message("sources: " ${ARG_SOURCES})
	message("link: " ${ARG_LINKS})
	message("xml: " ${ARG_PLUGIN_XML})
	message("Trying to add module " ${library} " " ${ARG_PLUGIN_XML})

	if(NOT TARGET ${library})
		vtk_module_add_module(${library}
		  SOURCES
			${ARG_SOURCES}
		  HEADERS
			${ARG_HEADERS}
		  )
	
		# vtk_module_compile_options(${library} PRIVATE "")
	
		vtk_module_link(${library}
		  PUBLIC
			${VTK_LIBRARIES}
			${ARG_LINKS}
		  )
	
		install(
		  TARGETS
			${library}
		  EXPORT
			VESTECVTKTargets
		  RUNTIME DESTINATION
			bin
		  ARCHIVE DESTINATION
			lib
		  LIBRARY DESTINATION
			lib
		  )
	endif()
	

	# Build the ParaView plugin, if we're building with ParaView
	# and this wrapper library has a plugin
	
	#add_paraview_plugin(${plugin_name} ${ARG_PLUGIN_VERSION}
	#	SERVER_MANAGER_XML ${ARG_PLUGIN_XML}
	#	SERVER_MANAGER_SOURCES ${ARG_SOURCES})
	
	# message("Trying to add module " ${library} " " ${ARG_PLUGIN_XML})
	list(APPEND VESTEC_XMLS ${ARG_PLUGIN_XML})
	set(VESTEC_XMLS ${VESTEC_XMLS} PARENT_SCOPE)
	list(APPEND VESTEC_MODULES ${library})
	set(VESTEC_MODULES ${VESTEC_MODULES} PARENT_SCOPE)
	# list(APPEND VESTEC_SOURCES ${ARG_SOURCES})
	# set(VESTEC_SOURCES ${VESTEC_SOURCES} PARENT_SCOPE)
	# list(APPEND VESTEC_SOURCES ${ARG_SOURCES})
	# set(VESTEC_SOURCES ${VESTEC_SOURCES} PARENT_SCOPE)
	# list(APPEND VESTEC_LINKS ${ARG_LINK})
	# set(VESTEC_LINKS ${VESTEC_LINKS} PARENT_SCOPE)

	# message("INTO-funct Modules_SOURCES "  ${VESTEC_SOURCES})
	# message("INTO-funct Modules_LINKS "  ${VESTEC_LINKS})
	# message("INTO-funct Modules_XMLS "  ${VESTEC_XMLS})
	# message("")

	#target_link_libraries(${plugin_name} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})
	
	#if (MSVC)
	#	target_compile_definitions(${plugin_name} PUBLIC VESTEC_PLUGIN)				
	#endif()	

	#target_include_directories(${plugin_name} PUBLIC
	#	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

	#install(TARGETS ${plugin_name} DESTINATION "lib/plugins")
endfunction()

# message("BEFORE Modules_SOURCES "  ${VESTEC_SOURCES})
# message("BEFORE Modules_LINKS "  ${VESTEC_LINKS})
# message("BEFORE Modules_XMLS "  ${VESTEC_XMLS})
# message("")

vtk_module_find_modules(
  VESTEC_VTK_MODULE_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )

message("VESTEC_VTK_MODULE_FILE " ${VESTEC_VTK_MODULE_FILE})

vtk_module_scan(
  MODULE_FILES
    ${VESTEC_VTK_MODULE_FILE}
  PROVIDES_MODULES
    VESTEC_PROVIDED_MODULES
  WANT_BY_DEFAULT
    ON
  )

message("VESTEC_PROVIDED_MODULES " ${VESTEC_PROVIDED_MODULES})

# vtk_module_build(
#     MODULES
#       ${VESTEC_PROVIDED_MODULES}
#     )

file(GLOB VESTEC_PLUGIN_DIRS *)
foreach(VESTEC_PLUGIN ${VESTEC_PLUGIN_DIRS})
	if (IS_DIRECTORY ${VESTEC_PLUGIN} AND NOT ${VESTEC_PLUGIN} MATCHES "singlePlugin")
		# add_subdirectory(${VESTEC_PLUGIN})
		include(${VESTEC_PLUGIN}/CMakeLists.txt)
		# message("INTO-FOR Modules_SOURCES "  ${VESTEC_SOURCES})
		# message("INTO-FOR Modules_LINKS "  ${VESTEC_LINKS})
		# message("INTO-FOR Modules_XMLS "  ${VESTEC_XMLS})
		# message("")
		# get_filename_component(DNAME ${VESTEC_PLUGIN} NAME)
		# list(APPEND VESTEC_MODULES ${DNAME})
		# list(APPEND VESTEC_MODULES ${VESTEC_PLUGIN} ";")
		list(APPEND VESTEC_VTK_MODULE_FILES ${VESTEC_PLUGIN}/vtk.module)
	endif()
endforeach()


paraview_plugin_scan(
  PLUGIN_FILES
    "singlePlugin/paraview.plugin"
  PROVIDES_PLUGINS
	VESTEC_PLUGINS
  ENABLE_BY_DEFAULT ON
  )

  paraview_plugin_build(
	TARGET
	  VESTECPVTargets
	PLUGINS
	  ${VESTEC_PLUGINS}
	AUTOLOAD
	  ${VESTEC_PLUGINS}
	RUNTIME_DESTINATION
	  "${CMAKE_INSTALL_BINDIR}"
	LIBRARY_DESTINATION
	  "${CMAKE_INSTALL_LIBDIR}"
	LIBRARY_SUBDIRECTORY
	  "${PARAVIEW_PLUGIN_SUBDIR}"
	)

install(
  TARGETS
    VESTECPVTargets
  DESTINATION
    "${VESTEC_INSTALL_PLUGIN_DIR}"
  )

