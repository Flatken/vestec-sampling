set(VESTEC_XMLS             "" CACHE INTERNAL "")
set(VESTEC_MODULES          "" CACHE INTERNAL "")
set(VESTEC_VTK_MODULE_FILES "" CACHE INTERNAL "")


function(vestec_add_paraview_plugin library)
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION"
		"SOURCES;LINK" ${ARGN})

	# Build the ParaView plugin, if we're building with ParaView
	# and this wrapper library has a plugin
	# set(plugin_name "${library}Plugin")
	#add_paraview_plugin(${plugin_name} ${ARG_PLUGIN_VERSION}
	#	SERVER_MANAGER_XML ${ARG_PLUGIN_XML}
	#	SERVER_MANAGER_SOURCES ${ARG_SOURCES})
	
	message("Trying to add module " ${library} " " ${ARG_PLUGIN_XML})
	list(APPEND VESTEC_XMLS ${ARG_PLUGIN_XML})
	list(APPEND VESTEC_MODULES ${library})

	#target_link_libraries(${plugin_name} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})
	
	#if (MSVC)
	#	target_compile_definitions(${plugin_name} PUBLIC VESTEC_PLUGIN)				
	#endif()	

	#target_include_directories(${plugin_name} PUBLIC
	#	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

	#install(TARGETS ${plugin_name} DESTINATION "lib/plugins")
endfunction()

file(GLOB VESTEC_PLUGIN_DIRS *)
foreach(VESTEC_PLUGIN ${VESTEC_PLUGIN_DIRS})
	if (IS_DIRECTORY ${VESTEC_PLUGIN} AND NOT ${VESTEC_PLUGIN} MATCHES "singlePlugin")
		add_subdirectory(${VESTEC_PLUGIN})
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
  PLUGINS
    ${VESTEC_PLUGINS}
)

install(
  TARGETS
    TTKPVTargets
  DESTINATION
    "${TTK_INSTALL_PLUGIN_DIR}"
  )

