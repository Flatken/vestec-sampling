# Function to create the ParaView plugin for a VESTEC ParaView plugin

# Options:
# SOURCES: Specify the list of source files for the library
# LINK: Specify the link dependencies of the library
# PLUGIN_VERSION: Specify ParaView plugin version
# PLUGIN_XML: Specify ParaView plugin XML file
#
# ttk_add_paraview_plugin(<library_name>
#     SOURCES <source list>
#     PLUGIN_VERSION <version>
#     PLUGIN_XML <XML file for the plugin>
#     LINK <libraries to link>)
#
function(vestec_add_paraview_plugin library)
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION"
		"SOURCES;LINK" ${ARGN})

	# Build the ParaView plugin, if we're building with ParaView
	# and this wrapper library has a plugin
	set(plugin_name "${library}Plugin")
	add_paraview_plugin(${plugin_name} ${ARG_PLUGIN_VERSION}
		SERVER_MANAGER_XML ${ARG_PLUGIN_XML}
		SERVER_MANAGER_SOURCES ${ARG_SOURCES})

	target_link_libraries(${plugin_name} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})
	
	#if (MSVC)
	#	target_compile_definitions(${plugin_name} PUBLIC VESTEC_PLUGIN)				
	#endif()	

	target_include_directories(${plugin_name} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

	install(TARGETS ${plugin_name} DESTINATION "lib/plugins")
endfunction()

file(GLOB VESTEC_PLUGIN_DIRS *)
foreach(VESTEC_PLUGIN ${VESTEC_PLUGIN_DIRS})
	if (IS_DIRECTORY ${VESTEC_PLUGIN})
		add_subdirectory(${VESTEC_PLUGIN})
	endif()
endforeach()