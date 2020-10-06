set(VESTEC_XMLS             "" ) #CACHE INTERNAL "")
set(VESTEC_MODULES          "" ) #CACHE INTERNAL "")
set(VESTEC_VTK_MODULE_FILES "" ) #CACHE INTERNAL "")
# set(VESTEC_SOURCES             "" CACHE INTERNAL "")
# set(VESTEC_HEADERS             "" CACHE INTERNAL "")
# set(VESTEC_LINKS             "" CACHE INTERNAL "")

function(vestec_register_paraview_plugin moduleName moduleXML)
	message("[VESTEC] Register module: " ${moduleName} " " ${moduleXML})
	list(APPEND VESTEC_XMLS ${moduleXML})
	set(VESTEC_XMLS ${VESTEC_XMLS} PARENT_SCOPE)
	list(APPEND VESTEC_MODULES ${moduleName})
	set(VESTEC_MODULES ${VESTEC_MODULES} PARENT_SCOPE)
endfunction()

function(vestec_add_paraview_plugin library)
	message("[VESTEC] Add module: " ${library})
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION;"
		"SOURCES;HEADERS;LINKS" ${ARGN})
	# message("sources: " ${ARG_SOURCES})
	# message("link: " ${ARG_LINKS})
	# message("xml: " ${ARG_PLUGIN_XML})
	# message("Trying to add module " ${library} " " ${ARG_PLUGIN_XML})

	vtk_module_add_module(${library}
	  SOURCES
		${ARG_SOURCES}
	  HEADERS
		${ARG_HEADERS}
	 )
	
	vtk_module_compile_options(${library} PRIVATE "/W4")
	
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

# when build against ParaView,
# vtk_module_build is already called
# by the paraview_add_plugin
# vtk_module_build(
#     MODULES
#       ${VESTEC_PROVIDED_MODULES}
#     )

file(GLOB VESTEC_PLUGIN_DIRS *)
foreach(VESTEC_PLUGIN ${VESTEC_PLUGIN_DIRS})
	if (IS_DIRECTORY ${VESTEC_PLUGIN} AND NOT ${VESTEC_PLUGIN} MATCHES "singlePlugin")
		# add_subdirectory(${VESTEC_PLUGIN})
		include(${VESTEC_PLUGIN}/VestecFilterRegister.cmake)
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

message("SCAN")

paraview_plugin_scan(
  PLUGIN_FILES
    "singlePlugin/paraview.plugin"
  PROVIDES_PLUGINS
	VESTEC_PLUGINS
  ENABLE_BY_DEFAULT ON
  )

message("BUILD")

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

#  message("INIT MODULES")
#  set(VESTEC_INIT "1")
 
#  foreach(VESTEC_PLUGIN ${VESTEC_PLUGIN_DIRS})
# 	 if (IS_DIRECTORY ${VESTEC_PLUGIN} AND NOT ${VESTEC_PLUGIN} MATCHES "singlePlugin")
# 		 # add_subdirectory(${VESTEC_PLUGIN})
# 		 include(${VESTEC_PLUGIN}/VestecFilterAdd.cmake)
# 	 endif()
#  endforeach()

install(
  TARGETS
    VESTECPVTargets
  DESTINATION
    "${VESTEC_INSTALL_PLUGIN_DIR}"
  )

