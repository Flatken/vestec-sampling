set(VESTEC_XMLS             "" ) 
set(VESTEC_MODULES          "" )
set(VESTEC_VTK_MODULE_FILES "" )

function(vestec_register_paraview_plugin moduleName moduleXML)
	message("[VESTEC] Register module: " ${moduleName} " " ${moduleXML})
	# replace variables of original XML file and store generated file in build dir
	configure_file(${CMAKE_CURRENT_LIST_DIR}/${moduleXML} ${CMAKE_BINARY_DIR}/paraview/xmls/${moduleXML})

	# use generated xml file for pv plugin
	list(APPEND VESTEC_XMLS ${CMAKE_BINARY_DIR}/paraview/xmls/${moduleXML})
	set(VESTEC_XMLS ${VESTEC_XMLS} PARENT_SCOPE)
	list(APPEND VESTEC_MODULES ${moduleName})
	set(VESTEC_MODULES ${VESTEC_MODULES} PARENT_SCOPE)
endfunction()

function(vestec_add_paraview_plugin library)
	message("[VESTEC] Add module: " ${library})
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION;"
		"SOURCES;HEADERS;LINKS" ${ARGN})
	
	vtk_module_add_module(${library}
	  SOURCES
		${ARG_SOURCES}
	  HEADERS
		${ARG_HEADERS}
	 )
	
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

vtk_module_find_modules(
  VESTEC_VTK_MODULE_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )

vtk_module_scan(
  MODULE_FILES
    ${VESTEC_VTK_MODULE_FILE}
  PROVIDES_MODULES
    VESTEC_PROVIDED_MODULES
  WANT_BY_DEFAULT
    ON
  )

file(GLOB VESTEC_PLUGIN_DIRS *)
foreach(VESTEC_PLUGIN ${VESTEC_PLUGIN_DIRS})
	if (IS_DIRECTORY ${VESTEC_PLUGIN} AND NOT ${VESTEC_PLUGIN} MATCHES "singlePlugin")
		include(${VESTEC_PLUGIN}/VestecFilterRegister.cmake)
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

if(NOT DEFINED VESTEC_PYTHON_MODULE_DIR)
 if (IS_ABSOLUTE ${VTK_PYTHONPATH})
   file(RELATIVE_PATH
	 PYTHON_SITE_PACKAGES_SUFFIX
	 ${VTK_PREFIX_PATH}
	 ${VTK_PYTHONPATH}
	 )
 else()
   set(PYTHON_SITE_PACKAGES_SUFFIX
	 ${VTK_PYTHONPATH}
	 )
 endif()
 set(VESTEC_PYTHON_MODULE_DIR
   ${PYTHON_SITE_PACKAGES_SUFFIX}
   CACHE
   PATH
   "Where the TTK python module in installed"
   )
endif()

vtk_module_wrap_python(
    MODULES
      ${VESTEC_MODULES}
    WRAPPED_MODULES
      VESTEC_WRAPPED_MODULES
    INSTALL_EXPORT
      VESTECPython
    PYTHON_PACKAGE
      "vestec"
    MODULE_DESTINATION
      "${VESTEC_PYTHON_MODULE_DIR}"
    CMAKE_DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/vestec/python"
    LIBRARY_DESTINATION
      "${CMAKE_INSTALL_LIBDIR}"
    INSTALL_HEADERS
      ON
    TARGET
      VESTEC::Python
)
# Generate __init__.py
set(InitContent "from __future__ import absolute_import\n")
foreach(MODULE ${VESTEC_MODULES})
  string(APPEND InitContent "from .${MODULE} import *\n")
endforeach()
file(GENERATE
  OUTPUT
  	"${CMAKE_BINARY_DIR}/${VESTEC_PYTHON_MODULE_DIR}/vestec/__init__.py"
  CONTENT
	"${InitContent}"
  )
install(
  FILES
	"${CMAKE_BINARY_DIR}/${VESTEC_PYTHON_MODULE_DIR}/vestec/__init__.py"
  DESTINATION
	"${VESTEC_PYTHON_MODULE_DIR}/vestec"
  )

install(
  TARGETS
    VESTECPVTargets
  DESTINATION
    "${VESTEC_INSTALL_PLUGIN_DIR}"
)

