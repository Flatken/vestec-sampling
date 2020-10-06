# fix for typo in ParaView 5.7
# https://discourse.paraview.org/t/install-headers-for-vtk-modules-build-by-paraview-5-7-0/2721
set(_paraview_add_plugin_MODULES_FILES TRUE)

message("[VESTEC] [DEBUG] into plugins/singlePlugin/CMakeLists.txt")
message("Modules_NAMES "  ${VESTEC_MODULES})
message("MODULE_FILES " ${VESTEC_VTK_MODULE_FILES})
# message("Modules_SOURCES "  ${VESTEC_SOURCES})
# message("Modules_LINKS "  ${VESTEC_LINKS})
message("Modules_XMLS "  ${VESTEC_XMLS})
message("===")

paraview_add_plugin(VestecPlugins
  VERSION "1.0"
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  MODULES
    ${VESTEC_MODULES}
  MODULE_FILES
    ${VESTEC_VTK_MODULE_FILES}
  MODULE_ARGS # BUG: see above
    INSTALL_HEADERS ON
  SERVER_MANAGER_XML
    ${VESTEC_XMLS}
  )

if(NOT "${VESTEC_INSTALL_PLUGIN_DIR}" STREQUAL "")
  install(
    TARGETS
      VestecPlugins
    DESTINATION
      "${VESTEC_INSTALL_PLUGIN_DIR}"
    )
endif()

