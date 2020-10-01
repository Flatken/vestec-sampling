# fix for typo in ParaView 5.7
# https://discourse.paraview.org/t/install-headers-for-vtk-modules-build-by-paraview-5-7-0/2721
set(_paraview_add_plugin_MODULES_FILES TRUE)

message("Modules "  ${VESTEC_MODULES})
message("MODULE_FILES " ${VESTEC_VTK_MODULE_FILES})
paraview_add_plugin(VestecSampling
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

install(
  TARGETS
    VestecSampling
  DESTINATION
    "lib/plugins"
)

