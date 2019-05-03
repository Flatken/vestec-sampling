# ==============================================================================
# install doxygen documentation
# ==============================================================================

# configure
# ==============================================================================
configure_file("doxyfile.in" "${CMAKE_BINARY_DIR}/doc/doxyfile")

# build
# ==============================================================================
add_custom_target(
    doc
    ${DOXYGEN_EXECUTABLE} "${CMAKE_BINARY_DIR}/doc/doxyfile"
    WORKING_DIRECTORY     "${CMAKE_BINARY_DIR}/doc"
)

# install
# ==============================================================================
install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html DESTINATION share/doc OPTIONAL)
