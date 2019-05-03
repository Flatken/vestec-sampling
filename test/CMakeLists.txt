# ==============================================================================
# Perform unit tests
# ==============================================================================

set(LIBRARIES
    VistaAtmosphere
    pthread
    ${GTEST_BOTH_LIBRARIES}
	${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${GLEW_LIBRARIES}
    ${OPENSG_LIBRARIES}
    ${VISTACORELIBS_LIBRARIES}
)


set(ENABLED_TESTS
    test_Dummy
)

foreach(TEST ${ENABLED_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ${LIBRARIES})
    add_test(${TEST} ${TEST})
endforeach()
