
# Catalyst pipelines
file(GLOB CATALYST_SCRIPTS "catalyst/*.py") 
install(FILES ${CATALYST_SCRIPTS} DESTINATION "bin")

# start script files -----------------------------------------------------------------
if(WIN32)
    file(GLOB BAT_FILES "windows/*.bat") 
    file(GLOB CMD_FILES "windows/*.cmd") 
    install(FILES ${BAT_FILES} DESTINATION "bin")
    install(FILES ${CMD_FILES} DESTINATION "bin")
else()
    file(GLOB SH_FILES "linux/*.sh") 
    file(GLOB PBS_FILES "linux/*.pbs")
    install(FILES ${SH_FILES} DESTINATION "bin"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)
    install(FILES ${PBS_FILES} DESTINATION "bin"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)
endif()

