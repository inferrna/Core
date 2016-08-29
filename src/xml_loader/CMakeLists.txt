if(APPLE)
        SET(CMAKE_INSTALL_RPATH @loader_path)
else(APPLE)
        SET(CMAKE_INSTALL_RPATH "$ORIGIN/") #If multiple directories are needed this can be used for example: “$ORIGIN/:$ORIGIN/$
endif(APPLE)

include_directories(${YAF_INCLUDE_DIRS})

add_executable(thebounty-xml xml-loader.cc)
target_link_libraries(thebounty-xml yafaraycore)

install (TARGETS thebounty-xml RUNTIME DESTINATION ${YAF_BIN_DIR})
