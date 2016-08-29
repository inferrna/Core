include_directories(${YAF_INCLUDE_DIRS})

add_executable(thebounty-xml xml-loader.cc)
target_link_libraries(thebounty-xml yafaraycore)

install (TARGETS thebounty-xml RUNTIME DESTINATION ${YAF_BIN_DIR})
