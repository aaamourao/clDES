set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_MAN YES)

set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

configure_file(${doxyfile_in} ${doxyfile} @ONLY)

add_custom_target(doc
                 COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                 WORKING_DIRECTORY
                 ${CMAKE_CURRENT_BINARY_DIR}
                 COMMENT "Generating documentation with Doxygen"
                 VERBATIM
)
# doxygen_add_docs(
#         doxygen
#         ${PROJECT_SOURCE_DIR}
#         COMMENT "Generate clDES man pages"
# )
