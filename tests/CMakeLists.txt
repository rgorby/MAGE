file(GLOB_RECURSE testFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.pf")
set(testLibs baselib gamlib remixlib voltlib)
add_pfunit_ctest (allTests TEST_SOURCES ${testFiles} LINK_LIBRARIES ${testLibs} MAX_PES 1)

