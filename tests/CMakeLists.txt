file(GLOB_RECURSE testFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.pf")
set(testLibs baselib gamlib remixlib voltlib)
add_pfunit_ctest (allTests TEST_SOURCES ${testFiles} LINK_LIBRARIES ${testLibs} MAX_PES 1)

# copy input xml files and grid h5 files to binary folder for running
file(GLOB_RECURSE inputFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.xml")
foreach(inputFile ${inputFiles})
  get_filename_component(inputFileJustName ${inputFile} NAME)
  add_custom_command(
    OUTPUT "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${inputFileJustName}"
    COMMAND cmake -E copy
      "${inputFile}"
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${inputFileJustName}"
    DEPENDS "${inputFile}"
  )
  list(APPEND inputFilesDest "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${inputFileJustName}")
endforeach(inputFile)
add_custom_target(CopyInputs DEPENDS ${inputFilesDest})

file(GLOB_RECURSE gridFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h5")
foreach(gridFile ${gridFiles})
  get_filename_component(gridFileJustName ${gridFile} NAME)
  add_custom_command(
    OUTPUT "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${gridFileJustName}"
    COMMAND cmake -E copy
      "${gridFile}"
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${gridFileJustName}"
    DEPENDS "${gridFile}"
  )
  list(APPEND gridFilesDest "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${gridFileJustName}")
endforeach(gridFile)
add_custom_target(CopyGrids DEPENDS ${gridFilesDest})

add_dependencies(allTests CopyInputs CopyGrids)
