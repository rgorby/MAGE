file(GLOB caseTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/cases/*.pf")
set(caseTestLibs baselib gamlib remixlib voltlib)
add_pfunit_ctest (caseTests TEST_SOURCES ${caseTestFiles} LINK_LIBRARIES ${caseTestLibs} MAX_PES 1)

file(GLOB caseMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/cases_mpi/*.pf")
set(caseMpiTestLibs baselib gamlib remixlib voltlib basempilib gammpilib)
add_pfunit_ctest (caseMpiTests TEST_SOURCES ${caseMpiTestFiles} LINK_LIBRARIES ${caseMpiTestLibs} MAX_PES 64)

file(GLOB gamTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/gamera/*.pf")
set(gamTestLibs baselib gamlib)
add_pfunit_ctest (gamTests TEST_SOURCES ${gamTestFiles} LINK_LIBRARIES ${gamTestLibs} MAX_PES 1)

file(GLOB mixTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/remix/*.pf")
set(mixTestLibs baselib gamlib remixlib voltlib)
add_pfunit_ctest (mixTests TEST_SOURCES ${mixTestFiles} LINK_LIBRARIES ${mixTestLibs} MAX_PES 1)

file(GLOB voltTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/voltron/*.pf")
set(voltTestLibs baselib gamlib remixlib voltlib)
add_pfunit_ctest (voltTests TEST_SOURCES ${voltTestFiles} LINK_LIBRARIES ${voltTestLibs} MAX_PES 1)

file(GLOB baseMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/base_mpi/*.pf")
set(baseMpiTestLibs baselib basempilib)
add_pfunit_ctest (baseMpiTests TEST_SOURCES ${baseMpiTestFiles} LINK_LIBRARIES ${baseMpiTestLibs} MAX_PES 64)

file(GLOB gamMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/gamera_mpi/*.pf")
set(gamMpiTestLibs baselib gamlib basempilib gammpilib)
add_pfunit_ctest (gamMpiTests TEST_SOURCES ${gamMpiTestFiles} LINK_LIBRARIES ${gamMpiTestLibs} MAX_PES 64)

file(GLOB voltMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/voltron_mpi/*.pf")
set(voltMpiTestLibs baselib gamlib basempilib gammpilib voltlib voltmpilib remixlib)
add_pfunit_ctest (voltMpiTests TEST_SOURCES ${voltMpiTestFiles} LINK_LIBRARIES ${voltMpiTestLibs} MAX_PES 17)

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

# all tests in one
file(GLOB_RECURSE allTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.pf")
set(allTestLibs baselib gamlib remixlib voltlib basempilib gammpilib voltmpilib)
add_pfunit_ctest (allTests TEST_SOURCES ${allTestFiles} LINK_LIBRARIES ${allTestLibs} MAX_PES 64)

add_dependencies(allTests CopyInputs CopyGrids caseTests gamTests mixTests voltTests gamMpiTests baseMpiTests caseMpiTests voltMpiTests)

