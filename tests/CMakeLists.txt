include_directories(${PROJECT_SOURCE_DIR}/external/pfunit-mpi/mod)
set(PFUNIT_DIR ${PROJECT_SOURCE_DIR}/external/pfunit-mpi)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)
file(WRITE ${PROJECT_BINARY_DIR}/generated/testSuites.inc "")

include_directories(${PFUNIT_DIR}/mod)

set(_test_sources)
file(GLOB_RECURSE testFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.pf")
foreach(_test ${testFiles})
    message("Processing Test File: ${_test}")

    get_filename_component(_test_noExt ${_test} NAME_WE)

    add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/generated/${_test_noExt}.F90
        COMMAND python ${PFUNIT_DIR}/bin/pFUnitParser.py ${_test} ${PROJECT_BINARY_DIR}/generated/${_test_noExt}.F90
        DEPENDS ${_test}
        )

    set(_test_sources ${_test_sources} ${PROJECT_BINARY_DIR}/generated/${_test_noExt}.F90)

    file(APPEND ${PROJECT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${_test_noExt}_suite)\n")
endforeach()

set_source_files_properties(${PFUNIT_DIR}/include/driver.F90 PROPERTIES GENERATED 1)

include_directories(${PROJECT_BINARY_DIR}/generated)

add_executable(
    pftest_alltests
    ${PFUNIT_DIR}/include/driver.F90
    ${_test_sources}
    )

target_link_libraries(
    pftest_alltests
    ${PFUNIT_DIR}/lib/libpfunit.a
    baselib
    gamlib
    remixlib
    chimplib
    )

add_test(
    NAME pftest_alltests
    COMMAND mpirun -np 13 ${PROJECT_BINARY_DIR}/pftest_alltests
    )

