cmake_minimum_required(VERSION 3.9)
project(Kaiju Fortran)

# use the new version of CMP0074, this tells cmake to use <PACKAGE>_ROOT environment
#   variables when looking for packages
cmake_policy(SET CMP0074 NEW)

# add and search for pfunit
list(APPEND CMAKE_PREFIX_PATH "./external")
find_package(PFUNIT QUIET)
if(PFUNIT_FOUND)
	enable_testing()
endif()

#Set hostname & configuration defaults before going to compilers.cmake
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)

#Get the current working branch
execute_process(COMMAND git rev-parse --abbrev-ref HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_BRANCH
	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND git log -1 --format=%h WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_COMMIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE)

#Add defines for commit hash/branch
#add_compile_definitions(GITCOMMITHASH=${GIT_COMMIT_HASH})
#add_compile_definitions(GITBRANCH=${GIT_BRANCH})

#Set some configuration defaults
option(ENABLE_OMP "Enable OMP parallelization" ON )
option(ENABLE_MPI "Enable MPI parallelization" OFF)

#Vector blocking
set(bricksize 16 CACHE STRING "Vector brick size")
set_property(CACHE bricksize PROPERTY STRINGS 8 16 32)
add_compile_definitions(BRICKSIZE=${bricksize})
mark_as_advanced(bricksize)

#-------------
#Set default IC files
file(GLOB GAMICDEF src/gamera/ICs/null.F90)
set(GAMIC ${GAMICDEF} CACHE FILEPATH "Gamera IC file")

file(GLOB VOLTICDEF src/voltron/ICs/earthcmi.F90)
set(VOLTIC ${VOLTICDEF} CACHE FILEPATH "Voltron IC file")

file(GLOB TPICDEF src/chimp/tpICs/tpICstd.F90)
set(TPIC ${TPICDEF} CACHE FILEPATH "CHIMP TP IC file")

file(GLOB EBICDEF src/chimp/ebICs/ebICstd.F90)
set(EBIC ${EBICDEF} CACHE FILEPATH "CHIMP EB IC file")


#-------------
#Include user cmake file if present
include(cmake/user.cmake OPTIONAL)


#Include compile flag generator
include(cmake/compilers.cmake)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#-------------
#Print info
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)

message("-------------------------")
message("Configuration summary ...")
message("System: ${HOST}")
message("OS: ${CMAKE_SYSTEM_NAME}")
message("Processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message("Compiler: ${CMAKE_Fortran_COMPILER_ID}")
message("HDF5 Wrapper: ${HDF5_Fortran_COMPILER_EXECUTABLE}")
message("Version: ${GIT_COMMIT_HASH} / ${GIT_BRANCH}")
message("Build Type: ${CMAKE_BUILD_TYPE}")
message("Base Flags: ${CMAKE_Fortran_FLAGS}")
message("Build Flags: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message("-------------------------\n")
#message("\tBinary directory: ${PROJECT_BINARY_DIR}")

#-------------
#Start doing real work

#Start w/ base
add_subdirectory(src/base)

#-------------
#Kaiju: CHIMP
message("Adding CHIMP module ...")
add_subdirectory(src/chimp)
#Chimp executables, target "chimp"
add_custom_target(chimp ALL)
message("\tEB IC file is ${EBIC}")
message("\tTP IC file is ${TPIC}")
set(chimpX "project;psd;push;slice;trace")
foreach (EXE ${chimpX})
	message("\tAdding executable ${EXE}.x")
	add_executable(${EXE}.x src/drivers/${EXE}x.F90)
	target_link_libraries(${EXE}.x baselib chimplib)
	add_dependencies(chimp ${EXE}.x)
endforeach()

#-------------
#Kaiju: Gamera
message("Adding Gamera module ...")

#Add source
add_subdirectory(src/gamera)

#Print gamera info
message("\tBricksize is ${bricksize}")
message("\tIC file is ${GAMIC}")
add_custom_target(gamera ALL)
message("\tAdding executable gamera.x")
add_executable(gamera.x src/drivers/gamerax.F90)
target_link_libraries(gamera.x baselib gamlib)
add_dependencies(gamera gamera.x baselib)

#-------------
#Kaiju: ReMIX
message("Adding ReMIX module ...")
#Add source
add_subdirectory(src/remix)
add_custom_target(remix ALL)
#message("\tAdding executable gamera2remix.x")
#add_executable(gamera2remix.x src/drivers/gamera2remix.F90)
#target_link_libraries(gamera2remix.x baselib gamlib remixlib)
#add_dependencies(remix gamera2remix.x gamlib baselib)
message("\tAdding executable remix.x")
add_executable(remix.x src/remix/mix.F90)
target_link_libraries(remix.x baselib remixlib)
add_dependencies(remix remix.x baselib)

#-------------
#Kaiju: Voltron
message("Adding Voltron module ...")
message("\tIC file is ${VOLTIC}")
add_subdirectory(src/voltron)
add_custom_target(voltron ALL)
message("\tAdding executable voltron.x")
add_executable(voltron.x src/drivers/voltronx.F90)
target_link_libraries(voltron.x baselib gamlib remixlib voltlib)
add_dependencies(voltron voltron.x remixlib gamlib baselib)

if(PFUNIT_FOUND)
	#-------------
	#Kaiju: tests
	message("Adding tests module ...")
	#Add source
	add_subdirectory(tests)
	add_dependencies(allTests baselib gamlib remixlib voltlib)
endif()

