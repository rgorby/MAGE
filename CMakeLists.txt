cmake_minimum_required(VERSION 3.10)
project(Kaiju Fortran)
enable_testing()

#Set hostname & configuration defaults before going to compilers.cmake
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)

# Get the current working branch
execute_process(COMMAND git rev-parse --abbrev-ref HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_BRANCH
	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND git log -1 --format=%h WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_COMMIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE)

#Set some configuration defaults
option(ENABLE_OMP "Enable OMP parallelization" ON )
option(ENABLE_MPI "Enable MPI parallelization" OFF)

#Vector blocking
set(bricksize 16 CACHE STRING "Vector brick size")
set_property(CACHE bricksize PROPERTY STRINGS 8 16 32)
add_compile_definitions(BRICKSIZE=${bricksize})
mark_as_advanced(bricksize)

#Include compile flag generator
include(cmake/compilers.cmake)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#-------------
#Include user cmake file if present
include(cmake/user.cmake OPTIONAL)

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
set(chimpX "project;psd;push;slice")
foreach (EXE ${chimpX})
	message("\tAdding executable ${EXE}.x")
	add_executable(${EXE}.x src/drivers/${EXE}x.F90)
	target_link_libraries(${EXE}.x baselib chimplib)
	add_dependencies(chimp ${EXE}.x)
endforeach()

#-------------
#Kaiju: Gamera
message("Adding Gamera module ...")

#Get IC files
file(GLOB GAMICLIST src/gamera/ICs/*.F90)
file(GLOB GAMICDEF src/gamera/ICs/null.F90)
set(GAMIC ${GAMICDEF} CACHE FILEPATH "Gamera IC file")
set_property(CACHE GAMIC PROPERTY STRINGS ${GAMICLIST})

#Add source
add_subdirectory(src/gamera)

#Print gamera info
message("\tBricksize is ${bricksize}")
message("\tIC file is ${GAMIC}")
add_custom_target(gamera ALL)
message("\tAdding executable gamera.x")
add_executable(gamera.x src/drivers/gamerax.F90)
target_link_libraries(gamera.x baselib gamlib)
add_dependencies(gamera gamera.x)

#-------------
#Kaiju: ReMIX
message("Adding ReMIX module ...")
#Add source
add_subdirectory(src/remix)
add_custom_target(remix ALL)
message("\tAdding executable gamera2remix.x")
add_executable(gamera2remix.x src/drivers/gamera2remix.F90)
target_link_libraries(gamera2remix.x baselib interfaceslib gamlib remixlib)
add_dependencies(remix gamera2remix.x)
