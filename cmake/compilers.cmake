#Handle compiler/defaults

#-------------
#Handle libraries (HDF5/OMP/MPI)
find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
set(CMAKE_REQUIRED_INCLUDES ${HDF5_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
#Use set compiler (below) or link_libraries but not both
set(CMAKE_Fortran_COMPILER ${HDF5_Fortran_COMPILER_EXECUTABLE})

find_package(OpenMP COMPONENTS Fortran)
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_FLAGS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})

if (ENABLE_MPI)
	find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

#-------------
#Set default build to release
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

#-------------
#Set base release options
set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

#Do compiler specific options
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	set(dialect "-free -implicitnone")
	#Base
	string(APPEND CMAKE_Fortran_FLAGS " -fPIC")
	#Production
	set(PROD "-align array64byte -align rec32byte -no-prec-div -fast-transcendentals")
	#Debug
	set(DEBUG "-g -traceback -check bounds -check uninit -debug all -gen-interfaces -warn interfaces -fp-stack-check")

	#Now do OS-dep options
	if (CMAKE_SYSTEM_NAME MATCHES Darwin)
		string(APPEND CMAKE_Fortran_FLAGS " -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 -xHost")
	else()
		#If we're not doing Mac, then add IPO
		string(APPEND PROD " -ipo")
	endif()

	#Handle individual hosts
	if (HOST MATCHES cheyenne)
		string(APPEND PROD " -march=corei7 -axCORE-AVX2")
	endif()
	
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	set(dialect "-ffree-form -ffree-line-length-none -fimplicit-none")
	#Base
	string(APPEND CMAKE_Fortran_FLAGS " -fPIC")
	#Production
	set(PROD "-ffast-math")
	#Debug
	set(DEBUG "-fbacktrace -g -Warray-temporaries -Wall -Wfatal-errors  -finit-local-zero")
	#Now do machine-dep options
	if (CMAKE_SYSTEM_NAME MATCHES Darwin)
		string(APPEND CMAKE_Fortran_FLAGS " -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000")
	endif()
endif()

string(APPEND CMAKE_Fortran_FLAGS " ${dialect}")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " ${DEBUG}")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " ${PROD}")

if(ENABLE_OMP)
	string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
endif()
if(ENABLE_MPI)
	add_definitions(${MPI_Fortran_COMPILE_FLAGS})
	include_directories(${MPI_Fortran_INCLUDE_PATH})
	link_directories(${MPI_Fortran_LIBRARIES})
	string(APPEND CMAKE_Fortran_FLAGS " -mt_mpi")
	set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

