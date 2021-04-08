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
#Set minimum compiler versions
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 17.0)
		message("Fortran compiler too old!  What, were you gonna use punch cards?")
		message(FATAL_ERROR "ifort > 17.0 required")
	elseif( (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 17.0) AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 18.0) )
		message(WARNING "Compiler has incomplete F2008 features, Git hash/compiler information won't be saved to H5 files")
		add_compile_definitions(__INTEL_COMPILER_OLD)
	endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.0)
		message("Fortran compiler too old!  What, were you gonna use punch cards?")
		message(FATAL_ERROR "gfortran > 8.0 required")
	endif()
endif()

#-------------
#Set base release options
set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g")

#Do compiler specific options
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	set(dialect "-free -implicitnone")
	#Base
	string(APPEND CMAKE_Fortran_FLAGS " -fPIC")
	#Production
	set(PROD "-align array64byte -align rec32byte -no-prec-div -fast-transcendentals")
	#Debug
	set(DEBUG "-g -traceback -check bounds -check uninit -debug all -gen-interfaces -warn interfaces -fp-stack-check")
	set(PRODWITHDEBUGINFO "-O3 -g -traceback -debug all -align array64byte -align rec32byte -no-prec-div -fast-transcendentals")

	#Now do OS-dep options
	if (CMAKE_SYSTEM_NAME MATCHES Darwin)
		string(APPEND CMAKE_Fortran_FLAGS " -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 -xHost")
	else()
		#If we're not doing Mac, then add IPO
		string(APPEND PROD " -ipo")
		string(APPEND PRODWITHDEBUGINFO " -ipo")
	endif()

	#Handle individual hosts
	if (HOST MATCHES cheyenne)
		string(APPEND PROD " -march=corei7 -axCORE-AVX2")
		string(APPEND PRODWITHDEBUGINFO " -march=corei7 -axCORE-AVX2")
	endif()

	#Check Intel Fortran version
	if(NOT ALLOW_INVALID_COMPILERS AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "19")
		message(FATAL_ERROR "Intel Fortran compilers 19 or newer are not supported. Set the ALLOW_INVALID_COMPILERS variable to ON to force compilation at your own risk.")
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
string(APPEND CMAKE_Fortran_FLAGS_RELWITHDEBINFO " ${PRODWITHDEBUGINFO}")

if(ENABLE_OMP)
	string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
endif()
if(ENABLE_MPI)
	add_definitions(${MPI_Fortran_COMPILE_FLAGS})
	include_directories(${MPI_Fortran_INCLUDE_PATH})
	link_directories(${MPI_Fortran_LIBRARIES})
	if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
		string(APPEND CMAKE_Fortran_FLAGS " -mt_mpi")
	endif()
	# no matching flag for GNU
	set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

