file(GLOB base_mpi_srcs *.F90)
add_library(basempilib ${base_mpi_srcs})
target_link_libraries(basempilib baselib)

