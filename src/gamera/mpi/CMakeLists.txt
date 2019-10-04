file(GLOB gam_mpi_srcs *.F90)
add_library(gammpilib ${gam_mpi_srcs})
target_link_libraries(gammpilib basempilib baselib gamlib)
