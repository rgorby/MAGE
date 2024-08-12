file(GLOB volt_mpi_srcs *.F90 modelInterfaces/*.F90)
add_library(voltmpilib ${volt_mpi_srcs})
target_link_libraries(voltmpilib voltlib rcmlib remixlib chimplib gammpilib gamlib basempilib baselib)
