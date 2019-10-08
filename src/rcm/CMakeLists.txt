file(GLOB rcm_srcs *.f90 *.F)
add_library(rcmlib ${rcm_srcs})
target_link_libraries(rcmlib baselib)
