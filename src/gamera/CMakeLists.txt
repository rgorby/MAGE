file(GLOB gam_srcs *.F90 apps/*.F90)
add_library(gamlib ${gam_srcs})
target_link_libraries(gamlib baselib)
