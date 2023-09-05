file(GLOB gl_srcs *.F90 apps/*.F90)
add_library(gllib ${gl_srcs})
target_link_libraries(gllib baselib)
