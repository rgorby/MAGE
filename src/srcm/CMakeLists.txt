file(GLOB srcmlib_srcs *.F90)
add_library(srcmlib ${srcmlib_srcs})
target_link_libraries(srcmlib baselib)
