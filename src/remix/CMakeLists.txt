file(GLOB remix_srcs *.F90)
add_library(remixlib ${remix_srcs})
target_link_libraries(remixlib baselib)
