file(GLOB kaimaglib_srcs *.F90)
add_library(kaimaglib ${kaimaglib_srcs})
target_link_libraries(kaimaglib baselib)