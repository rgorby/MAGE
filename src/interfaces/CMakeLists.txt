file(GLOB interfaces_srcs *.F90)
add_library(interfaceslib ${interfaces_srcs})
target_link_libraries(interfaceslib remixlib gamlib baselib)
