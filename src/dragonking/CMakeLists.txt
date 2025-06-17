file(GLOB dragonking_srcs *.F90)
add_library(dragonkinglib ${dragonking_srcs})
target_link_libraries(dragonkinglib remixlib baselib)
