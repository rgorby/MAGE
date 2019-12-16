file(GLOB chimp_srcs *.F90)
add_library(chimplib ${chimp_srcs} ${EBIC} ${TPIC})
target_link_libraries(chimplib baselib)