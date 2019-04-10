file(GLOB base_srcs *.F90 model_defs/*.F90 model_types/*.F90)
add_library(baselib ${base_srcs})