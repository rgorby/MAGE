file(GLOB volt_srcs *.F90 gameraInterfaces/*.F90 remixInterfaces/*.F90)
add_library(voltlib ${volt_srcs})
target_link_libraries(voltlib remixlib gamlib baselib)
