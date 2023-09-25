file(GLOB volt_srcs *.F90 gameraInterfaces/*.F90 remixInterfaces/*.F90 chimpInterfaces/*.F90 modelInterfaces/*.F90)
add_library(voltlib ${volt_srcs} ${VOLTIC})
target_link_libraries(voltlib rcmlib remixlib chimplib gamlib baselib)
