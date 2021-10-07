#!/bin/sh

if [ ! -d "./codecov_prof" ] 
then
    echo "The directory codecov_prof does not exist. Please move to the folder where the test executables are, and be sure to run them first."
    exit
fi

profmerge -cov_dir ./codecov_prof -prof_dir ./codecov_prof -src_root . -verbose  # merge all code coverage files  
# -comp: exclude filenames that contain 'tst' using a compfile
# -ccolor: use green for fully covered instead of white
cd ./codecov_prof
codecov -comp ../../../tests/comp_file.txt -ccolor '#d7fad2'

