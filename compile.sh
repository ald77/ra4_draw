#! /bin/bash

exit_code=0;

if [ $# -ne 0 ] && [ "$1" == "clean" ]
then
    rm -rf run/*.exe bin/*.o bin/*.a bin/*.d *.out src/baby*.cpp inc/baby*.hpp
    ./run/remove_backups.sh
    exit_code=$?
else
    bad_file=$(mktemp -t compile_XXXXXXXXXXXXXXXX)

    make -j $(getconf _NPROCESSORS_ONLN) -k -r -R 2> $bad_file
    exit_code=$?
    
    if [[ $exit_code != 0 ]] ; then
	printf "\n\e[31m==========  ERRORS AND WARNINGS  ============\e[0m\n"
	cat $bad_file >&2
    else
	printf "\n\e[32mCompiled successfully without errors or warnings!\e[0m\n\n"
    fi
    
    rm -rf $bad_file
fi

exit $exit_code
