#!/bin/bash

# make sure everything is up to date
git pull
git submodule update --init --recursive

# compile
scons --cfg=build/config-linux-gcc.py -j8
if [ $? -ne 0 ]
then
    echo Compilation failed.
    echo If the configuration failed, make sure that you installed all prerequisites.
    echo Have a look at install_prerequisites_ubuntu.sh
    echo
    echo Please report back any issues you encountered.
    exit 1
fi
