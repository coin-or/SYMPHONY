#!/bin/sh


if  [ "$1" = "make" ];
then
    make -C ../ install
elif [ "$1" = "gdb" ];
then
    gdb --args ./src/symphony -F Datasets/lecture.mps -f test/symopt
elif [ "$1" = "ill" ];
then
    ./src/symphony -F Datasets/illcond1.mps -f test/symopt 
elif [ "$1" = "other" ];
then
    ./src/symphony -F Datasets/"$2".mps -f test/symopt 
else
    ./src/symphony -BACH -F Datasets/lecture.mps -f test/symopt 
fi

