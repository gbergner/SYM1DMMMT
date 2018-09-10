#!/bin/bash
# Very simple configure scpit, has to be extended.

pwd=$(pwd)
echo "Please select the build version of the code: 1) normal serial host code with gnu compiler; 2) device and host code using pgi compiler; 3) mpi matrix parallel version"
read t
if [ "x$t" == "x1" ]
then
cp Makefile.inc_in_gnu Makefile.inc
fi
if [ "x$t" == "x2" ]
then
cp Makefile.inc_in Makefile.inc
fi
if [ "x$t" == "x3" ]
then
cp Makefile.inc_in_mpi Makefile.inc
fi


sed -i "s&THESOURCE&${pwd}&g" Makefile.inc

echo "please check Makefile.inc for further options"
