#!/bin/bash
# Very simple configure scpit, has to be extended.
pwd=$(pwd)
cp Makefile.inc_in Makefile.inc
sed -i "s&THESOURCE&${pwd}&g" Makefile.inc

