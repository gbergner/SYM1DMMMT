# SYM1DMMMT

Current version is a development version of the code.

## Compile

The configuration is done with a small script ./configure.sh. This allows to select between different version of the code

1) serial version (one core CPU), forked from Masanori Hanada and updated 

2) GPU version (CUDA + OpenAcc)

3) matrix parallel (MPI), forked from Masanori Hanada and updated

There are still a view settings that have to be manually edited in the Makefile.inc

## Compile time parameters

In the current version of the code several parameters have to be set manually. These are in the file

staticparameters.f90 

The MPI code has a separate file

matrix_parallel/size_parallel.h

## Testruns

These can be found in the subfolders example*.
