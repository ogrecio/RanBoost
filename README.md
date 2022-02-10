# RanBoost

   **RanBoost** is fortran code to implement Random Boosting in a genome-enabled prediction framework. 
   Full description of the algorith is available from the original [paper](https://doi.org/10.3168/jds.2012-5630).

# Introduction
This manual describes how to use the program RanBoost-L2Boost, which is focused on the analysis of genomic data using Random Boosting. L2Boost can analyze continuos and categorical traits. 

The code is written in fortran, with GNU GPL license. The program is compiled to run in all kind of platforms (windows, linux, mac, ..).


# Purpose
The user is encouraged to consult statistical details in the original paper. This is a developing software, therefore any feedback on possible bugs, problems, running errors or suggestions are welcome and encouraged.

# 1 - How to execute RandomBoosting

Download the executable file for your operating system:
MacOSX .
Linux .

It is also possible to clone the repository and compile the source code (e.g. with [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries)):
```
mkdir GWP
cd GWP
git clone https://github.com/ogrecio/RanBoost.git
cd RanBoost
gfortran src/L2Boost.f90 -o L2Boost.x
```

BLasso must be run in a command line mode. The execution is simple, just type in the command line the following:

```
./L2Boost.x 
```

RanBoost needs a parameter file named 'parBOOST' 


