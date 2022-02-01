# RanBoost

   **RanBoost** is fortran code to implement Random Boosting in a genome-enabled prediction framework. 
   Full description of the algorith is available from the original [paper](https://doi.org/10.3168/jds.2012-5630).

# Introduction
This manual describes how to use the program BLasso, which is focused on the analysis of genomic data using Bayesian LASSO. BLasso can analyze continuos and categorical traits. 

The code is written in fortran, with GNU GPL license. The program is compiled to run in all kind of platforms (windows, linux, mac, ..).


# Purpose
This manual does not aim to thoroughly describe the methodology behind Bayesian LASSO, but to be a self-explanatory guide to implement BLasso in user's own data. Statistical details fro Bayesian Lasso can be found in [1] and [2]. The user is encouraged to consult them for details. This is a developing software, therefore any feedback on possible bugs, problems, running errors or suggestions are welcome and encouraged.

# 1 - How to execute BLasso

Download the executable file for your operating system:
MacOSX .
Linux .

It is also possible to clone the repository and compile the source code (e.g. with [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries)):
```
mkdir GWP
git clone https://github.com/ogrecio/RanBoost.git
cd RanBoost
gfortran src/RanBoost.f90 -O3 -o RanBoost.x

```

BLasso must be run in a command line mode. The execution is simple, just type in the command line the following:

```
./RanBoost 
```

RanBoost needs a parameter file named 'xxx.bl' 

>Caution:  
>  store different runs in different folders, or BLasso may overwrite previous results.
The program will prompt the residual variance and posterior mean of the intercept every 500 iterations.

