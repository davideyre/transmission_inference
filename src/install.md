# Installation

Details are provided here for how to compile and develop using Xcode 10.1. 

## OpenMP support

Details on openMP support in Xcode can be found here - https://github.com/Homebrew/homebrew-core/issues/3266, see 2nd post

Summary of steps to install, assuming already have homebrew installed:

```
brew update
brew install llvm
```

Then in Xcode under Build Settings:

* Add a new user-defined setting CC with the value /usr/local/opt/llvm/bin/clang
* Add -fopenmp to Other C Flags
* Add /usr/local/opt/llvm/lib/clang/7.0.0/include to Header Search Paths
* Add /usr/local/opt/llvm/lib to Library Search Paths
* Set Enable Modules (C and Objective-C) to No.
* Set Index-While-Building to No

Under Build Phases:

* Add /usr/local/opt/llvm/lib/libiomp5.dylib to Link Binary With Libraries

Done. You can now #include <omp.h> and start using #pragma omp ... in your source code.

Can turn on and off with -fopenmp flag.

More details here - http://bisqwit.iki.fi/story/howto/openmp/ on openMP in general.


## RMath

To install standalone Rmath on mac OS, see section 2 and 9 in https://cran.r-project.org/doc/manuals/r-release/R-admin.html:

* download R source
* `tar -xvzf R-3.3.3.tar.gz`
* `rm R-3.3.3.tar.gz`
* in directory ./configure
* `cd src/nmath/standalone`
* `make`
* `make install`

Alternatively `brew install R`  seems to provide the necessary files.

May also need:

* brew install gcc
* brew install pcre

In Xcode:
 * Add `/usr/local/include` to header search path and `/usr/local/lib` to library search path
 * Add -lRMath to other linker flags

Install on linux:
* download R source, and in main directory untar the archiave
* `./configure --prefix=/home/local/GEL/davide/bin/R` or similar 
* `cd src/nmath/standalone/`
* `make; make install`


## Armadillo support

This is not used at present, but previous instructions are here for completeness:

On mac:
```
brew tap homebrew/science
brew install cmake pkg-config
brew install armadillo
```

then test:
`clang++ -I/usr/local/include -L/usr/local/lib test.cpp -o test -O2 -larmadillo -llapack -lblas`

for Xcode:
* add to Header Search Path: /usr/local/include
* add to Library Search Path: /usr/local/lib
* add to other linker flags: -larmadillo -llapack -lblas
 * set optimisation flag: Other C++ flags -O2

## Instructions for compiling outside of Xcode

### Compile - mac
```
cd /Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/src
g++ -o ../bin/transmission_test -std=c++11 -I /usr/local/include/ -L /usr/local/lib -lRmath -O2 *.cpp 
```

then test:
`./transmission_test -s 4598174 -i 5000 -p /Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/simulation_97674744`


### Compile on linux - ResComp

```
module purge
module load R/3.2.2
module load gcc/5.4.0
cd /users/bag/deyre/analysis/transmission_inference/src
g++ -o ../bin/transmission_test *.cpp -std=c++11 -lRmath -I /apps/well/R/3.2.2/lib64/R/include -L /apps/well/R/3.2.2/lib64 -lpthread -Wl,-rpath,/apps/well/R/3.2.2/lib64
```
last -Wl option required to ensure can find dynamic library for RMath

### Compile on linux - ARC

Load required modules:
```
module load R
module load gcc

cd /home/clme-transmission-modelling/gree0826/src
g++ -o /data/clme-transmission-modelling/gree0826/bin/transmission_test /data/clme-transmission-modelling/gree0826/src/*.cpp -std=c++11 -lRmath -I /system/software/linux-x86_64/R/3.3/include -L /system/software/linux-x86_64/R/3.3/lib64
```
run with:
```
/home/clme-transmission-modelling/gree0826/bin/transmission_test -s 4598174 -i 5000 -p /home/clme-transmission-modelling/gree0826/simulation/simulation_716365 > /home/clme-transmission-modelling/gree0826/log/simulation_716365.log 2>&1 &
```

