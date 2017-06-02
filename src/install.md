Installation
============

OpenMP
------

Details on openMP support in Xcode here - https://github.com/Homebrew/homebrew-core/issues/3266, see 2nd post

Summary of steps to install:

brew update
brew install llvm

Then in Xcode under Build Settings:

* Add a new user-defined setting CC with the value /usr/local/opt/llvm/bin/clang
* Add -fopenmp to Other C Flags
* Add /usr/local/opt/llvm/lib/clang/3.9.1/include to Header Search Paths
* Add /usr/local/opt/llvm/lib to Library Search Paths
* Set Enable Modules (C and Objective-C) to No.

Under Build Phases:

* Add /usr/local/opt/llvm/lib/libiomp5.dylib to Link Binary With Libraries

Done. You can now #include <omp.h> and start using #pragma omp ... in your source code.

Can turn on and off with -fopenmp flag

More details here - http://bisqwit.iki.fi/story/howto/openmp/ on openMP in general


RMath
-----

To install standalone Rmath on mac OS, see section 2 and 9 in https://cran.r-project.org/doc/manuals/r-release/R-admin.html:

* download R source
* tar -xvzf R-3.3.3.tar.gz
* rm R-3.3.3.tar.gz
* in directory ./configure
* cd src/nmath/standalone
* make
* make install

May also need:

* brew install gcc
* brew install pcre

In Xcode may need to add /usr/local/include to header search path and /usr/local/lib to library search path


Install on linux:

* download R source, and in main directory tar -xf then
* ./configure --prefix=/home/local/GEL/davide/bin/R
* cd src/nmath/standalone/
* make; make install


Compile - mac
-------------
cd /Users/davideyre/Dropbox/Epi_WGS_MCMC/cpp_inference_scratchpad/inference_test/inference_test
g++ -o transmission_test -std=c++11 -I /usr/local/include/ -L /usr/local/lib -lRmath *.cpp

then test:
./transmission_test -s 4598174 -i 5000 -p /Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/simulation_97674744

Compile on linux - GEL
----------------------

On GEL cannot install Rmath to standard location, therefore
g++ -o test *.cpp -std=c++11 -I /home/local/GEL/davide/bin/R/include -L /home/local/GEL/davide/bin/R/lib -Wl,-rpath,/home/local/GEL/davide/bin/R/lib -lRmath -lpthread

Compile on linux - ResComp
----------------------

module purge
module load R/3.3.2

#files in `R RHOME`/include and `R RHOME`/lib, e.g.


On GEL cannot install Rmath to standard location, therefore
module load gcc/5.4.0
cd /well/bag/deyre/analysis/transmission_inference/src
g++ -o ../transmission *.cpp -std=c++11 -lRmath -I /apps/well/R/3.3.2/lib64/R/include -L /apps/well/R/3.3.2/lib64/R/lib 

module purge
module load R
g++ -o ../transmission *.cpp -std=c++11 -lRmath -I /apps/well/R/3.2.5-openblas-0.2.18-omp-gcc4.8.2/lib64/R/include -L /apps/well/R/3.2.5-openblas-0.2.18-omp-gcc4.8.2/lib64/R/lib 


g++ -o ../transmission *.cpp -std=c++11 -I /apps/well/R/3.2.5-openblas-0.2.18-omp-gcc4.8.2/lib64/R/include -L /apps/well/R/3.2.5-openblas-0.2.18-omp-gcc4.8.2/lib64/R/lib -Wl,-rpath,/apps/well/R/3.2.5-openblas-0.2.18-omp-gcc4.8.2/lib64/R/lib -lRmath -lpthread


Compile on linux - ARC
----------------------

Load required modules
module load R
module load gcc

/home/clme-transmission-modelling/gree0826/src

g++ -o /data/clme-transmission-modelling/gree0826/bin/transmission_test /data/clme-transmission-modelling/gree0826/src/*.cpp -std=c++11 -lRmath -I /system/software/linux-x86_64/R/3.3/include -L /system/software/linux-x86_64/R/3.3/lib64

run with

/home/clme-transmission-modelling/gree0826/bin/transmission_test -s 4598174 -i 5000 -p /home/clme-transmission-modelling/gree0826/simulation/simulation_716365 > /home/clme-transmission-modelling/gree0826/log/simulation_716365.log 2>&1 &


