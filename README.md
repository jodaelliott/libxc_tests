# libxc_tests
Attemps to understand how to use libxc with FORTRAN


## Compilation of <code>libxc-6.0.0</code> with intel/2022

    wget https://gitlab.com/libxc/libxc/-/archive/6.0.0/libxc-6.0.0.tar.gz
    tar -zxvf libxc-6.0.0.tar.gz
    cd libxc-6.0.0/
    autoreconf -i
    
The step <code>autoreconf -i</code> generates the autotools configure script.

    export CC=gfortran
    export LD_LIBXC=/scratch/Applications/libxc/libxc/GCC
    ./configure --prefix=${LD_LIBXC}
    make
    make check
    make install

Code can be compiled, linking the library (as per https://stackoverflow.com/questions/74800483/undefined-reference-when-trying-to-link-libxc-to-fortran)

    gfortran -o example_basic.x example_basic.f90 -I${LD_LIBXC}/include -L${LD_LIBXC}/lib/ -lxcf90 -lxc 
