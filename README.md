# libxc_tests
Attemps to understand how to use libxc with FORTRAN


## Compilation of <code>libxc-6.0.0</code> with intel/2022

    wget https://gitlab.com/libxc/libxc/-/archive/6.0.0/libxc-6.0.0.tar.gz
    tar -zxvf libxc-6.0.0.tar.gz
    cd libxc-6.0.0/
    autoreconf -i
    
The step <code>autoreconf -i</code> generates the autotools configure script.

    export CC=icpc
    ./configure --prefix='/scratch/Applications/libxc/install-libxc-6.0.0_intel-2022'
    make
    make check
    make install
