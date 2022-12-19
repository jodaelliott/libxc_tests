# libxc_tests
Attemps to understand how to use libxc with FORTRAN


## Compilation of <code>libxc-6.0.0</code>

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

    gfortran -o qe_read qe_read.f90 /scratch/Applications/q-e_gcc_9.4.0/PW/src/libpw.a /scratch/Applications/q-e_gcc_9.4.0//Modules/libqemod.a /scratch/Applications/q-e_gcc_9.4.0//upflib/libupf.a /scratch/Applications/q-e_gcc_9.4.0//XClib/xc_lib.a /scratch/Applications/q-e_gcc_9.4.0//FFTXlib/libqefft.a -L/scratch/Applications/q-e_gcc_9.4.0//KS_Solvers -lks_solvers /scratch/Applications/q-e_gcc_9.4.0//LAXlib/libqela.a /scratch/Applications/q-e_gcc_9.4.0//UtilXlib/libutil.a /scratch/Applications/q-e_gcc_9.4.0//MBD/libmbd.a  -L/scratch/Applications/q-e_gcc_9.4.0//external/devxlib/src -ldevXlib  -L/dls_sw/apps/gcc/9.2.0/7/lib64 -llapack  -lblas  -L/scratch/Applications/q-e_gcc_9.4.0//FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys  -lfftw3  -lblas -I/scratch/Applications/q-e_gcc_9.4.0//external/devxlib/src -I/scratch/Applications/q-e_gcc_9.4.0//include -I/scratch/Applications/q-e_gcc_9.4.0//FoX/finclude  -I/scratch/Applications/q-e_gcc_9.4.0//upflib -I/scratch/Applications/q-e_gcc_9.4.0//XClib -I/scratch/Applications/q-e_gcc_9.4.0//Modules -I/scratch/Applications/q-e_gcc_9.4.0//FFTXlib -I/scratch/Applications/q-e_gcc_9.4.0//KS_Solvers -I/scratch/Applications/q-e_gcc_9.4.0//LAXlib -I/scratch/Applications/q-e_gcc_9.4.0//UtilXlib -I/scratch/Applications/q-e_gcc_9.4.0//MBD -I/scratch/Applications/q-e_gcc_9.4.0//FoX/finclude -I/scratch/Applications/q-e_gcc_9.4.0/dft-d3 -I/scratch/Applications/q-e_gcc_9.4.0/PW/src/
