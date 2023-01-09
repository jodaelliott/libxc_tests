# Alternate notes for work on Mac M1

I have exported the following environment variables:

    export LD_QE71=/Users/josh/Applications/ESPRESSO/q-e-qe-7.1
    export LD_LIBXC=/Users/josh/Applications/LIBXC/libxc-6.0.0/GCC

Compiling with <code>mpif90</code>, linking openmp using <code>-fopenmp</code>

    mpif90 -fopenmp -o qe_read.x qe_read.f90 -L${LD_QE71}/dft-d3 -ldftd3qe -L${LD_QE71}/PW/src/ -lpw -L${LD_QE71}/Modules/ -lqemod -L${LD_QE71}/upflib/ -lupf ${LD_QE71}/XClib/xc_lib.a -L${LD_QE71}/FFTXlib/src/ -lqefft -L${LD_QE71}/KS_Solvers -lks_solvers -L${LD_QE71}/LAXlib/ -lqela -L${LD_QE71}/UtilXlib/ -lutil -L${LD_QE71}/MBD/ -lmbd -L${LD_QE71}/external/devxlib/src -ldevXlib -L/usr/local/lib -llapack  -lblas -L${LD_QE71}/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys -I${LD_QE71}/external/devxlib/src -I${LD_QE71}/include -I${LD_QE71}/FoX/finclude -I${LD_QE71}/upflib -I${LD_QE71}/XClib -I${LD_QE71}/Modules -I${LD_QE71}/FFTXlib -I${LD_QE71}/KS_Solvers -I${LD_QE71}/LAXlib -I${LD_QE71}/UtilXlib -I${LD_QE71}/MBD -I${LD_QE71}/FoX/finclude -I${LD_QE71}/dft-d3 -I${LD_QE71}/PW/src/ -I${LD_QE71}/dft-d3

## Additional calls:

On Mac, at runtime <code>qe_read.x</code> throws up an error about calls to mpi routines. This is because the MPI environment was not initialized at the start of the program. To fix the the following code was added

    USE mp_global, ONLY : mp_startup

    CALL mp_startup( start_images=.TRUE., images_only=.TRUE.)
