PROGRAM example_basicF90

   USE xc_f90_lib_m, ONLY : xc_f90_version,                        & 
                            xc_f90_func_t,                         &   
                            xc_f90_func_info_t,                    &
                            xc_f90_func_get_info,                  &
                            xc_f90_func_info_get_name,             &
                            xc_f90_func_info_get_family,           &
                            xc_f90_func_info_get_kind,             &
                            xc_f90_func_init,                      &
                            xc_f90_func_end,                       &
                            xc_f90_lda_exc,                        &
                            xc_f90_gga_exc,                        &
                            XC_UNPOLARIZED,                        &
                            XC_LDA_X,                              & ! LDA Perdew-Zunger
                            XC_LDA_C_PZ,                           & ! LDA Perdew-Zunger
                            ! XC Family IDs
                            XC_FAMILY_LDA, XC_FAMILY_GGA,          &
                            XC_FAMILY_MGGA, XC_FAMILY_HYB_LDA,     &
                            XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA, &
                            ! XC Kind IDs
                            XC_EXCHANGE, XC_CORRELATION,           & 
                            XC_EXCHANGE_CORRELATION

   IMPLICIT NONE

   TYPE(xc_f90_func_t)      :: x_func, c_func ! Functional type
   TYPE(xc_f90_func_info_t) :: x_info, c_info ! Functional information

   ! Dimension of following variables is the number of grid points
   REAL(8), DIMENSION(:), ALLOCATABLE :: rho,&      ! Density
                                         rho_core,& ! Core Density
                                         ex,&       ! Exchange energy on grid point
                                         ec,&       ! Correlation energy on grid point
                                         sigma      ! Something to do with GGA, not used yet 
   REAL(8)                            :: Etxc       ! Total exchange-correlation energy

   CHARACTER(LEN=120) :: x_name, x_fam, x_part,&
                         c_name, c_fam, c_part

   INTEGER :: i, vmajor, vminor, vmicro, func_id = 1
   INTEGER :: func_x_id = XC_LDA_X, &
              func_c_id = XC_LDA_C_PZ

   ! These are variables that should be retriveable from the 
   ! Quantum ESPRESSO simulation:
   INTEGER(8) :: nr1 = 64,&        ! FFT grid along each direction
                 nr2 = 64,&
                 nr3 = 64,&
                 fftgp             ! Total number of grid points
   REAL(8)    :: omega = 3375.0000 ! Cell volume
  
   fftgp = nr1*nr2*nr3

   ! Allocate arrays
   ALLOCATE(rho(fftgp))
   ALLOCATE(rho_core(fftgp))
   ALLOCATE(ec(fftgp))
   ALLOCATE(ex(fftgp))
   ALLOCATE(sigma(fftgp))

   ! Read valence and core densities and sum them to compute total rho
   OPEN(UNIT=999,FILE='rho.dat')
   OPEN(UNIT=998,FILE='rho_core.dat')

   DO i = 1, 262144
      READ(999,*) rho(i)
      READ(998,*) rho_core(i)
   END DO

   CLOSE(998)
   CLOSE(999)

   rho(:) = rho(:) + rho_core(:)

   ! First we call the function to print the version of 
   ! libxc we are using and write the output to terminal
   CALL xc_f90_version(vmajor, vminor, vmicro)
   WRITE(*,100) vmajor, vminor, vmicro

   ! Let us get some information about the exchange correlation
   ! functional, and write that to the terminal as well
   CALL xc_f90_func_init(x_func, func_x_id, XC_UNPOLARIZED)
   CALL xc_f90_func_init(c_func, func_c_id, XC_UNPOLARIZED)
   x_info = xc_f90_func_get_info(x_func)
   c_info = xc_f90_func_get_info(c_func)
   x_name = xc_f90_func_info_get_name(x_info)
   c_name = xc_f90_func_info_get_name(c_info)

   ! To understand which family we have, we select on case by
   ! case basis
   SELECT CASE(xc_f90_func_info_get_family(x_info))
      CASE (XC_FAMILY_LDA);      WRITE(x_fam, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);      WRITE(x_fam, '(a)') "GGA"
      CASE (XC_FAMILY_MGGA);     WRITE(x_fam, '(a)') "Meta GGA"
      CASE (XC_FAMILY_HYB_LDA);  WRITE(x_fam, '(a)') "Hybrid LDA"
      CASE (XC_FAMILY_HYB_GGA);  WRITE(x_fam, '(a)') "Hybrid GGA"
      CASE (XC_FAMILY_HYB_MGGA); WRITE(x_fam, '(a)') "Hybrid Meta GGA"
   END SELECT

   SELECT CASE(xc_f90_func_info_get_family(c_info))
      CASE (XC_FAMILY_LDA);      WRITE(c_fam, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);      WRITE(c_fam, '(a)') "GGA"
      CASE (XC_FAMILY_MGGA);     WRITE(c_fam, '(a)') "Meta GGA"
      CASE (XC_FAMILY_HYB_LDA);  WRITE(c_fam, '(a)') "Hybrid LDA"
      CASE (XC_FAMILY_HYB_GGA);  WRITE(c_fam, '(a)') "Hybrid GGA"
      CASE (XC_FAMILY_HYB_MGGA); WRITE(c_fam, '(a)') "Hybrid Meta GGA"
   END SELECT

   ! To understand which part of the functional we have we also
   ! select on case by case basis
   SELECT CASE(xc_f90_func_info_get_kind(x_info))
      CASE (XC_EXCHANGE);             WRITE(x_part, '(a)') "Exchange"
      CASE (XC_CORRELATION);          WRITE(x_part, '(a)') "Correlation"
      CASE (XC_EXCHANGE_CORRELATION) 
         WRITE(x_part, '(a)') "Exchange-Correlation"
   END SELECT

   SELECT CASE(xc_f90_func_info_get_kind(c_info))
      CASE (XC_EXCHANGE);             WRITE(c_part, '(a)') "Exchange"
      CASE (XC_CORRELATION);          WRITE(c_part, '(a)') "Correlation"
      CASE (XC_EXCHANGE_CORRELATION) 
         WRITE(c_part, '(a)') "Exchange-Correlation"
   END SELECT

   ! Write the functional information to terminal
   WRITE(*,101) TRIM(x_name), TRIM(x_fam), TRIM(x_part)
   WRITE(*,101) TRIM(c_name), TRIM(c_fam), TRIM(c_part)

   ! At this point we can apply the functional, using the family to
   ! decide how to select the appropriate way to compute
   ! The exchange part:
   SELECT CASE(xc_f90_func_info_get_family(x_info))
      CASE (XC_FAMILY_LDA)
         CALL xc_f90_lda_exc(x_func, fftgp, rho(1), ex(1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(x_func, fftgp, rho(1), sigma(1),  ex(1))
   END SELECT
   ! The correlation part:
   SELECT CASE(xc_f90_func_info_get_family(c_info))
      CASE (XC_FAMILY_LDA)
         CALL xc_f90_lda_exc(c_func, fftgp, rho(1), ec(1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(c_func, fftgp, rho(1), sigma(1),  ec(1))
   END SELECT

   Etxc = 0.
   DO i = 1, fftgp
     Etxc = Etxc + 2.* (ex(i) + ec(i)) * rho(i)
   END DO

!   WRITE(*,*) Etxc
   WRITE(*,103) omega * Etxc/fftgp

   CALL xc_f90_func_end(x_func)

   ! Deallocate arrays
   DEALLOCATE(rho, rho_core, ex, ec, sigma)

100 FORMAT("Libxc version: ",I1,"."I1"."I1)
101 FORMAT("           XC: ",a/"       Family: ",a/"         Kind: ",a)
102 FORMAT(2F16.7)
103 FORMAT(/"Computed Exchange-Correlation Energy: ",F10.6,/)
END PROGRAM
