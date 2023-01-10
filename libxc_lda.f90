PROGRAM example_basicF90

   ! From Quantum ESPRESSO
   USE control_flags, ONLY : gamma_only
   USE io_files,      ONLY : prefix, restart_dir, tmp_dir
   USE io_rho_xml,    ONLY : read_scf
   USE lsda_mod,      ONLY : nspin
   USE mp_global,     ONLY : mp_startup
   USE fft_base,      ONLY : dfftp ! This is so we can allocate evc the same size as rho
   USE scf,           ONLY : rho
   USE cell_base,     ONLY : omega
   ! From libxc
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
                            XC_LDA_X,                              &
                            XC_LDA_C_PZ,                           &
                            ! XC Family IDs
                            XC_FAMILY_LDA, XC_FAMILY_GGA,          &
                            XC_FAMILY_MGGA, XC_FAMILY_HYB_LDA,     &
                            XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA, &
                            ! XC Kind IDs
                            XC_EXCHANGE, XC_CORRELATION,           & 
                            XC_EXCHANGE_CORRELATION

   IMPLICIT NONE

   TYPE(xc_f90_func_t)      :: xc_func_x, xc_func_c ! Functional type
   TYPE(xc_f90_func_info_t) :: xc_info              ! Functional information

   REAL(8) :: rho1(5)   = (/0.1, 0.2, 0.3, 0.4, 0.5/)
   REAL(8) :: sigma(5) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
   REAL(8) :: exc(5)

   REAL(8), ALLOCATABLE :: exc_espresso(:,:)

   REAL(8) :: etx, etc ! The total exchange energy

   CHARACTER(LEN=120) :: xc_name, xc_fam, xc_part
   CHARACTER(LEN=256) :: rho_file

   INTEGER :: i, vmajor, vminor, vmicro 
   INTEGER :: func_x_id = XC_LDA_X, &
              func_c_id = XC_LDA_C_PZ

   INTEGER(8) :: nnr

   CALL mp_startup( start_images=.TRUE., images_only=.TRUE.)

   prefix = 'ch4'
   tmp_dir = './example_sim/TMP/'

   CALL read_file()
   WRITE(rho_file, '(a)') TRIM(restart_dir()) // 'charge-density.dat'
   CALL read_scf( rho, nspin, gamma_only)

   ! Allocate the exchange-correlation energy arrray the same way as rho
   ! in scf_mod.f90
   nnr = dfftp%nnr
   ALLOCATE(exc_espresso(dfftp%nnr, nspin))

   ! First we call the function to print the version of 
   ! libxc we are using and write the output to terminal
   CALL xc_f90_version(vmajor, vminor, vmicro)
   WRITE(*,100) vmajor, vminor, vmicro

   ! Next we have to get rho from the completed simulation

   !------------------------------------------
   ! First we calculate the Exchange energy
   !------------------------------------------

   ! Let us get some information about the exchange 
   ! functional, and write that to the terminal as well
   CALL xc_f90_func_init(xc_func_x, func_x_id, XC_UNPOLARIZED)
   xc_info = xc_f90_func_get_info(xc_func_x)
   xc_name = xc_f90_func_info_get_name(xc_info)

   ! To understand which family we have, we select on case by
   ! case basis
   SELECT CASE(xc_f90_func_info_get_family(xc_info))
      CASE (XC_FAMILY_LDA);      WRITE(xc_fam, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);      WRITE(xc_fam, '(a)') "GGA"
      CASE (XC_FAMILY_MGGA);     WRITE(xc_fam, '(a)') "Meta GGA"
      CASE (XC_FAMILY_HYB_LDA);  WRITE(xc_fam, '(a)') "Hybrid LDA"
      CASE (XC_FAMILY_HYB_GGA);  WRITE(xc_fam, '(a)') "Hybrid GGA"
      CASE (XC_FAMILY_HYB_MGGA); WRITE(xc_fam, '(a)') "Hybrid Meta GGA"
   END SELECT

   ! To understand which part of the functional we have we also
   ! select on case by case basis
   SELECT CASE(xc_f90_func_info_get_kind(xc_info))
      CASE (XC_EXCHANGE);             WRITE(xc_part, '(a)') "Exchange"
      CASE (XC_CORRELATION);          WRITE(xc_part, '(a)') "Correlation"
      CASE (XC_EXCHANGE_CORRELATION) 
         WRITE(xc_part, '(a)') "Exchange-Correlation"
   END SELECT

   ! Write the functional information to terminal
   WRITE(*,101) TRIM(xc_name), TRIM(xc_fam), TRIM(xc_part)

   ! At this point we can apply the functional, using the family to
   ! decide how to select the appropriate way to compute
   SELECT CASE(xc_f90_func_info_get_family(xc_info))
      CASE (XC_FAMILY_LDA)
         CALL xc_f90_lda_exc(xc_func_x, nnr, rho%of_r(1,1), exc_espresso(1,1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(xc_func_x, 5_8, rho1(1), sigma(1),  exc(1))
   END SELECT

   etx = 0.

   !$OMP PARALLEL SHARED(etx)
   !$OMP DO REDUCTION(+:etx)
   DO i = 1, dfftp%nnr
      etx = etx + 2.*exc_espresso(i,1)*rho%of_r(i,1) 
   END DO
   !$OMP END DO
   !$OMP END PARALLEL

!   WRITE(*,103) SUM(exc_espresso)
   WRITE(*,103) omega * etx/nnr

   CALL xc_f90_func_end(xc_func_x)

   DEALLOCATE(exc_espresso)
   ALLOCATE(exc_espresso(dfftp%nnr, nspin))

   !---------------------------------------------
   ! Second we calculate the Correlation energy
   !---------------------------------------------

   ! This gets information about the correlation 
   ! functional, and writes in back to the terminal
   CALL xc_f90_func_init(xc_func_c, func_c_id, XC_UNPOLARIZED)
   xc_info = xc_f90_func_get_info(xc_func_c)
   xc_name = xc_f90_func_info_get_name(xc_info)

   SELECT CASE(xc_f90_func_info_get_family(xc_info))
      CASE (XC_FAMILY_LDA);       WRITE(xc_fam, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);       WRITE(xc_fam, '(a)') "GGA"
      CASE (XC_FAMILY_MGGA);      WRITE(xc_fam, '(a)') "Meta GGA"
      CASE (XC_FAMILY_HYB_LDA);   WRITE(xc_fam, '(a)') "Hybrid LDA"
      CASE (XC_FAMILY_HYB_GGA);   WRITE(xc_fam, '(a)') "Hybrid GGA"
      CASE (XC_FAMILY_HYB_MGGA);  WRITE(xc_fam, '(a)') "Hybrid Meta GGA"
   END SELECT

   SELECT CASE(xc_f90_func_info_get_kind(xc_info))
      CASE (XC_EXCHANGE);             WRITE(xc_part, '(a)') "Exchange"
      CASE (XC_CORRELATION);          WRITE(xc_part, '(a)') "Correlation"
      CASE (XC_EXCHANGE_CORRELATION) 
         WRITE(xc_part, '(a)') "Exchange-Correlation"
   END SELECT

   ! Write the functional information to terminal
   WRITE(*,101) TRIM(xc_name), TRIM(xc_fam), TRIM(xc_part)

   ! At this point we can apply the functional, using the family to
   ! decide how to select the appropriate way to compute
   SELECT CASE(xc_f90_func_info_get_family(xc_info))
      CASE (XC_FAMILY_LDA)
         CALL xc_f90_lda_exc(xc_func_c, nnr, rho%of_r(1,1), exc_espresso(1,1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(xc_func_c, nnr, rho1(1), sigma(1),  exc(1))
   END SELECT

   etc = 0.

   !$OMP PARALLEL SHARED(etc)
   !$OMP DO REDUCTION(+:etc)
   DO i = 1, dfftp%nnr
      etc = etc + 2.* exc_espresso(i,1) * rho%of_r(i,1) 
   END DO
   !$OMP END DO
   !$OMP END PARALLEL

   WRITE(*,104) omega * etc/nnr

   CALL xc_f90_func_end(xc_func_c)

   DEALLOCATE(exc_espresso)

   CALL stop_run(0)
   CALL do_stop(0)


100 FORMAT("Libxc version: ",I1,"."I1"."I1)
101 FORMAT("           XC: ",a/"       Family: ",a/"         Kind: ",a)
102 FORMAT(2F16.7)
103 FORMAT("Integral of exchange energies   : ",F16.8)
104 FORMAT("Integral of correlation energies: ",F16.8)
END PROGRAM
