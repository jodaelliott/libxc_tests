PROGRAM example_basicF90

   ! From Quantum ESPRESSO
   USE control_flags, ONLY : gamma_only
   USE io_files,      ONLY : prefix, restart_dir, tmp_dir
   USE io_rho_xml,    ONLY : read_scf
   USE lsda_mod,      ONLY : nspin
   USE mp_global,     ONLY : mp_startup
   USE fft_base,      ONLY : dfftp ! This is so we can allocate evc the same size as rho
   USE scf,           ONLY : rho, rho_core, rhog_core
   USE cell_base,     ONLY : omega
   USE gvect,         ONLY : g, ngm
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
                            xc_f90_lda_vxc,                        &
                            xc_f90_gga_exc,                        &
                            xc_f90_gga_vxc,                        &
                            XC_UNPOLARIZED,                        &
                            XC_LDA_X,                              &
                            XC_LDA_C_PZ,                           &
                            XC_GGA_X_PBE,                          &
                            XC_GGA_C_PBE,                          &
                            ! XC Family IDs
                            XC_FAMILY_LDA, XC_FAMILY_GGA,          &
                            XC_FAMILY_MGGA, XC_FAMILY_HYB_LDA,     &
                            XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA, &
                            ! XC Kind IDs
                            XC_EXCHANGE, XC_CORRELATION,           & 
                            XC_EXCHANGE_CORRELATION

   IMPLICIT NONE

   TYPE(xc_f90_func_t)      :: xc_func_x, xc_func_c  ! Functional type
   TYPE(xc_f90_func_info_t) :: xc_info_x, xc_info_c  ! Functional information

   REAL(8) :: exc(5)

   REAL(8), ALLOCATABLE :: ex(:,:), ec(:,:), sigma(:), vxc(:,:)
   REAL(8), ALLOCATABLE :: vrho_x(:,:), vrho_c(:,:), vsigma_x(:,:), vsigma_c(:,:)

   REAL(8) :: etx, etc, etxc, vtxc ! The total exchange energy

   COMPLEX(8), ALLOCATABLE :: rhogsum(:)
   REAL(8), ALLOCATABLE :: grho(:,:,:)

   CHARACTER(LEN=120) :: xc_name_x, xc_fam_x, xc_part_x
   CHARACTER(LEN=120) :: xc_name_c, xc_fam_c, xc_part_c
   CHARACTER(LEN=256) :: rho_file

   INTEGER :: i, k, vmajor, vminor, vmicro 
   INTEGER :: func_x_id = XC_GGA_X_PBE, &
              func_c_id = XC_GGA_C_PBE

   INTEGER(8) :: nnr

   CALL mp_startup( start_images=.TRUE., images_only=.TRUE.)

   prefix = 'ch4'
   tmp_dir = './example_sim/TMP/'

   CALL read_file()
   WRITE(rho_file, '(a)') TRIM(restart_dir()) // 'charge-density.dat'
   ! Next we have to get rho from the completed simulation
   CALL read_scf( rho, nspin, gamma_only)
   CALL set_rhoc()

   ! Allocate the exchange-correlation energy arrray the same way as rho
   ! in scf_mod.f90
   nnr = dfftp%nnr
   ALLOCATE(ex(dfftp%nnr, nspin))
   ALLOCATE(ec(dfftp%nnr, nspin))
   ALLOCATE(vrho_x(dfftp%nnr, nspin))
   ALLOCATE(vrho_c(dfftp%nnr, nspin))
   ALLOCATE(vsigma_x(dfftp%nnr, nspin))
   ALLOCATE(vsigma_c(dfftp%nnr, nspin))
   ALLOCATE(vxc(dfftp%nnr, nspin))

   ! First we call the function to print the version of 
   ! libxc we are using and write the output to terminal
   CALL xc_f90_version(vmajor, vminor, vmicro)
   WRITE(*,100) vmajor, vminor, vmicro

   !------------------------------------------------------
   ! We check if we need to calculate the gradient of rho
   ! by getting info on the functional chosen
   !------------------------------------------------------

   ! Let us get some information about the exchange 
   ! functional, and write that to the terminal as well
   ! ...Exchange
   CALL xc_f90_func_init(xc_func_x, func_x_id, XC_UNPOLARIZED)
   xc_info_x = xc_f90_func_get_info(xc_func_x)
   xc_name_x = xc_f90_func_info_get_name(xc_info_x)
   ! ...Correlation
   CALL xc_f90_func_init(xc_func_c, func_c_id, XC_UNPOLARIZED)
   xc_info_c = xc_f90_func_get_info(xc_func_c)
   xc_name_c = xc_f90_func_info_get_name(xc_info_c)

   ! To understand which family we have, we select on case by
   ! case basis
   ! ...Exchange
   SELECT CASE(xc_f90_func_info_get_family(xc_info_x))
      CASE (XC_FAMILY_LDA);      WRITE(xc_fam_x, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);      WRITE(xc_fam_x, '(a)') "GGA"
      CASE (XC_FAMILY_MGGA);     WRITE(xc_fam_x, '(a)') "Meta GGA"
      CASE (XC_FAMILY_HYB_LDA);  WRITE(xc_fam_x, '(a)') "Hybrid LDA"
      CASE (XC_FAMILY_HYB_GGA);  WRITE(xc_fam_x, '(a)') "Hybrid GGA"
      CASE (XC_FAMILY_HYB_MGGA); WRITE(xc_fam_x, '(a)') "Hybrid Meta GGA"
   END SELECT
   ! ...Correlation
   SELECT CASE(xc_f90_func_info_get_family(xc_info_c))
      CASE (XC_FAMILY_LDA);       WRITE(xc_fam_c, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);       WRITE(xc_fam_c, '(a)') "GGA"
      CASE (XC_FAMILY_MGGA);      WRITE(xc_fam_c, '(a)') "Meta GGA"
      CASE (XC_FAMILY_HYB_LDA);   WRITE(xc_fam_c, '(a)') "Hybrid LDA"
      CASE (XC_FAMILY_HYB_GGA);   WRITE(xc_fam_c, '(a)') "Hybrid GGA"
      CASE (XC_FAMILY_HYB_MGGA);  WRITE(xc_fam_c, '(a)') "Hybrid Meta GGA"
   END SELECT

   ! To understand which part of the functional we have we also
   ! select on case by case basis
   ! ...Exchange
   SELECT CASE(xc_f90_func_info_get_kind(xc_info_x))
      CASE (XC_EXCHANGE);             WRITE(xc_part_x, '(a)') "Exchange"
      CASE (XC_CORRELATION);          WRITE(xc_part_x, '(a)') "Correlation"
      CASE (XC_EXCHANGE_CORRELATION) 
         WRITE(xc_part_x, '(a)') "Exchange-Correlation"
   END SELECT
   ! ...Correlation
   SELECT CASE(xc_f90_func_info_get_kind(xc_info_c))
      CASE (XC_EXCHANGE);             WRITE(xc_part_c, '(a)') "Exchange"
      CASE (XC_CORRELATION);          WRITE(xc_part_c, '(a)') "Correlation"
      CASE (XC_EXCHANGE_CORRELATION) 
         WRITE(xc_part_c, '(a)') "Exchange-Correlation"
   END SELECT

   ! Write the functional information to terminal
   ! ...Exchange
   WRITE(*,101) TRIM(xc_name_x), TRIM(xc_fam_x), TRIM(xc_part_x)
   ! ...Correlation
   WRITE(*,101) TRIM(xc_name_c), TRIM(xc_fam_c), TRIM(xc_part_c)

   ! If the X or the C is GGA, then we have to compute the gradient of the density
   IF (xc_fam_x == 'GGA' .OR. xc_fam_c == 'GGA') THEN

      ALLOCATE(rhogsum(ngm))
      ALLOCATE(grho(3,dfftp%nnr,1))
      ALLOCATE(sigma(dfftp%nnr))
      
      DO k = 1, ngm
         rhogsum(k) = rho%of_g(k,1) ! ...This should work for now but needs updating 
      END DO

      CALL fft_gradient_g2r(dfftp, rhogsum, g, grho)

      DO i = 1, dfftp%nnr
         !sigma(i) = grho(1,i,1)*grho(1,i,1)
         sigma(i) = grho(1,i,1)**2 + grho(2,i,1)**2 + grho(3,i,1)**2
      END DO

      DEALLOCATE(rhogsum)
      DEALLOCATE(grho)
   END IF

   !------------------------------------------
   ! First we calculate the Exchange energy
   !------------------------------------------

   ! ...At this point we can apply the functional, using the family to
   ! ...decide how to select the appropriate way to compute
   SELECT CASE(xc_f90_func_info_get_family(xc_info_x))
      CASE (XC_FAMILY_LDA)
         CALL xc_f90_lda_exc(xc_func_x, nnr, rho%of_r(1,1), ex(1,1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(xc_func_x, nnr, rho%of_r(1,1), sigma(1), ex(1,1))
         CALL xc_f90_gga_vxc(xc_func_x, nnr, rho%of_r(1,1), sigma(1), vrho_x(1,1), vsigma_x(1,1))
   END SELECT

   etx = 0.

   !$OMP PARALLEL SHARED(etx)
   !$OMP DO REDUCTION(+:etx)
   DO i = 1, dfftp%nnr
      etx = etx + 2.*ex(i,1)*rho%of_r(i,1) 
   END DO
   !$OMP END DO
   !$OMP END PARALLEL

   WRITE(*,103) omega * etx/nnr

   CALL xc_f90_func_end(xc_func_x)

   !---------------------------------------------
   ! Second we calculate the Correlation energy
   !---------------------------------------------

   ! At this point we can apply the functional, using the family to
   ! decide how to select the appropriate way to compute
   SELECT CASE(xc_f90_func_info_get_family(xc_info_c))
      CASE (XC_FAMILY_LDA)
         CALL xc_f90_lda_exc(xc_func_c, nnr, rho%of_r(1,1), ec(1,1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(xc_func_c, nnr, rho%of_r(1,1), sigma(1),  ec(1,1))
         CALL xc_f90_gga_vxc(xc_func_c, nnr, rho%of_r(1,1), sigma(1),  vrho_c(1,1), vsigma_c(1,1))
   END SELECT

   etc = 0.

   !$OMP PARALLEL SHARED(etc)
   !$OMP DO REDUCTION(+:etc)
   DO i = 1, dfftp%nnr
      etc = etc + 2.* ec(i,1) * rho%of_r(i,1) 
   END DO
   !$OMP END DO
   !$OMP END PARALLEL

   WRITE(*,104) omega * etc/nnr

   !---------------------------------------------
   ! Combine exchange and correlation parts and
   ! apply any gradient corrections
   !---------------------------------------------
   
   vtxc = 0.
   etxc = 0.
   DO i = 1, dfftp%nnr
      vxc(i,1) = 2. * (vrho_x(i,1) + vrho_c(i,1))
      vtxc = vtxc + vxc(i,1)*rho%of_r(i,1)
   END DO   

   etxc = omega * (etc+etx)/nnr
   vtxc = omega * vtxc/nnr

!   etxc = 0.
!   WRITE(*,105) etxc
!   CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, vxc )   

   WRITE(*,105) etxc

   CALL xc_f90_func_end(xc_func_c)

   DEALLOCATE(ex, ec, vrho_x, vrho_c, vsigma_x, vsigma_c, vxc)

   CALL stop_run(0)
   CALL do_stop(0)

100 FORMAT("Libxc version: ",I1,"."I1"."I1)
101 FORMAT("           XC: ",a/"       Family: ",a/"         Kind: ",a)
102 FORMAT(2F16.7)
103 FORMAT("Integral of exchange energies   : ",F16.8)
104 FORMAT("Integral of correlation energies: ",F16.8)
105 FORMAT("Exchange correlation energy     : ",F16.8)

END PROGRAM
