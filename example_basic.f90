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
                            XC_GGA_C_PBE,                          &
                            ! XC Family IDs
                            XC_FAMILY_LDA, XC_FAMILY_GGA,          &
                            XC_FAMILY_MGGA, XC_FAMILY_HYB_LDA,     &
                            XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA, &
                            ! XC Kind IDs
                            XC_EXCHANGE, XC_CORRELATION,           & 
                            XC_EXCHANGE_CORRELATION

   IMPLICIT NONE

   TYPE(xc_f90_func_t)      :: xc_func ! Functional type
   TYPE(xc_f90_func_info_t) :: xc_info ! Functional information

   REAL(8) :: rho(5)   = (/0.1, 0.2, 0.3, 0.4, 0.5/)
   REAL(8) :: sigma(5) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
   REAL(8) :: exc(5)

   CHARACTER(LEN=120) :: xc_name, xc_fam, xc_part

   INTEGER :: i, vmajor, vminor, vmicro, func_id = XC_GGA_C_PBE

   ! First we call the function to print the version of 
   ! libxc we are using and write the output to terminal
   CALL xc_f90_version(vmajor, vminor, vmicro)
   WRITE(*,100) vmajor, vminor, vmicro

   ! Let us get some information about the exchange correlation
   ! functional, and write that to the terminal as well
   CALL xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
   xc_info = xc_f90_func_get_info(xc_func)
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
         CALL xc_f90_lda_exc(xc_func, 5_8, rho(1), exc(1))
      CASE (XC_FAMILY_GGA)
         CALL xc_f90_gga_exc(xc_func, 5_8, rho(1), sigma(1),  exc(1))
   END SELECT

   ! Write the results to the terminal
   WRITE(*,'(a)') "# Rho Exc"
   DO i = 1, 5
      WRITE(*,102) rho(i), exc(i)
   END DO

   CALL xc_f90_func_end(xc_func)

100 FORMAT("Libxc version: ",I1,"."I1"."I1)
101 FORMAT("           XC: ",a/"       Family: ",a/"         Kind: ",a)
102 FORMAT(2F16.7)
END PROGRAM
