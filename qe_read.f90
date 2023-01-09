PROGRAM qe_read

! To be able to compute the XC energy functional
! we will need to provide the density (rho) and the 
! gradient of the density (sigma) to libxc. 

! For simplicity of the demonstration this code will 
! read these quantities from a completed calculation.

   USE control_flags, ONLY : gamma_only
   USE io_files,      ONLY : prefix, restart_dir, tmp_dir
   USE io_rho_xml,    ONLY : read_scf
   USE lsda_mod,      ONLY : nspin
   USE mp_global,     ONLY : mp_startup
   USE scf,           ONLY : rho, rho_core
   USE fft_base,      ONLY : dfftp

   IMPLICIT NONE

   CHARACTER(LEN=256) :: rho_file
   INTEGER            :: grid_point

   CALL mp_startup( start_images=.TRUE., images_only=.TRUE. )

   prefix = 'ch4'
   tmp_dir = './example_sim/TMP/'

   ! Grab information from the completed calculation using read_file,
   ! this is under PW/src and so libpw needs to be linked correctly
   CALL read_file()   
   CALL set_rhoc()
   WRITE(*,'(a)') 'I Read the completed calculation'

   ! Read the charge density stored in charge-density.dat
   ! rho is TYPE(scf_type), print array dimensions to see what
   ! has been allocated
   WRITE(rho_file, '(a)') TRIM(restart_dir()) // 'charge-density.dat'
   CALL read_scf( rho, nspin, gamma_only)

   ! Write the valence and core density to file
   OPEN(UNIT=999, FILE='rho.dat')
   OPEN(UNIT=998, FILE='rho_core.dat')

   DO grid_point = 1, dfftp%nnr
      WRITE(999,*) rho%of_r(grid_point,1)
      WRITE(998,*) rho_core(grid_point)
   END DO

   CLOSE(998)
   CLOSE(999)

END PROGRAM
