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
   ! USE read_file
   USE scf,           ONLY : rho

   IMPLICIT NONE

   CHARACTER(LEN=256) :: rho_file

   prefix = 'ch4'
   tmp_dir = './example_sim/TMP/'

   ! Grab information from the completed calculation using read_file,
   ! this is under PW/src and so libpw needs to be linked correctly
   CALL read_file()   

   ! Read the charge density stored in charge-density.dat
   ! rho is TYPE(scf_type), print array dimensions to see what
   ! has been allocated
   WRITE(rho_file, '(a)') TRIM(restart_dir()) // 'charge-density.dat'
   CALL read_scf( rho, nspin, gamma_only)

   WRITE(*,*) 'Rho in Real-Space: ', SIZE(rho%of_r), SUM(rho%of_r)
   WRITE(*,*) 'Rho in G-Space: ', SIZE(rho%of_g), SUM(rho%of_g)

END PROGRAM
