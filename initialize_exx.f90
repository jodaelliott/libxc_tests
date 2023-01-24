PROGRAM initialize_exx

! Code which shows how to initialize the elements needed for an exact exchange
! calculation from a previous SCF run with pw.x
! Initial version of the code mirrors the skeleton procedure of a calculation='nscf' run
! with pw.x

   USE environment,          ONLY : environment_start
   USE read_input,           ONLY : read_input_file
   USE command_line_options, ONLY : input_file_
   USE check_stop,           ONLY : check_stop_init, check_stop_now

   USE klist,                ONLY : nks, nkstot
   USE ener,                 ONLY : ef, ef_up, ef_dw

   USE input_parameters,     ONLY : outdir

   USE control_flags, ONLY : gamma_only, lbands
   USE io_files,      ONLY : prefix, restart_dir, tmp_dir
   USE io_rho_xml,    ONLY : read_scf
   USE lsda_mod,      ONLY : nspin
   USE mp_global,     ONLY : mp_startup
   USE scf,           ONLY : rho, rho_core
   USE fft_base,      ONLY : dfftp
   USE exx,           ONLY : exxinit, use_ace, aceinit
   USE xc_lib,        ONLY : stop_exx
   USE wvfct,         ONLY : nbnd, et
   USE basis,          ONLY : starting_wfc, starting_pot

   IMPLICIT NONE

   REAL(8)            :: ef_scf, ef_scf_up, ef_scf_dw
   CHARACTER(LEN=256) :: rho_file
   INTEGER            :: grid_point

   !
   ! ... Init. parallelization, (N.B. only npools is really relevant here)
   !
   CALL mp_startup( )
   !
   ! ... Start clocks, set code and version, clear previous CRASH and print info on run.
   !
   CALL environment_start('EXX Test')
   !
   ! ... Read input data (currently nscf run - this should be replaced with a read_file)
   !
   CALL read_input_file( 'PW', input_file_)
   !
   ! ...At this point we can set/update any input keywords using the input_parameters module.
   !
   starting_wfc = 'file'
   starting_pot = 'file'
   prefix = 'ch4'
   outdir = './example_sim/TMP/'
   !use_ace = .TRUE. 

   !
   ! ... Copy input parameters in input file format to internal format
   !
   CALL iosys()
   !
   ! ... Create an instance to allow the code to terminate when it has done its job.
   !
   CALL check_stop_init()
   !
   ! ... Important routine that needs to be understood. This calls all sorts of things also
   ! ... called by XSpectra including DFT+U(+V) stuff and divide_et_impera, as well as 
   ! ... setup_exx which has rules that it must be called after setup_para, but before init_run.
   ! ... understanding this is probably going to really help in implementation in XSpectra.
   !
   CALL setup()
   CALL init_run()

   CALL c_bands_nscf()

   CALL poolrecover(et, nbnd, nkstot, nks)
   CALL sum_band()

   ef_scf = ef
   ef_scf_up = ef_up
   ef_scf_dw = ef_dw

   IF (lbands) THEN
      CALL weights_only()
   ELSE
      CALL weights()
   END IF

   !
   ! ... Make sure exx_is_active is .FALSE. to carry out all required initializations
   !
   CALL stop_exx()
   !
   ! ... Initialize everything needed for the EXX calculation
   !  
   CALL exxinit( .FALSE., nbnd)
   !
   ! ... Initialize everything needed for the ACE operators
   !
   CALL aceinit( .FALSE. )


END PROGRAM
