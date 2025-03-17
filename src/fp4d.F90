!  Main program 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!          4-Dimension Relativistic Landau Fokker-Planck code
!                      for Heating & Current Drive           
!
!  FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444        4444 DDDDDDDDDDDD
!  FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444        4444 DDDDDDDDDDDDD
!  FFFF             PPPP        PPPP 4444        4444 DDDD      DDDD
!  FFFF             PPPP        PPPP 4444        4444 DDDD       DDDD
!  FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444444444444444 DDDD        DDDD
!  FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444444444444444 DDDD        DDDD
!  FFFF             PPPP                         4444 DDDD       DDDD
!  FFFF             PPPP                         4444 DDDD      DDDD
!  FFFF             PPPP                         4444 DDDDDDDDDDDDD
!  FFFF             PPPP                         4444 DDDDDDDDDDDD
!
!  This version is alpha 1.1                     Date: Mar. 17. 2025
!  
!  Contact: Hyeonjun   Lee (hyeonjun@hanyang.ac.kr)
!              Yunho Jeong (   yunho@hanyang.ac.kr)
!            Jungpyo   Lee ( jungpyo@hanyang.ac.kr)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                     |  Library Dependency  |
!----------------------------------------------------------------------
!  Name                   Version                    
!----------------------------------------------------------------------
!  Intel OneAPI    |          ?
!  HDF5            |       1.12.2
!  NetCDF          |       4. 4.5
!
!                  FP4D was developed with the libraries listed above!                                    
!----------------------------------------------------------------------
!----------------------------------------------------------------------    
program FP4D
  ! parallelization
  use omp_module
  use mpi
  use my_mpi_module
  use blacs_mod

  !read input & write output
  use readInput
  use readHDF5Input
  use readNetCDFInput
  use writeHDF5Output

  ! init & post
  use FP4D_globals
  use FP4D_init
  use FP4D_timer_lib
  use FP4D_Post
  use FP4D_equilibrium
  use f_module

  ! explicit
  use RHS_module
  use quasi_f_module

  ! implicit
  use col_f_module
  use imp_RHS

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none
  integer :: isp, unit, ierr, j
  integer :: loss_ok, ips, idt=0

  integer :: is, ip, ith, ir, iz

  character(len=20) :: timestamp
  real(kind=8) :: wtime, smooth

  logical :: RHS_update, Potent_update, Col_update, hdf5_update, Remove_homogeneous_update

  ! Initialize MPI environment
  call my_mpi_init

  ! Determine the number of OpenMP threads
  sml_nthreads = 1
#ifdef _OPENMP
  !$omp parallel
  !$omp master
  sml_nthreads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel
#endif

  if (sml_mype==0) then
    print*,'----------------------------------------------------------------------'
    print*,'----------------------------------------------------------------------'
    print*,'      4-Dimension Relativistic Landau Fokker-Planck code'
    print*,'                 for Heating and Current Drive          '
    print*,''
    print*,'FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444        4444 DDDDDDDDDDDD'
    print*,'FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444        4444 DDDDDDDDDDDDD'
    print*,'FFFF             PPPP        PPPP 4444        4444 DDDD      DDDD'
    print*,'FFFF             PPPP        PPPP 4444        4444 DDDD       DDDD'
    print*,'FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444444444444444 DDDD       DDDD'
    print*,'FFFFFFFFFFFFFFFF PPPPPPPPPPPPPPPP 4444444444444444 DDDD       DDDD'
    print*,'FFFF             PPPP                         4444 DDDD       DDDD'
    print*,'FFFF  Alpha 1.1  PPPP      Hanyang Univ.      4444 DDDD      DDDD'
    print*,'FFFF             PPPP  Nuclear Fusion and     4444 DDDDDDDDDDDDD'
    print*,'FFFF             PPPP  Plasma Computation LAB 4444 DDDDDDDDDDDD'
    print*,'----------------------------------------------------------------------'
    print*,'----------------------------------------------------------------------'
  endif

#ifdef _OPENACC
  if (sml_mype==0) then
    print*,'OpenACC recognized'
  endif
#endif


  !-----------------------------------------------------------------------
  ! START INIT.
  !
  ! INITIALIZATION ORDER
  !    :: readnamelist -> check input -> setup globalvars -> eq_init -> f_init -> the others
  
  ! make 'out.fp4d.run'
  call FP4D_runfile_create()

  ! read 'input.namelist'
  call readNamelistInput("input.namelist")

  ! check parameter consistency
  call FP4D_check()

  ! set the mpi grid
  call get_mpi_grid() !JPL  

  ! set the global index
  call setup_FP4D_globalvars()

  ! Initilize the geometry using eqdsk and profile
  call eq_init()

  ! Initialize the distribution function
  call f_init()

  ! Initialize NB source
  if (opt_NB.eq.1) then
    if (opt_restart.eq.0) then
      call setup_hdf(NB_num)
      call openInputFile('output.h5', '/BEAM_DEPO_INFO')
      call readHDF5()
      call closeInputFile()
    endif
  endif

  ! Initialize Q(f) from TORIC
  if (opt_QLD.eq.2) then
    call openNCDF()
  endif

  ! Initilize the partiacle source in & out
  call stsl_init(sml_totalpe,sml_nthreads,sml_mype,&
                 sml_source_den,sml_source_t_ev,sml_source_width)  ! steady NB, RF, alpha source

  ! Initialize the collision module & QLD & quasi-neutrality module

  ! f_dsmu = f_smu_max/real(f_nmu,8)
  ! f_dvp  = f_vp_max/real(f_nvp,8)
  ! f_dth  = sml_2pi/real(nth,8)

  write(*,*) 'col var check',  sml_isp,sml_nsp,npsi,nth,f_nmu,f_nvp,sml_dt,f_vp_max,f_smu_max,pwrscale,sml_mype

  call col_f_setup (is_s,is_e,npsi,nth, f_nmu,f_nvp, sml_dt, f_vp_max,f_smu_max, pwrscale,sml_mype)
  call impRHS_setup(is_s,is_e,npsi,nth, f_nmu,f_nvp, sml_dt, f_vp_max,f_smu_max, sml_mype)

  call RHS_f_setup (sml_mype)

  if ( (opt_Potential.ne.0) .and. (opt_Potential.ne.-1) ) then
    ! if 'opt_Potential' is not zero,
    ! automatically determine whether adiabatic electron case or not.
    call quasi_f_setup(sml_mype)
  endif
  ! END INIT.
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! output handle : dataxx.h5, out.fp4d.FSAx_xx
  !
  call allocate_time_var()

  ! remove existing data**.h5 if opt_restart == 0
  if (sml_mype==0) then      
    if (opt_restart.eq.0) then
      call system('rm -f data*.h5')
      print *, "Removed existing data**.h5 file"
    endif
  endif

  ! load latest data**.h5 if opt_restart != 0
  if (opt_restart.ne.0) then
    call checkPrevData(hdf5_name)
  endif

  ! open the out.fp4d.run & write time
  call FP4D_output_open()
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! START Main loop: Time stepping
  !
  call mpi_barrier(sml_comm,ierr)

  call FP4D_timer_init()

  call timer_lib_in('TOTAL')
  idt = 0

  do sml_istep=0,sml_nstep

    sml_time = sml_time + sml_dt

    !-----------------------------------------------------------------------
    ! determine whether to upate or not
    !
    ! RHS_update
    if (sml_istep .lt. start_step_RHS) then
      RHS_update = .true.
    else
      RHS_update = mod(sml_istep, (step_Col*cycle_RHS)) .le. (step_Col*on_RHS)
    endif

    ! Col_update
    if (sml_istep .eq. start_step_Col) then
      call col_dt_gt(step_Col)

      ! YOUKNOW ******************************************************
      call QLD_dt_gt(step_Col)
      ! YOUKNOW ******************************************************

      write(*,*) 'end col_dt_gt ' 
    endif

    if (sml_istep .lt. start_step_Col) then
      Col_update = .true.
    else
      Col_update = mod(sml_istep, step_Col) .eq. 0
    endif

    ! hdf5_update
    hdf5_update = mod(sml_istep, out_step) .eq. 0

    ! Potent_update :: follow the RHS in current...
    Potent_update = RHS_update

    ! Remove_homogeneous_update :: follow the Col_update
    
    Remove_homogeneous_update = Col_update
    !-----------------------------------------------------------------------


    ! ! TEST START

    ! write(*,*) 'it:',sml_istep, RHS_update, Col_update, hdf5_update

    ! RHS_update = .false.
    ! Col_update = .false.
    ! hdf5_update = .false.
    ! Potent_update = .false.

    ! ! TEST END

    !-----------------------------------------------------------------------
    ! START time-advance
    !
    ! RHS update

    call timer_lib_in('RHS')

    ! smooth turn of the drift term
    if (smooth_tau.eq.0) then
      smooth = 1.0
    else
      if (sml_istep .lt. smooth_tau*4) then
        smooth = 1.0 - exp(-real(sml_istep, 8)/smooth_tau)
      else
        smooth = 1.0
      endif
    endif
    if (Col_update) then
      write(*,*) 'smooth', smooth
    endif

    if (RHS_update) then

      call cal_PS (f0_M, f1, g1, dfdt_p)
      call cal_RF (f0_M, f1,     dfdt_q)
      k1 = (dfdt_p + dfdt_d*smooth + df_sts*sml_dt)
      k1 = k1 + MERGE (dfdt_q, 0D0, sml_istep .ge. start_Q) ! start_Q (default : 1)

      select case (opt_time_solver)
        ! explicit time solver
        ! (0) : default : Euler 
        ! (1) : RK45
        ! (2) : RK45 + update potentail each stage
        case(0)
          f1 = f1 + k1
          Ff = f0_M + f1

        case(1,-1)
          call cal_PS (f0_M, f1+k1/2.0, g1+k1/2.0, dfdt_p)
          call cal_RF (f0_M, f1+k1/2.0,            dfdt_q)
          k2 = (dfdt_p + dfdt_d*smooth + df_sts*sml_dt)
          k2 = k2 + MERGE (dfdt_q, 0D0, sml_istep .ge. start_Q)

          call cal_PS (f0_M, f1+k2/2.0, g1+k2/2.0, dfdt_p)
          call cal_RF (f0_M, f1+k2/2.0,            dfdt_q)
          k3 = (dfdt_p + dfdt_d*smooth + df_sts*sml_dt)
          k3 = k3 + MERGE (dfdt_q, 0D0, sml_istep .ge. start_Q)

          call cal_PS (f0_M, f1+k3, g1+k3, dfdt_p)
          call cal_RF (f0_M, f1+k3,        dfdt_q)
          k4 = (dfdt_p + dfdt_d*smooth + df_sts*sml_dt)
          k4 = k4 + MERGE (dfdt_q, 0D0, sml_istep .ge. start_Q)

          f1 = f1 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
          Ff = f0_M + f1

        case default
          write(*,*) 'error :: opt_RK45'
      end select

      ! sink for background species
      if ((opt_APS.eq.1).or.(opt_NB.eq.1)) then
        f1 = f1 - df_stl * sml_dt
        Ff = f0_M + f1
      endif

    end if ! RHS_update

    call timer_lib_out('RHS')



    ! Col update
    call timer_lib_in('Col')
    
    if (Col_update) then
    
      ! output : dfdt_c
      call fm_collision(Ff,g1,dfdt_c,sml_nthreads,&
                        sml_istep,&
                        num_n,num_T,num_p,FSAn,FSAT,FSAp,&
                        qld_B,qld_C,qld_F,f0_M_tar,col_power,col_dens,col_inner_info)

      
      if (opt_QLD .eq. 3) then
        call fm_QLD(Ff,g1,dfdt_q,sml_nthreads,&                 
                  qld_B,qld_C,qld_F,f0_M_tar)
      endif

      f1 = f1 + dfdt_c
      Ff = f0_M + f1
    else
      dfdt_c = 0D0
    endif

    call timer_lib_out('Col')

    ! opt_remove_homogeneous
    call timer_lib_in('Remove_homogeneous')

    if (Remove_homogeneous_update) then
      if (opt_remove_homogeneous .eq. 1) then
        f1_tmp = f1
        ! remove maxwellian
        call f1_Maxwellian(f1_tmp, f1)
        Ff = f0_M + f1
      endif
    endif

    call timer_lib_out('Remove_homogeneous')

    ! if (sml_istep.eq.sup_step2) then
    !   call impRHS_dt_gt(sup_val2)
    ! endif

    ! if (opt_Epar.eq.1) then
    !   if (sml_istep.lt.sup_step2) then
    !     call cal_impRHS(Ff,g1,dfdt_e,sml_nthreads,eq_den,eq_t_ev,sml_mass,sml_charge,sml_Epar,&
    !           eq_J,eq_invJ,eq_intJ,eq_B,eq_dBdth,eq_dphidth,sml_istep)
    !   elseif (sml_istep.ge.sup_step2) then
    !     if ((mod(sml_istep,sup_val2)).eq.1) then
    !       call cal_impRHS(Ff,g1,dfdt_e,sml_nthreads,eq_den,eq_t_ev,sml_mass,sml_charge,sml_Epar,&
    !         eq_J,eq_invJ,eq_intJ,eq_B,eq_dBdth,eq_dphidth,sml_istep)
    !     endif
    !   endif
    ! endif

    ! Potent update
    call timer_lib_in('Potent')

    if (Potent_update) then
      if (opt_Potential.eq.0) then
        g1 = f1
      elseif (opt_Potential.eq.-1) then
        g1 = f1
        call find_potential (g1, f0_M, f1_tmp, eq_phi1, eq_dphidth) ! 'f1_tmp' has trash val.
            ! originally this subroutine calcualte the phi & f1, from g1
      else
        ! input   : g1, f1, f0_M, eq_phi1_tmp
        ! output  : g1_k4, eq_phi1, eq_dphidth
        eq_phi1_tmp = eq_phi1 ! important for overlap the eq_phi1

        call fm_quasi(g1, g1_k4, f1, f0_M, sml_nthreads,&
                    eq_phi1_tmp, eq_phi1, eq_dphidth)
        g1 = g1_k4
      endif
    endif

    call timer_lib_out('Potent')
    !
    ! END time-advance
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! record distribution to save in hdf5
    if (hdf5_update) then
      ! RECORD TIME
      if (sml_mype.eq.0) then
        call cal_timestamp(timestamp)
        write(*,'(A, I0, A)') 'Time Step ',sml_istep,' calculated! @ ' // trim(timestamp)
      endif

      if (out_h5.eq.1) then
        ! idt = 0 ... ndstep
        call record_time_data (idt)
        idt=idt+1 ! idt = 1 ... ndstep + 1
      endif
    endif
    !-----------------------------------------------------------------------

  enddo ! do sml_istep=0,sml_nstep

  call mpi_barrier(sml_comm,ierr)

  call timer_lib_out('TOTAL')
  !
  ! END Main loop: Time stepping
  !-----------------------------------------------------------------------

  ! save in hdf5
  call timer_lib_in('hdf5')

  if (out_h5.eq.1) then
    call print_msg('write data to the HDF5 file')
    call openOutputFile(hdf5_name)
    call writeHDF5t (idt)
    call closeOutputFile()
  endif

  call timer_lib_out('hdf5')

  call mpi_barrier(sml_comm,ierr)

  ! Post-Processing
  call timer_lib_in('Post')

  if (sml_mype==0) then
    call print_msg('Post-Processing')
    call Post_f_setup(sml_mype)
    call vel_integral()
    call Post_record()
    call Post_finalize()
  endif

  call timer_lib_out('Post')

  !-----------------------------------------------------------------------
  ! Finalization
  !
  ! Close output files & time record
  call FP4D_output_close()

  ! Deallocate collision memory
  call col_f_finalize
  ! Deallocate distribution function memory
  call f_finalize
  ! Deallocate 
  call quasi_f_finalize
  ! Deallocate 
  call imp_finalize
   
  ! Print execution time
  if (sml_mype==0) then
    wtime = timer_lib_time('TOTAL')
    write(*,'(a,i20,a)') 'Execution time: ', int(wtime), ' s'
  endif

  call my_mpi_finalize()

end program FP4D
