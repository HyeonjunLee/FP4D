module readInput

  use FP4D_globals
 
  implicit none

contains

  subroutine readNamelistInput(filename)

  implicit none

  character(len=*), intent(in) :: filename
  integer :: myrank !temprary before mpi
  integer :: fileUnit, didFileAccessWork
  integer :: isp
  real (kind=8) :: l_ie

  namelist / opt_param / &
    opt_Col,opt_PS,opt_0M,opt_APS,opt_QLD,opt_NB,opt_Epar,opt_global, opt_restart, &
    opt_PS_aux, opt_QLD_aux, opt_NB_aux, opt_profile_rho, opt_QLD_dist, opt_time_solver, opt_time_species, &
    opt_Potential, &
    special_tau_W, special_tau_He3, &
    opt_rotation, omega_rot_in, smooth_tau, &
    opt_logaritm_NEO, opt_remove_homogeneous
  namelist / time_param / &
    sml_dt, sml_nstep, start_step_Col, step_Col, start_step_RHS, on_RHS, cycle_RHS
  namelist / phy_param / &
    sml_mass, sml_charge, sml_den, &
    sml_t_ev, sml_flow, sml_isp, sml_nsp, &
    dlnn_dr, dlnT_dr, sml_Epar
  namelist / APS_param / &
    sml_source_den, sml_source_t_ev, sml_source_width, APS_sp
  namelist / NB_param / &
    NB_sp, NB_pwrscale, NB_sink_opt, NB_sink_sp, start_Q, NB_fast_den_factor, NB_num
  namelist / QLD_param / &
    RF_kpar, RF_omega, RF_W0, &
    RF_width, RF_n_harmonic, mx_sp, &
    pwrscale, &
    RF_vr1, RF_vr2, RF_vz1, RF_vz2

  namelist / f_param / &
    f_nvp, f_nmu, f_vp_max, f_smu_max, &
    f_velo_aniso
  namelist / PS_param / &
    itmax,fw,fr,sup_step2,sup_val2
  namelist / eqb_param/ &
    sml_eqb_opt,&
    deltapsi
  namelist / col_f_param / &
    col_f_nthreads, neg_frac, den_converge, ene_converge, mom_converge, col_iteration
  namelist / out_param / &
    out_h5, out_step

  namelist / geo_param / &
    npsi,nth,ifd0,ifd1,op_geo_info,ipe,xx0,xx1,q0,q1,q2,q3,q4,&
    a0,R0,B0,ase_q0,ase_a0,ase_R0,ase_B0,ase_kappa,ase_del,ase_shift
     
  fileUnit=11
  open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
    print *,"Proc ",myRank,": Error opening ", filename
    stop
  else     
    read(fileUnit, nml=opt_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the opt_param namelist in it."
      stop
    end if


    read(fileUnit, nml=time_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the time_param namelist in it."
      stop
    end if

    read(fileUnit, nml=phy_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the phy_param namelist in it."
      stop
    end if

    read(fileUnit, nml=APS_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the APS_param namelist in it."
      stop
    end if

    read(fileUnit, nml=NB_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the NB_param namelist in it."
      stop
    end if

    read(fileUnit, nml=QLD_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the QLD_param namelist in it."
      stop
    end if
 
    read(fileUnit, nml=eqb_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the eqb_param namelist in it."
      stop
    end if
     
    read(fileUnit, nml=f_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the f_param namelist in it."
      stop
    end if

    read(fileUnit, nml=PS_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the PS_param namelist in it."
      stop
    end if
 
    read(fileUnit, nml=col_f_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the col_f_param namelist in it."
      stop
    end if

    read(fileUnit, nml=out_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the out_param namelist in it."
      stop
    end if

    read(fileUnit, nml=geo_param, iostat=didFileAccessWork)
    if (didFileAccessWork /= 0) then
      print *,"Proc ",myRank,": Error!  I was able to open the file ", filename, &
        " but not read data from the geo_param namelist in it."
      stop
    end if
  end if

  write(*,*) 'hell-youknow-R0', R0

  close(unit = fileUnit)

  ! Write out all parameters
  if(sml_mype==0) then
    print *,"Successfully read parameters from namelist in ", filename, "."
    open(unit=15,file='input.used',status='replace')
    write(15,nml=opt_param)
    write(15,nml=time_param)
    write(15,nml=phy_param)
    write(15,nml=APS_param)
    write(15,nml=NB_param)
    write(15,nml=QLD_param)
    write(15,nml=eqb_param)
    write(15,nml=f_param)
    write(15,nml=PS_param)
    write(15,nml=col_f_param)
    write(15,nml=out_param)
    write(15,nml=geo_param)
    close(15)
  endif


  ! Convert input variables from normalized to SI units
  do isp=0,sml_nsp_max
    sml_mass(isp) = sml_mass(isp)*sml_prot_mass
    sml_charge_num (isp) = sml_charge(isp)
    sml_charge(isp) = sml_charge(isp)*sml_e_charge
  enddo
  ! Compute electron collision time using the ion-electron Coulomb logarithm
  l_ie  = 24D0-log(sqrt(sml_den(0)*1D-6)/sml_t_ev(0))
!  sml_tau_e = sqrt(sml_t_ev(0)**3)/(2.91D-6*sml_den(0)*1D-6*l_ie)
  sml_tau = 3.606D26*sqrt(sml_mass(0)*sml_t_ev(0)**3)/(sml_den(0)*l_ie) 
  if (opt_time_solver.ne.-1) then
     ! if opt_time_solver==-1, each species has different time-step
     sml_dt = sml_dt*sml_tau
  endif

  ! ndstep = (sml_nstep -1)/out_step+1

  if (mod(sml_nstep,out_step) .ne. 0) then
    STOP "STOP due to wrong [out_step]"
  endif

  ! output time step -> 0:ndstep
  ndstep = (sml_nstep)/out_step

  if (sml_mype.eq.0) then
    write (*,'(a,(ES12.4E2),a,(ES12.4E2))') '  tau [s] :',sml_tau,'  sml_dt [s]:',sml_dt
  endif
  end subroutine readNamelistInput

end module readInput

