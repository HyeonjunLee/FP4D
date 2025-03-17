#define HAVE_PARALLEL_HDF5

module writeHDF5Output

  use FP4D_globals
  use HDF5
  use mpi

  implicit none

  integer, private :: HDF5Error, ierr
  integer(HID_T), private :: HDF5FileID, parallelID
  integer(HID_T), dimension(:), allocatable, private :: groupIDs

  interface writeVariable
    module procedure writeVariable_integer
    module procedure writeVariable_scalar
  !  module procedure writeVariable_real
    module procedure writeVariable_1d
    module procedure writeVariable_2d
    module procedure writeVariable_3d
    module procedure writeVariable_4d
    module procedure writeVariable_5d
    module procedure writeVariable_6d
    module procedure writeVariable_7d
  end interface writeVariable

  interface rank
  !   module procedure rank_real
    module procedure rank_integer
    module procedure rank_character
    module procedure rank_scalar
    module procedure rank_1d
    module procedure rank_1d_integer
    module procedure rank_2d
    module procedure rank_3d
    module procedure rank_4d
    module procedure rank_5d
    module procedure rank_6d
  end interface rank

   real(kind=8), dimension(:,:,:,:,:,:), allocatable ::       &
     Ff_t, dfdt_col, dfdt_Q_t, dfdt_par, g1_t, &
     f1_t, dfdt_Epar, dfdt_impar
   real(kind=8), dimension(:,:,:), allocatable  :: phi1t, phi1t_k1, phi1t_k2, phi1t_k3
   real(kind=8), dimension(:,:,:,:), allocatable  :: den_t, ene_t, mom_t
   real(kind=8), dimension(:,:,:), allocatable  :: FSAn_t, FSAT_t, FSAp_t
      
   real(kind=8), dimension(:,:,:,:,:), allocatable  :: col_power_t ! youknow
   real(kind=8), dimension(:,:,:,:,:), allocatable  :: col_dens_t ! youknow
   real(kind=8), dimension(:,:,:,:,:,:,:), allocatable  :: col_inner_info_t ! youknow

   real(kind=8), dimension(:), allocatable :: time
contains

! -----------------------------------------------------------------------------------

  subroutine allocate_time_var()

    implicit none

    allocate(Ff_t       (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &
             g1_t       (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &
             f1_t       (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &

             dfdt_col   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &
             dfdt_Epar  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &
             dfdt_Q_t   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &
             dfdt_par   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &
             dfdt_impar (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep), &

             phi1t     (ith_s:ith_e, ip_s:ip_e, 0:ndstep),&
             phi1t_k1  (ith_s:ith_e, ip_s:ip_e, 0:ndstep),&
             phi1t_k2  (ith_s:ith_e, ip_s:ip_e, 0:ndstep),&
             phi1t_k3  (ith_s:ith_e, ip_s:ip_e, 0:ndstep),&
             den_t(ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep),&
             ene_t(ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep),&
             mom_t(ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep),&
             FSAn_t(ip_s:ip_e, is_s:is_e, 0:ndstep),&
             FSAT_t(ip_s:ip_e, is_s:is_e, 0:ndstep),&
             FSAp_t(ip_s:ip_e, is_s:is_e, 0:ndstep),&
             ! col_power_t (ip_s:ip_e, ith_s:ith_e, is_s:is_e, f_nsp-0+1, 0:ndstep), & ! youknow
             ! col_power_t (ith_s:ith_e,ip_s:ip_e,  is_s:is_e, f_nsp-0+1,0:ndstep), & ! JPL reorder
             col_power_t (ith_s:ith_e,ip_s:ip_e,  is_s:is_e, tar_s:tar_e,0:ndstep), & ! JPL reorder
             col_dens_t (ith_s:ith_e,ip_s:ip_e,  is_s:is_e, tar_s:tar_e,0:ndstep), & ! JPL reorder
             col_inner_info_t (ith_s:ith_e, ip_s:ip_e,  is_s:is_e, tar_s:tar_e, 0:2, 1:col_iteration, 0:ndstep), & ! JPL reorder
     
             time(0:ndstep))
    Ff_t=0D0; g1_t=0D0; dfdt_col=0D0;
    f1_t=0D0; dfdt_Epar=0D0; dfdt_Q_t=0D0
    dfdt_par=0D0; dfdt_impar=0d0;
    phi1t=0D0;
    phi1t_k1=0D0; phi1t_k2=0D0; phi1t_k3=0D0; 
    col_power_t=0D0;
    col_dens_t=0D0;
    col_inner_info_t=0D0;
    den_t=0D0; ene_t=0D0; mom_t=0D0
    FSAn_t=00;FSAT_t=0.0;FSAp_t=0.0

  end subroutine allocate_time_var

! -----------------------------------------------------------------------------------

  subroutine record_time_data(it)
    implicit none

    integer, intent(in) :: it

    ! Time Depedent Variables : Each Time Step!
    Ff_t(:,:,:,:,:,it) = Ff(:,:,:,:,:)
    g1_t(:,:,:,:,:,it) = g1(:,:,:,:,:)
    f1_t(:,:,:,:,:,it) = f1(:,:,:,:,:)
    dfdt_col(:,:,:,:,:,it)  = dfdt_c(:,:,:,:,:)
    dfdt_Epar(:,:,:,:,:,it) = dfdt_e(:,:,:,:,:)
    dfdt_Q_t(:,:,:,:,:,it) = dfdt_q(:,:,:,:,:)
    dfdt_par(:,:,:,:,:,it) = dfdt_p(:,:,:,:,:)
    dfdt_impar(:,:,:,:,:,it) = dfdt_imp(:,:,:,:,:)
    phi1t   (:,:,it) = eq_phi1   (:,:)
    phi1t_k1(:,:,it) = eq_phi1_k1(:,:)
    phi1t_k2(:,:,it) = eq_phi1_k2(:,:)
    phi1t_k3(:,:,it) = eq_phi1_k3(:,:)
    den_t(:,:,:,it) = num_n(:,:,:)
    ene_t(:,:,:,it) = num_T(:,:,:)
    mom_t(:,:,:,it) = num_p(:,:,:)
    FSAn_t(:,:,it) = FSAn(:,:)
    FSAT_t(:,:,it) = FSAT(:,:)
    FSAp_t(:,:,it) = FSAp(:,:)
    col_power_t(:,:,:,:,it) = col_power(:,:,:,:)
    col_dens_t(:,:,:,:,it) = col_dens(:,:,:,:)
    col_inner_info_t(:,:,:,:,:,:,it) = col_inner_info(:,:,:,:,:,:)
    time(it)=sml_time

  end subroutine record_time_data

! -----------------------------------------------------------------------------------

  subroutine openOutputFile(outputFilename)

    implicit none
    character (len=*), intent(in) :: outputFilename

    ! write(*,*) 'start openOutputFile', sml_mype

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call h5open_f(HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error initializing HDF5."
      stop
    end if
    ! print *, "HDF5 h5open_f status (Rank ", sml_mype, "): ", HDF5Error

    ! Initialize some stuff related to parallel file access
    call h5pcreate_f(H5P_FILE_ACCESS_F, parallelID, HDF5Error)
    ! print *, "HDF5 h5pcreate_f status (Rank ", sml_mype, "): ", HDF5Error
    
    call h5pset_fapl_mpio_f(parallelID, MPI_COMM_WORLD, MPI_INFO_NULL, HDF5Error)
    ! print *, "HDF5 h5pset_fapl_mpio_f status (Rank ", sml_mype, "): ", HDF5Error

    call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error, access_prp=parallelID)
    ! for serial version
    !  call h5fcreate_f(outputFilename, H5F_ACC_TRUNC_F, HDF5FileID, HDF5Error)

    ! print *, "HDF5 h5fcreate_f status (Rank ", sml_mype, "): ", HDF5Error

    if (HDF5Error < 0) then
      print *,"Error opening HDF5 output file."
      stop
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! write(*,*) 'end openOutputFile', sml_mype


  end subroutine openOutputFile
! -----------------------------------------------------------------------------------
  subroutine writeHDF5t(it)
    implicit none

    integer, intent(in) :: it

    write(*,*) 'start writeHDF5t', sml_mype

    !Time Independent Variables : First Time Step!

    call writeVariable_scalar(f_smu_max,"vperp_max",it)
    call writeVariable_scalar(f_vp_max,"vpar_max",it)
    call writeVariable_scalar(f_velo_aniso,"velo_aniso",it)
    call writeVariable_integer(sml_isp,"ispecies",it)
    call writeVariable_integer(sml_nsp,"nspecies",it)
    call writeVariable_integer(f_nmu,"Nvperp",it)
    call writeVariable_integer(f_nvp,"Nvpar",it)
    call writeVariable_scalar(sml_tau,"tau",it)
    call writeVariable_integer(sml_nstep,"nstep",it)
    call writeVariable_scalar(sml_prot_mass,"proton_mass",it)
    call writeVariable_scalar(sml_e_charge,"electron_charge",it)
    call writeVariable_scalar(sml_pi,"pi",it)
    call writeVariable_scalar(sml_2pi,"2pi",it)
    call writeVariable_scalar(sml_ev2j,"eV2J",it)
    call writeVariable_integer(sml_nsp_max,"nsp_max",it)
    call writeVariable_integer(step_Col,"step_Col",it)
    call writeVariable_integer(sup_val2,"sup_val2",it)
    call writeVariable_integer(start_step_Col,"start_step_Col",it)
    call writeVariable_integer(sup_step2,"sup_step2",it)
    call writeVariable_scalar(sml_j2ev,"J2eV",it)
    call writeVariable_scalar(sml_dt,"dt",it)
    call writeVariable_scalar(eq_dtheta,"dth",it)
    call writeVariable_integer(sml_nstep,"time_step",it)
    call writeVariable_scalar(sml_source_width,"source_width",it)
    call writeVariable_integer(NB_sp,"psource",it)
    call writeVariable_integer(APS_sp,"NBsource",it)
    call writeVariable_scalar(neg_frac,"neg_frac",it)
    call writeVariable_real(sml_den,"sml_n",it)
    call writeVariable_real(sml_t_ev,"sml_t_ev",it)
    call writeVariable_real(sml_flow,"flow_init",it)
    call writeVariable_real(sml_mass,"mass",it)
    call writeVariable_real(sml_charge,"charge",it)
    call writeVariable_real(dlnn_dr,"dlnn_dr",it)
    call writeVariable_real(dlnT_dr,"dlnT_dr",it)
    call writeVariable(eq_t_ev,"eq_t_ev",it)
    call writeVariable(eq_den,"eq_den",it)
    call writeVariable(eq_lambda,"eq_lambda",it)
    call writeVariable(eq_tau,"eq_tau",it)
    ! call writeVariable(eq_nu,"eq_nu",it)
    call writeVariable_real(sml_source_den,"source_density",it)
    call writeVariable_real(sml_source_t_ev,"source_temperature",it)
    call writeVariable_real(RF_W0,"W0",it)
    call writeVariable_real(RF_omega,"RF_frequency",it)
    call writeVariable_real(RF_width,"RF_width",it)
    call writeVariable_real(RF_kpar,"kpar",it)
    call writeVariable_integer(nth,"ntheta",it)
    call writeVariable_integer(npsi,"npsi",it)
    call writeVariable_real(eq_theta,"theta",it)
    call writeVariable_2D(eq_J,"J",it)
    call writeVariable_real(eq_psin,"psi",it) 
    call writeVariable_2D(eq_Bt,"Bt",it)
    call writeVariable_2D(eq_Bp,"Bp",it)
    call writeVariable_2D(eq_B,"B",it)
    call writeVariable_2D(eq_I,"I",it)
    call writeVariable_2D(eq_dBdth,"dBdth",it)
    call writeVariable_2D(eq_dBdpsi,"dBdpsi",it)
    call writeVariable_2D(eq_invJ,"invJ",it)
    call writeVariable_2D(eq_dlnn,"eq_dlnn",it)
    call writeVariable_2D(eq_dlnT,"eq_dlnT",it)
    call writeVariable(eq_omega,"omega",it)
    call writeVariable_5D(dfdt_d,"dfdt_d",it)
    call writeVariable(df_sts,"df_sts",it)
    call writeVariable(df_stl,"df_stl",it)
    call writeVariable_real(sml_Epar,"Epar",it)
    call writeVariable(vth,"vth",it)
    call writeVariable(eq_Volume,"eq_Volume",it) !youknow
    call writeVariable(vel_volume,"vel_volume",it) !youknow
    call writeVariable(ndstep,"ndstep",it)
    call writeVariable(out_step,"out_step",it)
    call writeVariable(f0_M,"f0_M",it)
    call writeVariable(f0_M_tar,"f0_M_tar",it)
    call writeVariable(geo_R,"geo_R",it)      !ym
    call writeVariable(geo_Z,"geo_Z",it)      !ym
    call writeVariable(geo_psi,"geo_psi",it)  !ym
    call writeVariable(geo_dRdth,"geo_dRdth",it)      !ym
    call writeVariable(sts_int,"sts_int",it)
    ! rotation
    if (opt_rotation .ne. 0) then
      call writeVariable(omega_rot, "omega_rot",  it)
      call writeVariable(phi_rot,   "phi_rot",    it)
      call writeVariable(phi_rot_avg,   "phi_rot_avg",    it)

      call writeVariable(phi_rot_star,   "phi_rot_star",    it)
      call writeVariable(phi_rot_star_avg,   "phi_rot_star_avg",    it)

      call writeVariable(lambda_rot,      "lambda_rot",     it)
      call writeVariable(dlambda_rot_dth,"dlambda_rot_dth", it)
      
      call writeVariable(den_rot,   "den_rot",    it)

      call writeVariable(small_n0,   "small_n0",    it)

    endif

    ! Time Depedent Variables : Each Time Step!

    call writeVariable_6d(Ff_t,"Ff",it)
    call writeVariable_6d(dfdt_col,"dfdt_col",it)
    call writeVariable_6d(dfdt_Epar,"dfdt_Epar",it)
    call writeVariable_6d(dfdt_Q_t,"dfdt_Q",it)
    call writeVariable_6d(dfdt_par,"dfdt_par",it)
    call writeVariable_6d(dfdt_impar,"dfdt_impar",it)    
    call writeVariable_6d(g1_t,"g1",it)
    call writeVariable_6d(f1_t,"f1",it)
    call writeVariable(phi1t,"phi1",it)
    call writeVariable(phi1t_k1,"phi1_k1",it)
    call writeVariable(phi1t_k2,"phi1_k2",it)
    call writeVariable(phi1t_k3,"phi1_k3",it)

    call writeVariable_5D(col_power_t,"col_power",it)
    call writeVariable_5D(col_dens_t,"col_dens",it)

    call writeVariable_7D(col_inner_info_t,"col_inner_info",it)

    call writeVariable_4d(den_t,"den_t",it)
    call writeVariable_4d(ene_t,"ene_t",it)
    call writeVariable_4d(mom_t,"mom_t",it)
    call writeVariable(FSAn_t,"FSAn_t",it)
    call writeVariable(FSAp_t,"FSAp_t",it)
    call writeVariable(FSAT_t,"FSAT_t",it)
    
    call writeVariable(time,"time",it)
    call writeVariable(sml_time,"end_time",it)
    
    ! MPI_BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    write(*,*) 'end writeHDF5t', sml_mype

  end subroutine writeHDF5t
! -----------------------------------------------------------------------------------

  subroutine closeOutputFile()

    implicit none

    print *, "HDF5  complete - start"

    call h5fclose_f(HDF5FileID, HDF5Error)
    print *, "HDF5 h5fclose_f status (Rank ", sml_mype, "): ", HDF5Error

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call h5pclose_f(parallelID, HDF5Error)
    print *, "HDF5 h5pclose_f status (Rank ", sml_mype, "): ", HDF5Error

    call h5close_f(HDF5Error)
    print *, "HDF5 h5close_f status (Rank ", sml_mype, "): ", HDF5Error

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    print *, "HDF5  complete"

    deallocate(Ff_t)
    deallocate(g1_t)
    deallocate(f1_t)

    deallocate(dfdt_col)
    deallocate(dfdt_Epar)
    deallocate(dfdt_Q_t)
    deallocate(dfdt_par)
    deallocate(dfdt_impar)

    deallocate(phi1t)
    deallocate(phi1t_k1)
    deallocate(phi1t_k2)
    deallocate(phi1t_k3)

    deallocate(den_t)
    deallocate(ene_t)
    deallocate(mom_t)
    deallocate(FSAn_t)
    deallocate(FSAp_t)
    deallocate(FSAT_t)

    deallocate(col_power_t)
    deallocate(col_dens_t)
    deallocate(col_inner_info_t)

    deallocate(time)

  end subroutine closeOutputFile
! -----------------------------------------------------------------------------------

  subroutine writeVariable_scalar(var, varname, i)

    real(kind=8), intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(0) :: dimensions
    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
       call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
       call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
       call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
       call h5dclose_f(dsetID, HDF5Error)
       call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeVariable_scalar

  subroutine writeVariable_integer(var, varname, i)

    integer, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(0) :: dimensions

    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
       call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (smy_mype.eq.0) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_integer


  subroutine writeVariable_real(var, varname, i)

    real(kind=8), dimension(:), intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID,HDF5Error)
    if (sml_mype.eq.0) then   
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then   
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID,HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_real


  subroutine writeVariable_1d(var, varname, i)

    real(kind=8), dimension(:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      endif
      dimensions = shape(var)
    end if
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      endif
      dimensions = shape(var)
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_1d

  subroutine writeVariable_2d(var, varname, i)

    real(kind=8), dimension(:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(2) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
    dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then   
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_2d

  subroutine writeVariable_3d(var, varname, i)

    real(kind=8), dimension(:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(3) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_3d

  subroutine writeVariable_4d(var, varname, i)

    real(kind=8), dimension(:,:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(4) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_4d

  subroutine writeVariable_5d(var, varname, i)

    real(kind=8), dimension(:,:,:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(5) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_5d

  subroutine writeVariable_6d(var, varname, i)

    real(kind=8), dimension(:,:,:,:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(6) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      endif
      dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif 
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_6d


  subroutine writeVariable_7d(var, varname, i)

    real(kind=8), dimension(:,:,:,:,:,:,:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(7) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      endif
      dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
    endif 
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeVariable_7d

  subroutine writeIntegerVariable_1d(var, varname, i)

    integer, dimension(:), allocatable, intent(in) :: var
    integer, intent(in) :: i
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: ierror

#ifdef HAVE_PARALLEL_HDF5
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var)
    endif
    call MPI_Bcast(dimensions,2*size(dimensions),MPI_INTEGER,0,MPI_COMM_WORLD,ierror) ! 2*size(dimensions) since there seems to be no MPI datatype for long integers in Fortran, while the HSIZE_T kind is a long integer
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      if (.not. allocated(var)) then
        print *,"Tried to write unallocated variable:",varname
        stop
      end if
      dimensions = shape(var) 
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_DOUBLE, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    endif
#endif

  end subroutine writeIntegerVariable_1d

  subroutine writeIntegerNoGroup(var,varname)

    Integer, intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: ierror

    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, H5T_NATIVE_INTEGER, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, H5T_NATIVE_INTEGER, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeIntegerNoGroup

  subroutine writeStringNoGroup(var,varname)
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dspaceID, dsetID, stringType
    integer(HSIZE_T), dimension(0) :: dimensions
    integer(HSIZE_T) :: stringLength
    integer :: ierror

    stringLength = len(var)
    dimensions = shape(var)

#ifdef HAVE_PARALLEL_HDF5
    call h5tcopy_f(H5T_FORTRAN_S1, stringType, HDF5Error)
    call h5tset_size_f(stringtype, stringLength, HDF5Error)
    call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
    call h5dcreate_f(HDF5FileID, varname, stringType, dspaceID, dsetID, HDF5Error)
    if (sml_mype.eq.0) then
      call h5dwrite_f(dsetID, stringType, var, dimensions, HDF5Error)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call h5dclose_f(dsetID, HDF5Error)
    call h5sclose_f(dspaceID, HDF5Error)
#else
    if (sml_mype.eq.0) then
      call h5tcopy_f(H5T_FORTRAN_S1, stringType, HDF5Error)
      call h5tset_size_f(stringtype, stringLength, HDF5Error)
      call h5screate_simple_f(rank(var), dimensions, dspaceID, HDF5Error)
      call h5dcreate_f(HDF5FileID, varname, stringType, dspaceID, dsetID, HDF5Error)
      call h5dwrite_f(dsetID, stringType, var, dimensions, HDF5Error)
      call h5dclose_f(dsetID, HDF5Error)
      call h5sclose_f(dspaceID, HDF5Error)
    end if
#endif

  end subroutine writeStringNoGroup

!  integer function rank_real(A)
!    real(kind=8), dimension(:), intent(in) :: A
!    rank_real=size(shape(A))
!    return
!  end function rank_real

  integer function rank_character(A)
    character(len=*), intent(in) :: A
    rank_character=size(shape(A))
    return
  end function rank_character

  integer function rank_integer(A)
    integer, intent(in) :: A
    rank_integer=size(shape(A))
    return
  end function rank_integer

  integer function rank_scalar(A)
    real(kind=8), intent(in) :: A
    rank_scalar=size(shape(A))
    return
  end function rank_scalar

  integer function rank_1d(A)
    real(kind=8), allocatable, intent(in) :: A(:)
    rank_1d=size(shape(A))
    return
  end function rank_1d

  integer function rank_1d_integer(A)
    integer, dimension(:), allocatable, intent(in) :: A
    rank_1d_integer=size(shape(A))
    return
  end function rank_1d_integer

  integer function rank_2d(A)
    real(kind=8), dimension(:,:), allocatable, intent(in) :: A
    rank_2d=size(shape(A))
    return
  end function rank_2d

  integer function rank_3d(A)
    real(kind=8), dimension(:,:,:), allocatable, intent(in) :: A
    rank_3d=size(shape(A))
    return
  end function rank_3d

  integer function rank_4d(A)
    real(kind=8), dimension(:,:,:,:), allocatable, intent(in) :: A
    rank_4d=size(shape(A))
    return
  end function rank_4d


  integer function rank_5d(A)
    real(kind=8), dimension(:,:,:,:,:), allocatable, intent(in) :: A
    rank_5d=size(shape(A))
    return
  end function rank_5d

  integer function rank_6d(A)
    real(kind=8), dimension(:,:,:,:,:,:), allocatable, intent(in) :: A
    rank_6d=size(shape(A))
    return
  end function rank_6d
end module writeHDF5Output

