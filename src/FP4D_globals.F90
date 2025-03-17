!> Module for global variables and important simulation parameters
module FP4D_globals
  !---------------------------------------------------------------
  ! Constant parameters
  !  
  real (kind=8), parameter :: sml_prot_mass = 1.6720D-27            !< proton mass (MKS)
  real (kind=8), parameter :: sml_e_charge  = 1.6022D-19            !< elementary charge (MKS)
  real (kind=8), parameter :: sml_ev2j      = sml_e_charge          !< eV --> J conversion
  real (kind=8), parameter :: sml_j2ev      = 1D0/sml_e_charge      !< J --> eV conversion
  real (kind=8), parameter :: sml_pi        = 3.1415926535897932    !< pi
  real (kind=8), parameter :: sml_2pi       = 6.2831853071795862D0  !< 2 pi
  real (kind=8), parameter :: sml_cvel      = 299792458
  integer, parameter :: sml_nsp_max         = 3                     !< Maximum number of species minus 1
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Input parameters - opt_param
  !  
    ! Physics Term On & Off or Option => 0: Off 1: On
  integer :: opt_Col = 1 !< Collision Operatoer in RHS
  integer :: opt_PS  = 1 !< f1 Parallel Streaming Term in LHS 
  integer :: opt_0M  = 1 !< Maxwellian Magnetic Drift Term in RHS
  integer :: opt_APS = 0 !< Alpha Particle Source Term in RHS
  integer :: opt_QLD = 0 !< Quasi-Linear Diffusion Term in RHS
  integer :: opt_NB  = 0 !< Neutral Beam Heating Term in RHS
  integer :: opt_Epar= 0
  integer :: opt_global = 0 !< Radially Global Simulation
  integer :: opt_restart = 0

    ! youknow
  integer :: opt_PS_aux  = 1 !< f1 Parallel Streaming Term in LHS 
  integer :: opt_QLD_aux = 1 !< f1 Parallel Streaming Term in LHS 
  integer :: opt_NB_aux = 1  !< f_module.F90 | NB source normalization 
  integer :: opt_profile_rho = 1
  integer :: opt_QLD_dist = 0 ! explicit Q(f) in RHS_module.F90 !(default) 0 : f0, 1: f0+f1

  integer :: opt_time_solver = 0 ! 0: Euler, 1: RK, -1: steady state
  integer :: opt_time_species = 99 ! only for 'opt_time_solver==-1'.. set 'b' in tau_ab
      ! set the default as 99 to get error if don't set this value

  integer :: opt_logaritm_NEO = 0 ! 0: defalut(NRL Plasma Formulary 2013, same with XGCa) 1:NEO (all species fixed )

  integer :: opt_remove_homogeneous = 0 ! 0: defalut 1:remove homogeneous sol using 


    ! youknow - rotation
  real (kind=8) :: omega_rot_in = 0
  integer :: opt_rotation = 0 ! 0: off, 1: anlytic phi_rot, 2: NEO

  real (8) :: smooth_tau = 0.0
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Input parameters - time_param
  !  
  real (kind=8) :: sml_dt = 0.25D0 !< Collision time step in units of the electron collision time
  real (kind=8) :: sml_tau = 0D0 !< collision time of 1st species
  integer :: sml_nstep = 10        !< Number of time steps

  integer :: start_step_Col = 1000000
  integer :: step_Col  = 1

  integer :: start_step_RHS = 1  ! start point for time stepping
  integer :: on_RHS     = 1      ! (1) : start (1) : start time-stepping (2) : 10000000,
  integer :: cycle_RHS  = 1      ! (1) : start (1) : start time-stepping (2) : 10000000,
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Input parameters - phy_param
  !  
  integer :: sml_isp = 1           !< Index of first species (0 is for electrons, 1 for main ion species, >1 for impurities)
  integer :: sml_nsp = 1           !< Index of last species
  
  real (kind=8), dimension(0:sml_nsp_max) :: sml_den    = (/ 1D19,9.95003D18,3.33D16,1.667D16 /)
  !^ Density of each species make sure that the sum of species 1-sml_nsp is equal to the density of species 0
  real (kind=8), dimension(0:sml_nsp_max) :: sml_t_ev   = (/ 1D3,1D3,1D3,1D3 /)                   !< Temperature of each species
  real (kind=8), dimension(0:sml_nsp_max) :: sml_flow   = (/ 0.0,0.0,0.0,0.0 /)                   !< Flow of each species
  real (kind=8), dimension(0:sml_nsp_max) :: sml_mass   = (/ 5.448D-4,2D0,7D0,12D0 /)             !< Mass of each species in units of the proton mass
    ! sml_charge     :: Ze
  real (kind=8), dimension(0:sml_nsp_max) :: sml_charge = (/ -1D0,1D0,3D0,6D0 /)                  !< Charge of each species, in units of the elementary charge
    ! sml_charge_num :: Z
  real (kind=8), dimension(0:sml_nsp_max) :: sml_charge_num = (/ -1D0,1D0,3D0,6D0 /)                  !< Charge number
  real (kind=8), dimension(0:sml_nsp_max) :: dlnn_dr= (/-1D0, -1D0, -1D0, -1D0/)
  real (kind=8), dimension(0:sml_nsp_max) :: dlnT_dr= (/-1D0, -1D0, -1D0, -1D0/)
  real (kind=8), dimension(0:sml_nsp_max) :: sml_Epar = (/-1D0, -1D0, -1D0, -1D0/)
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Input parameters - APS_param & NB_param
  !  
    ! Auxiliary heating variables
  integer :: APS_sp = 2             !< Alpha Particle Source Species
  integer :: NB_sp = 1              !< Ionized NB source Species

  real (kind=8), dimension(0:sml_nsp_max) :: sml_source_den=(/1D19,9.95003D18,3.33D16,1.667D16/)  !< Source density of each species 
  real (kind=8), dimension(0:sml_nsp_max) :: sml_source_t_ev= (/1D3,1D3,1D3,1D3 /)                !< Source temperature of each species
  real (kind=8), dimension(0:sml_nsp_max) :: RF_W0= (/0D0, 0D0, 0D0, 0D0/)          
  real (kind=8), dimension(0:sml_nsp_max) :: RF_omega= (/0D0, 0D0, 0D0, 0D0/)
  real (kind=8), dimension(0:sml_nsp_max) :: RF_width= (/0.1D0, 0.1D0, 0.1D0, 0.1D0/)
  integer, dimension(0:sml_nsp_max) :: RF_n_harmonic=(/0,0,0,0/)
  real (kind=8), dimension(0:sml_nsp_max) :: RF_kpar=(/1D0, 1D0, 1D0, 1D0/)
  real (kind=8) :: sml_source_width = 1D0
  real (kind=8) :: pwrscale = 1D0 ! youknow - RF implicit (TORIC)

  real (kind=8) :: NB_pwrscale = 1D0  ! youknow
  integer :: NB_sink_opt = 0          ! youknow 
  integer :: NB_sink_sp  = 0          ! youknow

  real (kind=8) :: NB_fast_den_factor = 1D-5  ! youknow 

  integer :: NB_num = 20000

  integer :: start_Q = 1              ! youknow
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Input parameters - geo_param
  !  
  real (kind=8) :: deltapsi = 1D0 ! youknow
  integer :: sml_eqb_opt = 0 !Equilibrium Option 

    ! FROM 'eqmodule.F90' -- eqdsk
  !!!! inputs !!!!!!
  public :: xx0,xx1,q0,q1,q2,q3,q4,a0,R0,B0,ase_q0,ase_a0,ase_R0,ase_B0
  public :: ase_kappa,ase_del,ase_shift
  public :: npsi,nth,ifd0,ifd1,op_geo_info
  public :: ipe ! youknow

  ! integer
  integer :: npsi = 1
  integer :: nth = 32
  integer :: ifd0 = 3
  integer :: ifd1 = 27
  integer :: op_geo_info = 0
  integer :: ipe = 1      ! youknow - | 0 : print on | 1: print off 

  ! real
  real :: xx0 = 0.4
  real :: xx1 = 0.5
  real :: q0 = 1.4
  real :: q1 = 0
  real :: q2 = 0
  real :: q3 = 0
  real :: q4 = 0
  real :: a0 = 0.5
  real :: R0 = 1.8
  real :: B0 = 1.8
  real :: ase_q0 = 1.4
  real :: ase_a0 = 0.5
  real :: ase_R0 = 0.48
  real :: ase_B0 = 1.3
  real :: ase_kappa = 1.9
  real :: ase_del = 1.5
  real :: ase_shift = 0
  !---------------------------------------------------------------



  !---------------------------------------------------------------
  ! Input parameter - ??
  !  
  real (kind=8) :: neg_frac = 0.3D0
  ! integer :: sup_val = 10
  integer :: sup_val2=1
  ! integer :: sup_step = 10000
  integer :: sup_step2=1

  integer :: mx_sp = 2              !< RF momentum, energy sink species
  real (kind=8) :: fw = 1.0
  real (kind=8) :: fr = 3.0

  integer :: itmax=1000
  !---------------------------------------------------------------



  !---------------------------------------------------------------
  ! Input parameter - Numerical Resolution
  !  
  ! integer :: f_isp = 0                           !< Index of first species (copied from sml_isp)
  ! integer :: f_nsp = 1                           !< Index of last species (copied from sml_nsp; total number of species: f_nsp-f_isp+1)
  integer :: f_nmu = 44                          !< Grid size in v_perp direction --> N_perp = f_nmu+1
  integer :: f_nvp = 44                          !< Grid size in v_para direction --> N_para = 2*f_nvp+1
  ! integer :: f_npsi = 1
  ! integer :: f_nth = 10

  integer :: size_f0
  real (kind=8) :: f_smu_max = 4D0               !< Perp. velocity cutoff in units of the thermal velocity
  real (kind=8) :: f_vp_max = 4D0                !< Parallel velocity cutoff in units of the thermal velocity
  ! real (kind=8) :: f_dsmu                        !< Grid resolution in v_perp direction
  ! real (kind=8) :: f_dvp                         !< Grid resolution in v_para direction
  ! real (kind=8) :: f_dth
  real (kind=8) :: f_velo_aniso = 0.3D0          !< Velocity anisotropy of the initial distribution function
                                                  !! (T_para = f_velo_aniso*f_t_ev) 

  ! Resolution
  integer :: nspecies
  ! integer :: npsi, nth  ! Defined part of eqmodule.F90
  integer :: nthm1
  integer :: nvr, nvrm1
  integer :: nvz, nvzm1

  ! range of index
  integer :: is_s , is_e
  integer :: tar_s, tar_e
  integer :: ip_s , ip_e
  integer :: ith_s, ith_e
  integer :: ir_s , ir_e
  integer :: iz_s , iz_e
  !---------------------------------------------------------------

  ! Kareny - 1981
  ! Set the range of D_rr,RF in 'normalized velocity'
    ! except the range, the value of D_rr,RF = 0
  real (kind=8) :: RF_vr1 = 0, RF_vr2 = 1
  real (kind=8) :: RF_vz1 = 0, RF_vz2 = 1


  !---------------------------------------------------------------
  ! 'quasi_f_module.F90'
  !
    ! adiabatic electron
  real(kind=8), allocatable :: n0_ae (:,:)
  real(kind=8), allocatable :: T0_ae (:)

  integer :: opt_Potential = 10
    ! If opt_Potential is 0, 
      ! g1=f1
    ! If opt_Potential is not 0,
      ! it is automatically determine @ quasi_f_setup
  integer :: opt_Potential_converge
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! FROM 'FP4D_equilibrium.F90'
  !
  ! 1D-variable
  real (kind=8), allocatable :: eq_psin(:),eq_intJ(:),eq_theta(:),eq_dpsidrho(:)
  real (kind=8), allocatable :: eq_Volume(:,:) ! real space volume element ! for rho not psi. rho :: sqrt((psi-psi0)/(psib-psi0))
   
  ! 2D-variable
  ! FROM eq_data < arr_B.txt
  ! dimension (0:Nth, 1:Npsi)
  real (kind=8), allocatable, dimension(:,:) :: eq_B, eq_Bt,  eq_Bp
  real (kind=8), allocatable, dimension(:,:) :: eq_dBdth, eq_dBdpsi
  real (kind=8), allocatable, dimension(:,:) :: eq_I, eq_J, eq_invJ

  real (kind=8), allocatable :: geo_R(:,:),  geo_Z(:,:), geo_psi(:,:)
  real (kind=8), allocatable :: geo_dRdth(:,:)

  ! rotation
          ! Belli_2009_PPCF eq(13)
          ! only for (e,i) case or (e,i,W) w/ trace limit
          ! isp=0 -> e , isp=1 -> D
  real (kind=8), allocatable :: geo_FSA_R2  (:)

  ! phi_rot      :: phi_0_tilda refer : eq(10) in [Belli_2009]
  ! phi_rot_star :: phi_0_tilda refer : eq(10) in [Belli_2009]
  real (kind=8), allocatable :: phi_rot          (:,:) ! Belli_2009_PPCF eq(13) ... temporary --> only for (e,i) case
  real (kind=8), allocatable :: phi_rot_star     (:,:)
  real (kind=8), allocatable :: phi_rot_avg      (:)   
  real (kind=8), allocatable :: phi_rot_star_avg (:)   

  real (kind=8), allocatable :: omega_rot   (:)    
  real (kind=8), allocatable :: lambda_rot  (:,:,:)
  real (kind=8), allocatable :: dlambda_rot_dth  (:,:,:)
  real (kind=8), allocatable :: den_rot     (:,:) ! eq_den(ip,is) -> den_rot(ip,is)
  ! den_rot = N0(psi)  
  real (kind=8), allocatable :: small_n0     (:,:,:) 

  real (kind=8), allocatable :: eq_omega  (:,:,:)
  real (kind=8), allocatable :: eq_phi1   (:,:),  eq_dphidth(:,:)
  real (kind=8), allocatable :: eq_phi1_k1(:,:),  eq_dphidth_k1(:,:) ! youknow
  real (kind=8), allocatable :: eq_phi1_k2(:,:),  eq_dphidth_k2(:,:)
  real (kind=8), allocatable :: eq_phi1_k3(:,:),  eq_dphidth_k3(:,:)
  real (kind=8), allocatable :: eq_phi1_tmp(:,:)

  real (kind=8), allocatable :: prf_t_ev  (:,:),  prf_dlnT_dr(:,:)
  real (kind=8), allocatable :: prf_den   (:,:),  prf_dlnn_dr(:,:)
  real (kind=8), allocatable :: eq_profile(:,:)
  real (kind=8), allocatable :: eq_t_ev   (:,:),  eq_dlnT(:,:)
  real (kind=8), allocatable :: eq_den    (:,:),  eq_dlnn(:,:)

  real (kind=8), allocatable :: eq_lambda(:,:,:), eq_tau(:,:,:),  eq_nu(:,:,:)

  real (kind=8) :: eq_drho
  real (kind=8) :: eq_dtheta
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! From 'f_module.F90'
  !
  real (kind=8), allocatable :: vth(:,:)

  ! vel_grid
  real (kind=8) :: norm_dr, norm_dz
  real (kind=8), allocatable :: norm_vr (:), norm_vz (:) ! norm_vr (0)=norm_dr/3 ! vel_volume & RHS module (PS, QL)          
  real (kind=8), allocatable :: norm_vr0(:)              ! norm_vr0(0)=0         ! only for defining maxwellian + Kareny_1981
  real (kind=8), allocatable :: vel_volume(:,:,:) ! volume element. not total volume
  ! real (kind=8), allocatable :: conv_factor(:,:) --> before add, check the overlap

  ! drift component
  real (kind=8), allocatable, dimension(:,:,:) :: drift_psi, drift_theta
  real (kind=8), allocatable, dimension(:,:,:) :: drift_rot_1, drift_rot_2

  ! distribution function
  real (kind=8), allocatable :: Ff      (:,:,:,:,:)     !< Distribution function
  real (kind=8), allocatable :: tmp_Ff  (:,:,:,:,:)     !< Distribution function !youknow makes this var to incarnate the option of ineg in cql3d
  real (kind=8), allocatable :: dfdt_c  (:,:,:,:,:)  !< Change dt*C(f) of distribution function
  real (kind=8), allocatable :: dfdt_q  (:,:,:,:,:)  !< Change dt*Q(f) of distribution function
  real (kind=8), allocatable :: dfdt_p  (:,:,:,:,:)  !< Change dt*PS(f) of distribution function
  real (kind=8), allocatable :: f1      (:,:,:,:,:)  !< Change dt*PS(f) of distribution function
  real (kind=8), allocatable :: f1_an   (:,:,:,:,:)  
  real (kind=8), allocatable :: f1_neo  (:,:,:,:,:)  
  real (kind=8), allocatable :: dfdt_d  (:,:,:,:,:)  !< Change dt*MD(f) of distribution function
  real (kind=8), allocatable :: df_sts  (:,:,:,:,:) !! Source distribution !Hyeonjun Lee
  real (kind=8), allocatable :: sts_int (:,:,:)
  real (kind=8), allocatable :: df_stl  (:,:,:,:,:) !! Sink distribution !Hyeonjun Lee
  real (kind=8), allocatable :: qld_B   (:,:,:,:), qld_C(:,:,:,:),qld_F(:,:,:,:)!youknow_TORIC
  real (kind=8), allocatable :: dfdt_e  (:,:,:,:,:)
  real (kind=8), allocatable :: dfdt_imp(:,:,:,:,:)
  real (kind=8), allocatable :: k1(:,:,:,:,:), k2(:,:,:,:,:), k3(:,:,:,:,:), k4(:,:,:,:,:)
  real (kind=8), allocatable :: frv(:,:,:,:)
  real (kind=8), allocatable :: g1      (:,:,:,:,:)
  real (kind=8), allocatable :: g1_k1   (:,:,:,:,:) ! youknow
  real (kind=8), allocatable :: g1_k2   (:,:,:,:,:)
  real (kind=8), allocatable :: g1_k3   (:,:,:,:,:)
  real (kind=8), allocatable :: g1_k4   (:,:,:,:,:)
  real (kind=8), allocatable :: f0_M    (:,:,:,:,:) !<
  real (kind=8), allocatable :: f1_tmp  (:,:,:,:,:)

  real (kind=8), allocatable :: f0_M_tar(:,:,:,:,:) ! youknow
  real (kind=8), allocatable :: g1_tar  (:,:,:,:,:) ! youknow
  real (kind=8), allocatable :: f1_tar  (:,:,:,:,:) ! youknow

  ! collisional power deposition
  real (kind=8), allocatable :: col_power(:,:,:,:) ! youknow ! Npsi, Nth, colliding sp., target sp.
  real (kind=8), allocatable :: col_dens(:,:,:,:) ! youknow ! Npsi, Nth, colliding sp., target sp.

  real (kind=8), allocatable :: col_inner_info(:,:,:,:,:,:) ! youknow ! Npsi, Nth, colliding sp., target sp.

  ! use in collision module
  real (kind=8), allocatable :: num_n(:,:,:),num_T(:,:,:),num_p(:,:,:)
  real (kind=8), allocatable :: FSAn(:,:), FSAT(:,:),FSAp(:,:)
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! From 'FP4D_Post.F90'
  !
  real (kind=8), allocatable, dimension(:,:,:) :: pflux, eflux, mflux, pflux_col, eflux_col, mflux_col
  real (kind=8), allocatable, dimension(:,:) :: FSA_pflux, FSA_eflux, FSA_mflux, FSA_pflux_col, FSA_eflux_col, FSA_mflux_col
  real (kind=8), allocatable, dimension(:,:) :: jpar, jparB
  real (kind=8), allocatable, dimension(:) :: FSA_jpar, FSA_jparB
  real (kind=8), allocatable, dimension(:,:,:) :: post_dens, post_NUpar, post_NUparB
  real (kind=8), allocatable, dimension(:,:) :: FSA_post_dens, FSA_post_NUpar, FSA_post_NUparB
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! MPI variables
  !
  integer :: sml_totalpe   = 1  !< Number of MPI ranks
  integer :: sml_mype      = 0  !< MPI rank ID
  integer :: sml_comm      = 0  !< Global communicator (MPI_COMM_WORLD)
  integer :: sml_comm_null = 0  !< Self-communicator

  integer :: root      = 0      ! JPL: root processor inder
  integer :: np1, np2           ! JPL: number of mpi particion in psi and theta
  integer, dimension(:), allocatable :: istart1, iend1, istart2, iend2, ip1_mpi
  !JPL: psi (1) and theta (2) start and end index for each pe
  integer ::nptot_mpi
  !OpenMP
  integer :: sml_nthreads  = 1  !< Number of OpenMP threads (not an inputvariable)
  integer :: col_f_nthreads = 1        !< number of OpenMP threads used in the inner nested level
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! output file
  !
  character(len=80)  :: runfile_FP4Dout = 'out.fp4d.run'
  integer :: io_FP4Dout = 12

  character(len=80)  :: runfile_NEO_input = 'out.fp4d.neo_input'
  integer :: io_NEO = 13

  character (len=16) :: hdf5_name='data00.h5'
  !---------------------------------------------------------------

  ! input file
  character (len=*), parameter:: filename_profile='profile.dat'
  integer :: io_profile = 1002


  ! Time and time step
  real (kind=8) :: sml_time=0d0     !< Time of the simulation in seconds
  integer :: sml_istep          !< Time step index
  
  integer :: out_h5 = 1, out_step = 1, ndstep=1


  ! Collision

  ! youknow
  integer :: col_iteration = 10     ! youknow 

  real (kind=8) :: den_converge = 1D-8
  real (kind=8) :: ene_converge = 1D-8
  real (kind=8) :: mom_converge = 1D-8 ! do not used


  ! youknow - special option
  integer :: special_tau_W = 0    ! only  only for 'opt_time_solver==-1'. If this option is 1. then, we set the eq_tau(0,1,2) = eq_tau(0,1,1)
  integer :: special_tau_He3 = 0  ! only  only for 'opt_time_solver==-1'. If this option is 1. then, we set the eq_tau(0,1,3) = eq_tau(0,1,1)



  contains

  subroutine setup_FP4D_globalvars()

    ! Set the index

    nspecies = sml_nsp - sml_isp + 1
    is_s  = sml_isp ;   is_e  = sml_nsp
    tar_s = 0       ;   tar_e = sml_nsp

    ! npsi  = f_npsi
    ip_s  = 1       ;   ip_e  = npsi

    ! nth   = f_nth
    nthm1 = nth-1
    ith_s = 0       ;   ith_e   = nthm1   ! 0:nthm1

    nvr   = f_nmu+1
    nvrm1 = f_nmu
    ir_s  = 0       ;   ir_e   = nvrm1    ! 0:nvrm1

    nvz   = 2*f_nvp+1
    nvzm1 = nvz-1
    iz_s  = -nvzm1/2;   iz_e   = nvzm1/2  !-f_nvp : f_nvp


    write(*,'(A12, I5, I5)') 'is idx' , is_s,  is_e
    write(*,'(A12, I5, I5)') 'tar idx', tar_s, tar_e
    write(*,'(A12, I5, I5)') 'ip idx' , ip_s,  ip_e
    write(*,'(A12, I5, I5)') 'ith idx', ith_s, ith_e
    write(*,'(A12, I5, I5)') 'ir idx' , ir_s,  ir_e
    write(*,'(A12, I5, I5)') 'iz idx' , iz_s,  iz_e


    ! Set the dtheta??


    ! f_vp_max, f_smu_max
    ! psi min psi max .... and so on...
    ! dt
    ! tau
    ! mesh_r, mesh_z
    ! vth


  end subroutine setup_FP4D_globalvars

end module FP4D_globals
