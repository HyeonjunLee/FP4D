#ifdef USE_ASYNC
#define ASYNC(istream)  async(istream)
#else
#define ASYNC(istream)
#endif

module col_f_module

use FP4D_globals, only : istart1, istart2
use FP4D_globals, only : iend1, iend2
use FP4D_globals, only : sml_mype
use FP4D_globals, only : col_f_nthreads
use FP4D_globals, only : opt_Col, opt_QLD
use FP4D_globals, only: den_converge, ene_converge, col_iteration

integer :: col_f_nvr, col_f_nvz      !< The number of "GRIDS" in vr and vz, respectively
integer :: col_f_nvrm1, col_f_nvzm1  !< To enhance understandings of compiler and conseq. to avoid temp. copy
integer :: col_f_nvrm1_nvzm1         !< To enhance understandings of compiler and conseq. to avoid temp. copy
integer :: col_f_nth, col_f_nthm1
integer :: col_f_npsi
real (kind=8) :: col_f_vp_max        !< Max v_para of velocity-grid
real (kind=8) :: col_f_smu_max       !< Max v_perp of velocity-grid
real (kind=8) :: col_f_dvp           !< v_para grid spacing
real (kind=8) :: col_f_dsmu          !< v_perp grid spacing
real (kind=8) :: col_f_dth
integer :: col_f_ntotal_v            !< set in col_f_setup
real (kind=8) :: col_f_dt            !< Time step for collision operation
real (kind=8) :: QLD_f_dt            !< Time step for QLD operation ! youknow
real (kind=8) :: rf_pwrscale ! youknow

integer, allocatable, dimension(:,:) :: index_map_LU
integer, allocatable, dimension(:) :: LU_cvalues, LU_rowindx, LU_colptr
integer :: LU_n, LU_nnz, LU_nrhs, LU_ldb

type col_f_sp_type
  integer :: spid    !! [setup_col_sp] species ID
  integer :: nodeid  !! [setup_col_sp] (mesh) node ID
  real (kind=8) :: T_avg  !! [setup_col_sp]
  real (kind=8) :: mass_au !! [setup_col_sp]
  real (kind=8) :: charge_eu , Ze!! [setup_col_sp]
  real (kind=8) :: vpar_beg   !! [setup_col_sp]
  real (kind=8) :: den_lambda_gamma !! [setup_col_sp]
  real (kind=8) :: conv_factor !! [get_local_f<setup_col_sp]
  real (kind=8), allocatable, dimension(:,:) :: pdf_n !![get_local_f<setup_col_sp] Given initial PDF
  real (kind=8), allocatable, dimension(:,:) :: pdf_np1!![get_local_f<setup_col_sp] To be used as PDF at next Picard step. BUT has a role of "df" at end
  real (kind=8), allocatable, dimension(:,:,:) :: ED   ! 5 x (col_f_nvr-1) x (col_f_nvz-1)

  ! JPL :: find collisional power partition
  real (kind=8), allocatable, dimension(:) :: power
  real (kind=8), allocatable, dimension(:) :: col_dens
  real (kind=8), allocatable, dimension(:,:,:) :: inner_info ! youknow

     !!!! ## ----- For the purpose of diagnosis (print_snapshot)
     !!real (kind=8) :: numeric_T_par, numeric_T_perp
     !!real (kind=8) :: numeric_dens

     !! Below is from obsolete col_f_core_type
  real (kind=8) :: numeric_vth2, numeric_T, &  !![col_fm_core_init<col_fm_core]
                   numeric_vtheq2, vth !! numeric_vtheq2 is for compatability with col_f_core_type. numeric_Teq is not needed.
  real (kind=8) :: mass, mesh_dz, mesh_dr, &  !![setup_col_sp]
                   dens, ens, mom ! [col_fm_core_init<col_fm_core]
                    
  real (kind=8), allocatable, dimension(:) :: mesh_r, mesh_r_half, &  !![setup_col_sp]
                                                 mesh_z, mesh_z_half, &  !![setup_col_sp]
                                                 local_center_volume, &  !![setup_col_sp]
                                                 vol !![setup_col_sp]
     real (kind=8), allocatable, dimension(:,:) :: delta_r, delta_z
end type col_f_sp_type

  type(col_f_sp_type), allocatable :: col_sp_cell_all(:,:)

  integer :: col_f_sp_s, col_f_sp_e               !! [col_f_setup] the (starting, ending) index of species to be collided
                                                  !!              i.e. species to be collided = sp[col_f_sp_s, col_f_sp_e]
                                                  !!              You can play with this parameter to choose part of species to be collided
  
  integer :: tar_f_sp_s, tar_f_sp_e               !! [col_f_setup] the (starting, ending) index of species to be collided
  ! colliding species : [col_f_sp_s, col_f_sp_e]
  ! traget    species : [tar_f_sp_s, tar_f_sp_e]



  integer :: col_f_sp_num                         !! [col_f_setup] the number of species to be collided ( = col_f_sp_index_e-col_f_sp_index_s+1 )
  integer :: col_f_i_nthreads                     !! [col_f_setup] the number of inner OMP threads
  integer :: col_f_o_nthreads                     !! [col_f_setup] the number of outer OMP threads

  !! <BIG MEMORY ALLOCATION>
  !! Rough estimate on (minimal) memory requirements for N spcies / 1 outer thread
  !! N x 1.10 Mbyte + N(N-1) x 19.8 Mbyte
  !! [1  Outer thread] N=1 : 1.10 Mbyte, N=2 : 61.6 Mbyte, N=3 : 122  Mbyte, N=4 : 242  Mbyte
  !! [8  Outer thread] N=1 : 8.80 Mbyte, N=2 : 493  Mbyte, N=3 : 977  Mbyte, N=4 : 1.89 Gbyte
  !! [16 Outer thread] N=1 : 17.6 Mbyte, N=2 : 986  Mbyte, N=3 : 1.90 Gbyte, N=4 : 3.78 Gbyte

  real(kind=8), allocatable, dimension(:,:,:) :: gammac_spall  !![setup_col_sp]
  !! ## ----- Memory for 1 self-interaction, M_s  : 8 byte x (31x31x63) x5 = 2.31Mbyte (Kernel default)
  !! ##                                           :        x (31x31x30) x5 = 1.10Mbyte (XGC default)
  !!   -  Self-interaction can be cut down to x3/5 memory using the conservation identity
  real (kind=8), dimension(:,:,:,:,:), allocatable :: M_s
  !! ## ----- Memory for 1 cross-interaction, M_ab : 8 byte x (31x63x31x63) x3 = 87.3Mbyte  (Kernel default)
  !! ##                                            :        x (31x30x31x30) x3 = 19.8Mbyte  (XGC default)
  real (kind=8), dimension(:,:,:,:,:,:), allocatable :: M_ab
  !!!!DIR$ ATTRIBUTES FASTMEM :: M_ab    !!!quad-cache mode will cause memory allocation error. quad-flat with numactl has not been tested yet

  !! <SOLVER-related>
  integer, parameter :: vpic_inner_iter_max = 20

  integer, private, parameter :: MAX_THREADS=48

  ! ------------------------------------------------------------
  ! set use_superlu = .true. to use sparse direct SuperLU solver
  ! set use_superlu = .false. to use lapack band solver
  ! ------------------------------------------------------------
  logical, parameter :: use_superlu = .false.

  integer ::  col_f_mat_list(MAX_THREADS)
  integer ::  col_f_vecb_list(MAX_THREADS)
  integer ::  col_f_vecx_list(MAX_THREADS)
  integer ::  col_f_ksp_list(MAX_THREADS)


  contains
    ! BLAS solver collision operator time-stepping
#include "bsolver.F90"

    !> Set up global variables of the collision operator
    !! @param[in]   isp   Index of 1st particle species, integer
    !! @param[in]   nsp   Index of last particle species, integer
    !! @param[in]   nmu   Grid size in v_perp direction, integer
    !! @param[in]   nvp   Grid size in v_parallel direction, integer
    !! @param[in]   dt    Time step, real(8)
    !! @param[in]   vp_max   v_para velocity cut-off
    !! @param[in]   smu_max   v_perp velocity cut-off
    !! @param[in]   dvp   v_para grid spacing
    !! @param[in]   dsmu   v_perp grid spacing
    !! @param[in]   mype   MPI rank ID
    !!
    subroutine col_f_setup(isp,nsp,npsi,nth,nmu,nvp,dt,vp_max,smu_max,pwrscale,mype)
      use FP4D_globals, only:sml_2pi

      implicit none
      integer, intent(in) :: isp, nsp, nmu, nvp, mype, nth, npsi
      real (kind=8), intent(in) :: dt, vp_max, smu_max
      real (kind=8), intent(in) :: pwrscale
      !
      integer :: i, j
      integer :: mat_pos_rel(9), mat_pos_rel_indx(9), LU_i_arr(9)
      integer :: elem_n, mat_pos, LU_i, LU_j, incr_LU
      integer, allocatable, dimension(:) :: LU_colptr_num
    
      col_f_sp_s = isp
      col_f_sp_e = nsp
      col_f_sp_num = col_f_sp_e - col_f_sp_s + 1

      ! youknow
      ! Adiabatic electron only target of collision

      ! target_f_sp_s always start from 0
      tar_f_sp_s = 0
      tar_f_sp_e = nsp

      col_f_i_nthreads = col_f_nthreads
      !call set_omp_outer_thread_num( col_f_o_nthreads )
      col_f_o_nthreads = 1  !JPL: use mpi instead of  out_openmp

      if(mype .eq. 0) then
          write(*,'(a,i4)') '[MULTI-COL] inner thread = ', col_f_i_nthreads
          write(*,'(a,i4)') '[MULTI-COL] outer thread = ', col_f_o_nthreads
      endif

      col_f_npsi        = npsi
      col_f_nth         = nth
      ! col_f_dth         = dth
      col_f_dth         = sml_2pi/real(nth,8)   ! same with eq_dtheta
      col_f_nthm1       = nth-1
      col_f_nvr         = nmu+1
      col_f_nvrm1       = nmu
      col_f_nvz         = nvp*2 + 1
      col_f_nvzm1       = col_f_nvz-1
      col_f_ntotal_v    = col_f_nvr * col_f_nvz
      col_f_nvrm1_nvzm1 = col_f_nvrm1 * col_f_nvzm1
      col_f_dt          = dt ! delta time for collision operation
      QLD_f_dt          = dt ! youknow
      col_f_vp_max      = vp_max
      col_f_smu_max     = smu_max
      ! col_f_dvp         = dvp
      ! col_f_dsmu        = dsmu
      col_f_dvp         = vp_max/real(nvp,8)    ! same with norm_dz
      col_f_dsmu        = smu_max/real(nmu,8)   ! same with norm_dr

      rf_pwrscale       = pwrscale

      ! 2013-02-23 =============================================== SUPER LU
      ! index_map_LU, LU_rowindx, LU_cvalues
      LU_n = col_f_ntotal_v !global
      LU_nnz = 4*4 + 6*2*(col_f_nvr-2)+6*2*(col_f_nvz-2)+9*(col_f_nvr-2)*(col_f_nvz-2) !global
      LU_nrhs = 1 !global
      LU_ldb = LU_n !global
      !$omp critical (alloc1)
      allocate(LU_rowindx(LU_nnz), LU_cvalues(LU_n), LU_colptr(LU_n+1), index_map_LU(9,LU_n))       !global
      allocate(LU_colptr_num(LU_n))    !local
      !$omp end critical (alloc1)

      LU_colptr_num = 0   ! number of elements in each column, local
      LU_rowindx = 0         ! global

      !below is time independent. Move to init_col in module
      LU_colptr(1) = 1
      do i=1,LU_n
          !for colptr
          LU_i = (i-1) / col_f_nvr+1
          LU_j = mod(i-1, col_f_nvr)+1
          if( (LU_i .eq. 1) .or. (LU_i .eq. col_f_nvz) ) then
              if( (LU_j .eq. 1) .or. (LU_j .eq. col_f_nvr) ) then
                  incr_LU = 4
              else
                  incr_LU = 6
              endif
          else
              if( (LU_j .eq. 1) .or. (LU_j .eq. col_f_nvr) ) then
                  incr_LU = 6
              else
                  incr_LU = 9
              endif
          endif
          LU_colptr(i+1) = LU_colptr(i)+incr_LU
      enddo

      !===============
      !  3--6--9
      !  |  |  |
      !  2--5--8
      !  |  |  |
      !  1--4--7
      !==============
      mat_pos_rel=(/-col_f_nvr-1,-col_f_nvr, -col_f_nvr+1, -1, 0, 1, col_f_nvr-1, col_f_nvr, col_f_nvr+1/)
      index_map_LU = 0
      do i=1,col_f_nvz
          do j=1, col_f_nvr
              mat_pos = j+(i-1)*col_f_nvr
      
              ! INDEXING
              if(i .eq. 1) then
                  if(j .eq. 1) then
                      !(J,I)=(1,1)
                      elem_n = 4
                      mat_pos_rel_indx(1:elem_n) =(/5,6,8,9/)
                   elseif (j .eq. col_f_nvr) then
                      !(J,I)=(Nr,1)
                      elem_n = 4
                      mat_pos_rel_indx(1:elem_n)=(/4,5,7,8/)
                   else
                      !(J,I)=(:,1)
                      elem_n = 6
                      mat_pos_rel_indx(1:elem_n)=(/4,5,6,7,8,9/)
                   endif
               elseif(i .eq. col_f_nvz) then
                   if(j .eq. 1) then
                      !(J,I)=(1,mesh_Nz)
                      elem_n = 4
                      mat_pos_rel_indx(1:elem_n)=(/2,3,5,6/)
                   elseif (j .eq. col_f_nvr) then
                      !(J,I)=(Nr,mesh_Nz)
                      elem_n = 4
                      mat_pos_rel_indx(1:elem_n)=(/1,2,4,5/)
                   else
                      !(J,I)=(:,mesh_Nz)
                      elem_n = 6
                      mat_pos_rel_indx(1:elem_n)=(/1,2,3,4,5,6/)
                   endif
               else
                   if(j .eq. 1) then
                      !(J,I) = (1,:)
                      elem_n = 6
                      mat_pos_rel_indx(1:elem_n)=(/2,3,5,6,8,9/)
                   elseif(j .eq. col_f_nvr) then
                      !(J,I) = (mesh_Nr,:)
                      elem_n = 6
                      mat_pos_rel_indx(1:elem_n)=(/1,2,4,5,7,8/)
                   else
                      !(J,I) = (:,:)
                      elem_n = 9
                      mat_pos_rel_indx(1:elem_n)=(/1,2,3,4,5,6,7,8,9/)
                   endif
               endif
               LU_i_arr(1:elem_n) = mat_pos+mat_pos_rel(mat_pos_rel_indx(1:elem_n))  ! I need to change LU_i to array
               index_map_LU(mat_pos_rel_indx(1:elem_n),mat_pos) = LU_colptr(LU_i_arr(1:elem_n))+LU_colptr_num(LU_i_arr(1:elem_n))
               LU_colptr_num(LU_i_arr(1:elem_n)) = LU_colptr_num(LU_i_arr(1:elem_n))+1
               LU_rowindx(index_map_LU(mat_pos_rel_indx(1:elem_n),mat_pos)) = mat_pos
               LU_cvalues(mat_pos) = index_map_LU(5, mat_pos)  !For implicit time marching
          enddo
      enddo
      !$omp critical (alloc1)
      deallocate(LU_colptr_num)
      !$omp end critical (alloc1)

      return

    end subroutine col_f_setup

    !> Deallocate global arrays of col_f_module
    subroutine col_f_finalize
      implicit none

      !$omp critical (alloc1)
      deallocate(LU_rowindx, LU_cvalues, LU_colptr, index_map_LU)
      !$omp end critical (alloc1)

      return

    end subroutine col_f_finalize

    !> Dummy routine for domain decomposition
    !> In full XGC, this would have some more logic
    subroutine get_mesh_node_range_for_this_mpirank(inode1_in, inode2_in, inode1, inode2, stride )
      implicit none
      integer, intent(in) :: inode1_in, inode2_in
      !
      integer :: inode1, inode2, stride

      inode1=inode1_in
      inode2=inode2_in
      stride = 1

      return

    end subroutine get_mesh_node_range_for_this_mpirank


    !> This routine decides how the workload is distributed
    !! between threads --> How many grid points for each thread
    subroutine get_mesh_node_range_for_threads_in_this_mpirank( inode1, inode2, stride,  &
                                              nthreads, isize, i_beg, i_end, i_stride )
      use omp_module, only : split_indices
      implicit none
      integer, intent(in) :: inode1, inode2, stride, nthreads
      integer :: isize, i_stride
      integer :: i_beg(size(col_f_mat_list)), i_end(size(col_f_mat_list))
      integer :: ith

#if defined(_OPENMP)
      isize = min(inode2-inode1+1,size(col_f_mat_list),nthreads)
      if (isize .gt. 0) then
         call split_indices(inode2-inode1+1,isize,i_beg,i_end)
         i_beg(1:isize) = i_beg(1:isize) + (inode1-1)
         i_end(1:isize) = i_end(1:isize) + (inode1-1)
      else
         i_beg(:)=2
         i_end(:)=1
      endif
      i_stride = stride
#else
      isize = 1
      i_beg(1) = inode1
      i_end(1) = inode2
      i_stride = stride
#endif

      return

    end subroutine get_mesh_node_range_for_threads_in_this_mpirank


    !> compute the number of threads used in the outer
    !! level OpenMP parallelization.
    !! The number of inner threads is set with an input
    !! parameter (col_f_nthreads)
    !! --> number of outer threads = OMP_NUM_THREADS/#(inner threads)
    !!
    subroutine set_omp_outer_thread_num( nthreads )
#ifdef _OPENMP
      use omp_lib
#endif
      implicit none
      integer :: nthreads
      logical :: is_omp_nested
      character(len=255) :: str

#ifndef _OPENMP
      nthreads=1
#else

      ! ----------------------
      ! Note omp_get_nested() may not work correctly
      ! with mpixlf95_r on BG/Q Mira
      ! ----------------------
      nthreads = omp_get_max_threads()
      nthreads = max(1, nthreads/col_f_nthreads )
      is_omp_nested = .false.
      str = ''
      call getenv("OMP_NESTED",str)
      if (len(trim(str)) >= 1) then
        if ((trim(str).eq."TRUE").or. &
            (trim(str).eq."true").or. &
            (trim(str).eq."1")) then
           is_omp_nested = .true.
        endif
      endif

      if (is_omp_nested) then
        if (sml_mype == 0) then
          !$omp  master
          write(*,'(a,L4)') 'omp_get_nested() ',omp_get_nested()
          write(*,'(a,i4)') 'omp_get_max_theads() ', omp_get_max_threads()
          write(*,'(a,i4)') 'omp_get_num_theads() ', omp_get_num_threads()
          write(*,'(a,i4)') 'nthreads = ', nthreads
          !$omp  end master
        endif
      endif
#endif

      return

    end subroutine set_omp_outer_thread_num


    !> This routine decides how the workload is distributed
    !! between threads --> How many grid points for each thread
    !> Allocate memory for Landau interaction matrix and
    !! the advection and diffusion operators
    !! This may be a lot of memory, especially for
    !! cases with multiple particle species. Therefore, they
    !! are allocated and deallocated with every call to the
    !! collision routine.
    !!
    !! @param[in]  isize   Dynamic memory size
    !!
    subroutine init_col_f_module_dynamic(isize)
      implicit none
      integer, intent(in) :: isize
      integer :: alloc_stat
      integer :: isp, ithread
    
      !!##### M_s, M_ab
      allocate(M_s( col_f_nvr-1, 5, (col_f_nvr-1)*(col_f_nvz-1),col_f_sp_s:col_f_sp_e,isize ),&
               stat=alloc_stat)

#ifdef _OPENACC
!$acc enter data pcreate(M_s)
#endif
      ! if((col_f_sp_num.ge.2).and.(opt_col.ne.2)) then
      if( ((col_f_sp_num.ge.2).or.(col_f_sp_s.ne.tar_f_sp_s)).and.(opt_col.ne.2)) then
        if(allocated(M_ab)) print *, 'BUG found. memory is allocated already'
        ! allocate(M_ab((col_f_nvr-1)*(col_f_nvz-1), 3, (col_f_nvr-1)*(col_f_nvz-1), &
        !                col_f_sp_s:col_f_sp_e-1, col_f_sp_s:col_f_sp_e, &
        !                isize ), stat=alloc_stat)
        ! youknow M_ab called M_ab(:,:,:,spj,spi)
        allocate(M_ab((col_f_nvr-1)*(col_f_nvz-1), 3, (col_f_nvr-1)*(col_f_nvz-1), &
                       tar_f_sp_s:tar_f_sp_e-1, col_f_sp_s:col_f_sp_e, &
                       isize ), stat=alloc_stat)                       
      else if (opt_col.eq.2) then
        if(allocated(M_ab)) print *, 'BUG found. memory is allocated already'
        allocate(M_ab((col_f_nvr-1)*(col_f_nvz-1), 3, (col_f_nvr-1)*(col_f_nvz-1), &
                       col_f_sp_s:col_f_sp_e, col_f_sp_s:col_f_sp_e, &
                       isize ), stat=alloc_stat)
        
#ifdef _OPENACC
!$acc enter data pcreate(M_ab)
#endif
      endif

      !!##### col_spall
      ! allocate(col_sp_cell_all(col_f_sp_s:col_f_sp_e, isize), stat=alloc_stat)
      allocate(col_sp_cell_all(tar_f_sp_s:tar_f_sp_e, isize), stat=alloc_stat)

      !!##### gammac_spall
      ! allocate(gammac_spall(col_f_sp_s:col_f_sp_e,col_f_sp_s:col_f_sp_e, isize), stat=alloc_stat)
      allocate(gammac_spall(tar_f_sp_s:tar_f_sp_e,col_f_sp_s:col_f_sp_e, isize), stat=alloc_stat)

      !!##### mesh_r, mesh_z, mesh_r_half, mesh_z_half, local_center_voume, pdf_n, pdf_np1, vol, delta_r, delta_z, ED
      do ithread=1,isize
        !  do isp=col_f_sp_s,col_f_sp_e
         do isp=tar_f_sp_s,tar_f_sp_e
             allocate(col_sp_cell_all(isp,ithread)%mesh_r(col_f_nvr), &
                      col_sp_cell_all(isp,ithread)%mesh_z(col_f_nvz), &
                      col_sp_cell_all(isp,ithread)%mesh_r_half(col_f_nvr-1), &
                      col_sp_cell_all(isp,ithread)%mesh_z_half(col_f_nvz-1), &
                      col_sp_cell_all(isp,ithread)%local_center_volume(col_f_nvr-1), &
                      col_sp_cell_all(isp,ithread)%pdf_n(col_f_nvr, col_f_nvz), &
                      col_sp_cell_all(isp,ithread)%pdf_np1(col_f_nvr, col_f_nvz), &
                      col_sp_cell_all(isp,ithread)%vol(col_f_nvr), &
                      col_sp_cell_all(isp,ithread)%power(tar_f_sp_s:tar_f_sp_e), & !JPL
                      col_sp_cell_all(isp,ithread)%col_dens(tar_f_sp_s:tar_f_sp_e), & !JPL
                      col_sp_cell_all(isp,ithread)%inner_info(tar_f_sp_s:tar_f_sp_e, 0:2, 1:col_iteration), & !youknow
                      ! col_sp_cell_all(isp,ithread)%delta_r(col_f_nvr-1,col_f_sp_s:col_f_sp_e), &
                      ! col_sp_cell_all(isp,ithread)%delta_z(col_f_nvz-1,col_f_sp_s:col_f_sp_e), &
                      col_sp_cell_all(isp,ithread)%delta_r(col_f_nvr-1,tar_f_sp_s:tar_f_sp_e), &
                      col_sp_cell_all(isp,ithread)%delta_z(col_f_nvz-1,tar_f_sp_s:tar_f_sp_e), &
                      col_sp_cell_all(isp,ithread)%ED(5, col_f_nvr-1,col_f_nvz-1),&
                      stat=alloc_stat)

              col_sp_cell_all(isp,ithread)%inner_info = 0D0
         enddo
      enddo

      return

    end subroutine init_col_f_module_dynamic


    !> Deallocate memory that is no longer needed after
    !! completion of the collision calculation
    !!
    !! @param[in]   isize   Dynamic memory size
    subroutine finalize_col_f_module_dynamic( isize )
      implicit none
      integer, intent(in) :: isize
      integer :: ithread, i
      integer :: alloc_stat
     
      deallocate(gammac_spall, stat=alloc_stat)

      do ithread=1,isize
        !  do i=col_f_sp_s,col_f_sp_e
         do i=tar_f_sp_s,tar_f_sp_e
            deallocate(col_sp_cell_all(i,ithread)%mesh_r,      col_sp_cell_all(i,ithread)%mesh_z,      &
                       col_sp_cell_all(i,ithread)%mesh_r_half, col_sp_cell_all(i,ithread)%mesh_z_half, &
                       col_sp_cell_all(i,ithread)%local_center_volume,&
                       col_sp_cell_all(i,ithread)%vol,col_sp_cell_all(i,ithread)%pdf_np1,&
                       col_sp_cell_all(i,ithread)%pdf_n,&
                       col_sp_cell_all(i,ithread)%delta_r, col_sp_cell_all(i,ithread)%delta_z, &
                       col_sp_cell_all(i,ithread)%ED,stat=alloc_stat)
         enddo
      enddo
      deallocate(col_sp_cell_all, stat=alloc_stat)
      deallocate(gammac_spall, stat=alloc_stat)

#ifdef _OPENACC
!$acc  exit data delete(M_s)
#endif
      deallocate(M_s,stat=alloc_stat)
      
      ! if(((col_f_sp_num.gt.1).and.(opt_col.ne.2)).or.(opt_col.eq.2)) then
      if (allocated(M_ab)) then
#ifdef _OPENACC
!$acc exit data delete(M_ab)
#endif
         deallocate(M_ab, stat=alloc_stat)
      endif

      return

    end subroutine finalize_col_f_module_dynamic


    !> Convert the distribution function to the units used
    !! internally by the collision operator
    !! The distribution function set up with f0_module is in normalized velocity space
    !! whereas the collision operator uses SI units.
    !!
    subroutine setup_col_sp(F_in,g_in,f0_M_tar_4d_in,inode, mass_in, t_ev_in,den_in,charge_in,&
                 col_spall, gammac,&
                 ith,ipsi)
      use FP4D_globals, only: sml_ev2j, sml_2pi, sml_pi, sml_e_charge, sml_prot_mass
      implicit none
      integer, intent(in) :: inode
      integer, intent(in) :: ith, ipsi
      real (kind=8), intent(in) :: F_in(-col_f_nvzm1/2:col_f_nvzm1/2,0:col_f_nthm1,  &
                                        0:col_f_nvrm1,col_f_sp_s:col_f_sp_e)
      real (kind=8), intent(in) :: g_in(-col_f_nvzm1/2:col_f_nvzm1/2,0:col_f_nthm1,  &
                                        0:col_f_nvrm1,col_f_sp_s:col_f_sp_e)
      real (kind=8), intent(in) :: f0_M_tar_4d_in(-col_f_nvzm1/2:col_f_nvzm1/2,0:col_f_nthm1,  &
                                        0:col_f_nvrm1,0         :col_f_sp_s)
      real (kind=8), dimension(0:col_f_sp_e) :: mass_in, t_ev_in, charge_in, den_in
      ! type(col_f_sp_type), intent(inout) :: col_spall(col_f_sp_s:col_f_sp_e)
      type(col_f_sp_type), intent(inout) :: col_spall(tar_f_sp_s:tar_f_sp_e)
      ! real (kind=8), intent(inout) :: gammac(col_f_sp_s:col_f_sp_e, col_f_sp_s:col_f_sp_e)
      real (kind=8), intent(inout) :: gammac(tar_f_sp_s:tar_f_sp_e, col_f_sp_s:col_f_sp_e)
      real (kind=8) :: mass, mesh_dz, mesh_dr, vth, smu_n, pi4bb0vth3_dd,  &
                       t_ev, den, charge, conv_factor

      integer :: isp, j, imu, ir, iz


      !! %PARAMETER SETTINGS BASED ON USER INPUT%
      ! do isp=col_f_sp_s,col_f_sp_e
      do isp=tar_f_sp_s,tar_f_sp_e

        col_spall(isp)%spid        = isp    !! this should be first for identification in any routine
        col_spall(isp)%nodeid      = inode  !! mainly for diagnosis in col_fm_core (redundant paramter in some sense)

        den    = den_in(isp)
        t_ev   = t_ev_in(isp)
        mass   = mass_in(isp)
        charge = charge_in(isp)
        conv_factor=1.D0/sqrt((t_ev*sml_2pi * sml_e_charge / mass)**3)

        col_spall(isp)%conv_factor = conv_factor
        do imu=0, col_f_nvrm1
          ! Local f with basic correction:
          ! Simply cut off all negative values
          
          if (isp.lt.col_f_sp_s) then ! isp < col_f_sp_s ==> only target species
            col_spall(isp)%pdf_n(imu+1,:)=f0_M_tar_4d_in(:,inode,imu,isp)*conv_factor
          else
            col_spall(isp)%pdf_n(imu+1,:)=max(F_in(:,inode,imu,isp),0.D0)*conv_factor
          endif
        enddo

        !! pdf_np1 = pdf_n
        do iz=1, col_f_nvz
          do ir=1, col_f_nvr
            col_spall(isp)%pdf_np1(ir,iz) = col_spall(isp)%pdf_n(ir,iz)
          enddo
        enddo
 
        !! vol
        vth=sqrt(t_ev*sml_ev2j/mass)
        col_spall(isp)%vth=vth
        pi4bb0vth3_dd=sml_2pi*vth*vth*vth*col_f_dsmu*col_f_dvp
        !! --- 0th component
        ! smu_n=col_f_dsmu/3D0
        smu_n=col_f_dsmu/3D0 ! youknow
        col_spall(isp)%vol(1) = 0.5D0*pi4bb0vth3_dd*smu_n
        !! --- 1st to (last-1)
        do imu=1, col_f_nvrm1-1
           smu_n=col_f_dsmu*real(imu,8)
           col_spall(isp)%vol(imu+1) = pi4bb0vth3_dd*smu_n
        enddo
        !! --- last component
        smu_n=col_f_dsmu*(real(col_f_nvrm1,8)-1D0/3D0) !rh imu or f0_nmu???
        col_spall(isp)%vol(col_f_nvr)=0.5D0*pi4bb0vth3_dd*smu_n
        !!vth_par  = sqrt(sp_t_par(i)  * sml_e_charge / mass)
        !!vth_perp = sqrt(sp_t_perp(i) * sml_e_charge / mass)
        mesh_dz = col_f_dvp *vth
        mesh_dr = col_f_dsmu*vth

        col_spall(isp)%T_avg            = t_ev
        col_spall(isp)%den_lambda_gamma = den
        col_spall(isp)%mass_au          = mass/sml_prot_mass
        col_spall(isp)%mass             = mass
        col_spall(isp)%charge_eu        = charge/sml_e_charge
        col_spall(isp)%Ze               = charge*sml_e_charge
        !! < Mesh-related quantities>
        col_spall(isp)%vpar_beg         = -col_f_vp_max*vth  !! (=lx)
        col_spall(isp)%mesh_dz          = mesh_dz      !! (=mesh_dz)
        col_spall(isp)%mesh_dr          = mesh_dr     !! (=mesh_dr)


        ! mesh_r, mesh_r_half, mesh_z, mesh_z_half
        do j=1, col_f_nvr-1
          col_spall(isp)%mesh_r(j)       = mesh_dr*(real(j,8)-1.0)
          col_spall(isp)%mesh_r_half(j)  = mesh_dr*(real(j,8)-0.5)
        end do
        col_spall(isp)%mesh_r(j)       = mesh_dr*(real(j,8)-1.0)
        do j=1, col_f_nvz-1
          col_spall(isp)%mesh_z(j)       = col_spall(isp)%vpar_beg + mesh_dz*(real(j,8)-1.0)
          col_spall(isp)%mesh_z_half(j)  = col_spall(isp)%vpar_beg + mesh_dz*(real(j,8)-0.5)
        end do
        col_spall(isp)%mesh_z(col_f_nvz)   = col_spall(isp)%vpar_beg + mesh_dz*(col_f_nvz-1)
        !! local_center_volume (mesh_r_half, mesh_dr, mesh_dz)
        !! : volume centered at a cell
        col_spall(isp)%local_center_volume = col_spall(isp)%mesh_r_half*mesh_dr*mesh_dz

     enddo !isp

      !gamma
      do isp=col_f_sp_s, col_f_sp_e
        !  do j=col_f_sp_s, col_f_sp_e
         do j=tar_f_sp_s, tar_f_sp_e

            call col_f_lambda_gamma_pair(col_spall(isp), col_spall(j), gammac(j,isp), ipsi)  !! C(i->j)

         enddo !j
      enddo

      return

    end subroutine setup_col_sp


    !> Computes the base collision frequency -- Eq. (6) in Hager et al.
    !!
    subroutine col_f_lambda_gamma_pair(sp_a, sp_b, gammac_loc, ipsi)  !! C(a->b)
    ! colliding : a, target : b
      use FP4D_globals, only : sml_e_charge,  sml_pi
      use FP4D_globals, only : opt_time_solver, opt_time_species
      use FP4D_globals, only : eq_tau,sml_istep
      use FP4D_globals, only : opt_logaritm_NEO, eq_den, eq_t_ev
      

      implicit none
      integer, intent(in) :: ipsi
      type(col_f_sp_type) :: sp_a, sp_b
      real (kind=8) :: gammac_loc
      !
      real (8) :: lambda, ti_ev, massi, densi, chargei, te_ev, masse, dense
      real (8) :: tmpr
      logical :: e_in, e_col, e_tar

      !! check if electron is in
      e_col = sp_a%spid .eq. 0
      e_tar = sp_b%spid .eq. 0
      e_in  = e_col .or. e_tar

      !
      if ( (opt_time_solver.eq.-1) .and. (opt_time_species.eq.99) ) then
        write(*,*) 'error :: set the opt_time_species'
      endif

      select case(opt_logaritm_NEO)
        case(0) ! default : following NRL Plasma Formulary 2013
          ! NRL Plasma Formulary 2011 version, p.34
          !(c) Mixed ion-ion collisions (here, same species ion-ion collisions)
          if(e_in) then
            !if electron is invovled
              if( e_col .and. e_tar ) then
                  !(a) Thermal electron-electron collisions
                  ! Note that pointers for colliding and target are same.
                  lambda = 2.35D1 - log(sqrt(sp_a%den_lambda_gamma*1D-6)*(sp_a%T_avg**(-1.25D0)))  &
                                  - sqrt(abs(1D-5+((log(sp_a%T_avg) -2D0)**2)/16D0))
                  gammac_loc = (sml_e_charge**4) * lambda / (sp_a%mass*((8.8542D-12)**2) * 8D0 * sml_pi)
              else
                  if(e_col) then
                    ti_ev   = sp_b%T_avg
                    massi   = sp_b%mass_au
                    densi   = sp_b%den_lambda_gamma
                    chargei = sp_b%charge_eu
                    te_ev   = sp_a%T_avg
                    masse   = sp_a%mass_au
                    dense   = sp_a%den_lambda_gamma
                  else
                    te_ev   = sp_b%T_avg
                    masse   = sp_b%mass_au
                    dense   = sp_b%den_lambda_gamma
                    ti_ev   = sp_a%T_avg
                    massi   = sp_a%mass_au
                    densi   = sp_a%den_lambda_gamma
                    chargei = sp_a%charge_eu
                  endif
                  !(b) Electron-ion (and ion-electron) collisions
                  if(ti_ev * masse/massi .lt. te_ev) then
                    if(te_ev .lt. 1D1*(chargei**2)) then
                        lambda = 2.3D1 - log(chargei*sqrt(dense*1D-6/(te_ev**3)))
                    else
                        lambda = 2.4D1 - log(sqrt(dense*1D-6)/te_ev)
                    endif
                  else
                    lambda = 3D1 - log((chargei**2)*sqrt(densi*1D-6/(ti_ev**3))/massi)
                  endif
          
                  tmpr = (sml_e_charge**4) * (chargei**2) * lambda /( ((8.8542D-12)**2) * 8D0 * sml_pi )
                  gammac_loc = tmpr / sp_a%mass ! colliding : a, target : b
              endif
          else ! if(e_in) then
            lambda = 2.3D1 - log( sp_a%charge_eu * sp_b%charge_eu * (sp_a%mass_au + sp_b%mass_au) &
                                      / (sp_a%mass_au*sp_b%T_avg + sp_b%mass_au*sp_a%T_avg) &
                                      * sqrt( sp_a%den_lambda_gamma*1D-6*(sp_a%charge_eu**2)/sp_a%T_avg &
                                            + sp_b%den_lambda_gamma*1D-6*(sp_b%charge_eu**2)/sp_b%T_avg) )
            tmpr = (sml_e_charge**4) * (sp_a%charge_eu**2) * (sp_b%charge_eu**2) * lambda /( ((8.8542D-12)**2) * 8D0 * sml_pi )
            gammac_loc = tmpr / sp_a%mass   ! colliding : a, target : b
          endif ! if(e_in) then

        case(1) ! following NEO (all species has ln Lambda_ei)

          ! 0th species must be electron
          lambda = 2.4D1 - log(sqrt(eq_den(ipsi,0)*1D-6)/eq_t_ev(ipsi,0))
          tmpr = (sml_e_charge**4) * (sp_a%charge_eu**2) * (sp_b%charge_eu**2) * lambda /( ((8.8542D-12)**2) * 8D0 * sml_pi )
          gammac_loc = tmpr / sp_a%mass   ! colliding : a, target : b          

          ! write(*,*) 'lambda', sp_a%spid, sp_b%spid, lambda

      end select

      ! write(*,*) 

      ! write(*,*) 'sp_a%spid', sp_a%spid
      ! write(*,*) 'sp_b%spid', sp_b%spid
      ! write(*,*) eq_tau(ipsi,opt_time_species,sp_a%spid)

      if (opt_time_solver.eq.-1) then
        ! colliding particle have to have the same collision time. --> fix the 'b' when define the tau_ab
        gammac_loc = gammac_loc * col_f_dt * eq_tau(ipsi,opt_time_species,sp_a%spid) ! eq_tau(psi, target, colliding)
      else
        gammac_loc = gammac_loc * col_f_dt
      endif

    end subroutine col_f_lambda_gamma_pair


    !> Top-level collision routine
    !! This routine performs the loop over configuration space grid points with
    !! OpenMP parallelization --> "Outer level" parallelism.
    !!
    subroutine fm_collision(f_in,g1_in,dfc_out,nthreads,&
                 sml_istep,int_n,int_T,int_p,FSAn,FSAT,FSAp,&
                 qld_B,qld_C,qld_F,f0_M_tar_in,col_power_out,col_dens_out,col_inner_info_out)
      use FP4D_globals, only: sml_nstep  !JPL

      use FP4D_globals, only: sml_totalpe

      use FP4D_globals, only: nptot_mpi, ip1_mpi,  npsi, nth, nvr, nvz, root
      ! youknow
      use FP4D_globals, only:sml_mass,sml_charge,RF_kpar
      use FP4D_globals, only: eq_den,eq_t_ev,eq_J,eq_intJ

      use FP4D_globals, only:is_s, is_e, tar_s, tar_e, ip_s, ip_e
      use FP4D_globals, only:ith_s, ith_e, ir_s, ir_e, iz_s, iz_e
#ifdef _OPENMP
      use omp_lib
#endif
      use mpi
      implicit none
      !
      ! youknow_TORIC : caution for order of dimension
      real (kind=8), intent(in) :: qld_B(1:col_f_npsi, 0:col_f_nthm1, 1:col_f_nvrm1, 1:col_f_nvzm1)
      real (kind=8), intent(in) :: qld_C(1:col_f_npsi, 0:col_f_nthm1, 1:col_f_nvrm1, 1:col_f_nvzm1)
      real (kind=8), intent(in) :: qld_F(1:col_f_npsi, 0:col_f_nthm1, 1:col_f_nvrm1, 1:col_f_nvzm1)
      
      integer, intent(in) :: nthreads,sml_istep
      real(kind=8), intent(in) :: f_in       (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)

      real(kind=8), intent(in) :: g1_in      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)
      real(kind=8), intent(in) :: f0_M_tar_in(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, tar_s:tar_e) ! 0:isp
                                                                    
      real(kind=8), intent(out) :: dfc_out   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)
      real(kind=8), intent(out) :: col_power_out (ith_s:ith_e, ip_s:ip_e, is_s:is_e, tar_s:tar_e)! youknow
      real(kind=8), intent(out) :: col_dens_out  (ith_s:ith_e, ip_s:ip_e, is_s:is_e, tar_s:tar_e) ! youknow
      real(kind=8), intent(out) :: col_inner_info_out(ith_s:ith_e, ip_s:ip_e, is_s:is_e, tar_s:tar_e, 0:2, 1:col_iteration) ! youknow


      real(kind=8), dimension(iz_s:iz_e, ith_s:ith_e, ir_s:ir_e, is_s:is_e) :: f,g1 !!! youknow :: order is different with f_in, g1_in
      real(kind=8), dimension(iz_s:iz_e, ith_s:ith_e, ir_s:ir_e, tar_s:tar_e) :: f0_M_tar_4d  ! 0:isp   !!! youknow :: order is different with f_in, g1_in

      real(kind=8), intent(out),dimension(ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: int_n, int_T,int_p
      real(kind=8), intent(out),dimension(ip_s:ip_e, is_s:is_e) :: FSAn, FSAT,FSAp 

      ! Hmm... adiabatic FSAn ? no!
      integer :: alloc_stat
      integer :: inode1, inode2, stride, isize, i_stride
      integer :: i_beg(size(col_f_mat_list)), i_end(size(col_f_mat_list))
      integer :: ipsi, ith, isp, imu, jsp
      integer :: ierr, ipe, ndata, tag,count,elements, i_1, i_2  !JPL
      integer :: ntemp_p, ip1
      integer :: newtype, newtype1, newtype2, newtype3
      integer :: newtype_send, newtype_recv

      integer, dimension(4) :: sizes, subsizes, starts
      integer, dimension(2) :: sizes2, subsizes2, starts2
      integer rstatus(MPI_STATUS_SIZE)
      logical :: update_ok


      ! write(*,*) 'size check'
      ! write(*,*) 'col_f_nvzm1', col_f_nvzm1
      ! write(*,*) 'col_f_nvrm1', col_f_nvrm1
      ! write(*,*) 'col_f_nthm1', col_f_nthm1
      ! write(*,*) 'col_f_npsi', col_f_npsi

      FSAn=0.0; FSAT=0.0; FSAp=0.0;
      dfc_out(:,:,:,:,:) = 0D0;
!      call get_mesh_node_range_for_this_mpirank(inode1_in, inode2_in, inode1, inode2, stride)
!      call get_mesh_node_range_for_threads_in_this_mpirank ( inode1, inode2, stride,  &
!                                        nthreads, isize, i_beg, i_end, i_stride )
      
!     JPL: Dec 05, 2024. use thread only on vperp and vpar, and mpi only on psi
!     and theta parallelization --> col_f_i_nthreads=col_f_o_nthreads, iszie=1

!     call get_mesh_node_range_for_this_mpirank(0, col_f_nthm1, inode1, inode2,
!     stride)
      inode1=1  !JPL
      inode2=1  !JPL
      stride=1
      call get_mesh_node_range_for_threads_in_this_mpirank ( inode1, inode2,stride,  &
                                        nthreads, isize, i_beg, i_end, i_stride)


      isize=col_f_npsi   !JPL:: isize is defined by npsi?


       call init_col_f_module_dynamic( isize )   !! allocating memory for M_s (and M_ss)
      
      !$omp parallel default (none)                                          &
      !$omp shared(FSAn, FSAT,FSAp, int_n, int_T, int_p) &
      !$omp shared(isize, sml_mype, col_f_npsi, sml_istep, eq_J, eq_intJ, col_f_dth) &
      !$omp shared(col_sp_cell_all, col_f_mat_list, col_f_ksp_list, col_f_vecb_list)&
      !$omp shared(col_f_vecx_list, col_f_sp_num, gammac_spall, M_s, M_ab)&
      !$omp shared(col_f_sp_s, col_f_sp_e, col_f_nvrm1,  f_in)&
      !$omp shared(tar_f_sp_s, tar_f_sp_e) &
      !$omp shared(sml_mass, sml_charge, eq_den, eq_t_ev) &
      !$omp shared(dfc_out, col_f_nthm1, f, g1, opt_Col, g1_in)&
      !$omp shared(f0_M_tar_in, f0_M_tar_4d, col_power_out, col_dens_out,  col_inner_info_out)             &
      !$omp shared(qld_B, qld_C, qld_F,rf_pwrscale) &
      !$omp shared(istart1,iend1,istart2,iend2) &  !JPL
      !$omp shared(col_iteration, i_1, i_2) &  !youknow
      !$omp private(ith, ipsi, update_ok, isp, imu, jsp) num_threads(col_f_o_nthreads)
      !$omp do
      do ipsi=istart1(sml_mype),iend1(sml_mype)    !  JPL: mpi for psi and th
        do ith=istart2(sml_mype)-1,iend2(sml_mype)-1

          do imu=0, col_f_nvrm1
            ! do isp = col_f_sp_s, col_f_sp_e
            do isp = tar_f_sp_s, tar_f_sp_e ! for define the f0_M_tar

              if (isp.lt.col_f_sp_s) then ! Array size is [0:isp] but [0:isp-1] 
                ! f0_M_tar_4d(:,ith,imu,isp)=f0_M_tar_in(:,ith,imu,ipsi,isp)
                f0_M_tar_4d(:,ith,imu,isp) = f0_M_tar_in(:,imu,ith,ipsi,isp) ! caution for order
              else
                ! f(:,ith,imu,isp)=f_in(:,ith,imu,ipsi,isp)
                ! g1(:,ith,imu,isp)=g1_in(:,ith,imu,ipsi,isp) 
                f (:,ith,imu,isp) = f_in (:,imu,ith,ipsi,isp)  ! caution for order
                g1(:,ith,imu,isp) = g1_in(:,imu,ith,ipsi,isp)  ! caution for order
              endif
            enddo  !isp
          enddo   !imu
        !enddo

        !do ith=0,col_f_nthm1  !JPL
          call setup_col_sp(f,g1,f0_M_tar_4d,ith, sml_mass, eq_t_ev(ipsi,:), eq_den(ipsi,:),sml_charge,&
                            col_sp_cell_all(:,ipsi), gammac_spall(:,:,ipsi),&
                            ith,ipsi)

           !call MPI_Barrier(MPI_COMM_WORLD, ierr)  !JPL
        !!##### TIME ADVANCE
          select case (opt_col)
            case(0,1)
              ! if (col_f_sp_num.eq.1) then !!## Single Species Collisions
              if ( (col_f_sp_num.eq.1) .and. (col_f_sp_s.eq.0) ) then !!## Single Species Collisions
                  call col_fm_core(ith,ipsi,col_sp_cell_all(:,ipsi),gammac_spall(:,:,ipsi),update_ok,sml_istep,&
                  qld_B(ipsi,ith,:,:), qld_C(ipsi,ith,:,:), qld_F(ipsi,ith,:,:), &       ! qld_B (1:npsi,0:nthm1,0:nmu,-nvp:nvp)
                  col_f_mat_list(ipsi),col_f_ksp_list(ipsi),col_f_vecb_list(ipsi),col_f_vecx_list(ipsi),M_s(:,:,:,:,ipsi))
              else !!## Multi-Species Collisions
                call col_fm_core(ith,ipsi,col_sp_cell_all(:,ipsi),gammac_spall(:,:,ipsi),update_ok,sml_istep,&
                  qld_B(ipsi,ith,:,:), qld_C(ipsi,ith,:,:), qld_F(ipsi,ith,:,:), &       ! qld_B (1:npsi,0:nthm1,0:nmu,-nvp:nvp)
                  col_f_mat_list(ipsi),col_f_ksp_list(ipsi),col_f_vecb_list(ipsi),col_f_vecx_list(ipsi),M_s(:,:,:,:,ipsi),M_ab(:,:,:,:,:,ipsi))
              end if
            case(2)  !!## Relativistic Collisions
              call col_fm_core(ith,ipsi,col_sp_cell_all(:,ipsi),gammac_spall(:,:,ipsi),update_ok,sml_istep,&
                 qld_B(ipsi,ith,:,:), qld_C(ipsi,ith,:,:), qld_F(ipsi,ith,:,:), &       ! qld_B (1:npsi,0:nthm1,0:nmu,-nvp:nvp)
                 col_f_mat_list(ipsi),col_f_ksp_list(ipsi),col_f_vecb_list(ipsi),col_f_vecx_list(ipsi),M_s(:,:,:,:,ipsi),M_ab(:,:,:,:,:,ipsi))
          end select
          do isp = col_f_sp_s, col_f_sp_e
            int_n(ith,ipsi,isp)=col_sp_cell_all(isp,ipsi)%dens
            int_T(ith,ipsi,isp)=col_sp_cell_all(isp,ipsi)%numeric_T
            int_p(ith,ipsi,isp)=col_sp_cell_all(isp,ipsi)%mom
            if (update_ok.and.(opt_col.ne.0)) then
              do imu=0, col_f_nvrm1
                dfc_out(:,imu,ith,ipsi,isp) = col_sp_cell_all(isp,ipsi)%pdf_np1(imu+1,:)
              enddo
            endif
            ! youknow :: Save 'col_dw' for each species 
            do jsp = tar_f_sp_s, tar_f_sp_e ! Caution :: colliding species [type var] has each target power
              col_power_out(ith,ipsi,isp,jsp) = col_sp_cell_all(isp,ipsi)%power(jsp)   !JPL: reorder col_power_out 
              col_dens_out(ith,ipsi,isp,jsp) = col_sp_cell_all(isp,ipsi)%col_dens(jsp)   !JPL: reorder col_power_out 

              do i_1 = 0, 2
              do i_2 = 1, col_iteration
                  col_inner_info_out(ith,ipsi,isp,jsp,i_1,i_2) = col_sp_cell_all(isp,ipsi)%inner_info(jsp,i_1,i_2)   !JPL: reorder col_power_out 
              enddo
              enddo

              ! write(*,*) ''
              ! write(*,*) 'STEP 2'
              ! write(*,*) 'ith isp jsp',ith,isp,jsp,'col_inner_info_out-0',col_inner_info_out(ith,ipsi,isp,jsp,0,:),col_dens_out(ith,ipsi,isp,jsp)
              ! write(*,*) 'ith isp jsp',ith,isp,jsp,'col_inner_info_out-1',col_inner_info_out(ith,ipsi,isp,jsp,1,:),'00000'
              ! write(*,*) 'ith isp jsp',ith,isp,jsp,'col_inner_info_out-2',col_inner_info_out(ith,ipsi,isp,jsp,2,:),col_power_out(ith,ipsi,isp,jsp)

            enddo
          enddo ! isp
        enddo   ! ith
      end do !JPL: ipsi
 
!$omp end do
!$acc wait
!$omp end parallel

!$acc wait



!       JPL: MPI transfer
!       send to root

      tag    = 2001
      sizes    = [col_f_nvz,col_f_nvr,nth,npsi]
      starts   = [0,0,0,0]     ! let's say we're looking at region "0"  

      ntemp_p=1
      do ip1=1,nptot_mpi
        ntemp_p=ntemp_p*(sml_mype-ip1_mpi(ip1))
      enddo
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! write(*,*) 'works 0',sizes
      if ((sml_mype .ne. root).and.(ntemp_p.eq.0)) then
         subsizes = [col_f_nvz,col_f_nvr,iend2(sml_mype)-istart2(sml_mype)+1,iend1(sml_mype)-istart1(sml_mype)+1]
! size of sub-region
         call MPI_Type_create_subarray(4, sizes, subsizes, starts,MPI_ORDER_FORTRAN, MPI_double, newtype, ierr)
              ! write(*,*) 'works 00',sml_mype,subsizes
         call MPI_Type_commit(newtype, ierr)
!         ndata=nvr*nvz*(iend2(sml_mype)-istart2(sml_mype)+1)*(iend1(sml_mype)-istart1(sml_mype)+1)

         do  isp = col_f_sp_s, col_f_sp_e
                ! write(*,*) 'works 01',sml_mype,istart2(sml_mype)-1,istart1(sml_mype)
            call MPI_Send(dfc_out(-col_f_nvzm1/2,0,istart2(sml_mype)-1,istart1(sml_mype),isp), &
                 1, newtype, root, tag, MPI_COMM_WORLD, ierr)
                ! write(*,*) 'works 011',sml_mype
         enddo
      else if (sml_mype .eq. root) then
         do ip1=2,nptot_mpi
           ipe=ip1_mpi(ip1)
              ! write(*,*) 'works 02',sml_mype,ipe
          !  ndata=nvr*nvz*(iend2(ipe)-istart2(ipe)+1)*(iend1(ipe)-istart1(ipe)+1)
           subsizes = [col_f_nvz,col_f_nvr,iend2(ipe)-istart2(ipe)+1,iend1(ipe)-istart1(ipe)+1]
            ! write(*,*) 'works 03',sml_mype,subsizes
           call MPI_Type_create_subarray(4, sizes, subsizes,starts,MPI_ORDER_FORTRAN,MPI_double, newtype, ierr)
           call MPI_Type_commit(newtype, ierr)
            ! write(*,*) 'works 04',sml_mype,subsizes
           do  isp = col_f_sp_s, col_f_sp_e
                !  write(*,*) 'works 50',sml_mype,isp, ipe, istart2(ipe),dfc_out(-col_f_nvzm1/2,0,10,1,0)
                call MPI_Recv(dfc_out(-col_f_nvzm1/2,0,istart2(ipe)-1,istart1(ipe),isp), &
                     1, newtype, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)

                !  write(*,*) 'works51',sml_mype,isp,dfc_out(-col_f_nvzm1/2,0,10,1,0)
           enddo
         enddo
      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

 do ipsi=1,isize    ! isize = outer thread num?????
        do ith=0,col_f_nthm1
             do isp = col_f_sp_s, col_f_sp_e
     
      ndata=col_f_nvz*col_f_nvr
      !write(*,*) 'works 91',sml_mype, ndata, ipsi,ith,isp

      call MPI_BCAST(dfc_out(-col_f_nvzm1/2,0,ith,ipsi,isp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
      !write(*,*) 'works 92',sml_mype, ndata, ipsi,ith,isp
enddo
enddo
enddo


! !     transfer dens..


      tag    = 2002
      sizes2    = [nth,npsi]
      starts2   = [0,0]     ! let's say we're looking at region "0"  

      ntemp_p=1
      do ip1=1,nptot_mpi
        ntemp_p=ntemp_p*(sml_mype-ip1_mpi(ip1))
      enddo
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! write(*,*) 'works 1',sizes
      if ((sml_mype .ne. root).and.(ntemp_p.eq.0)) then
         subsizes2 =[iend2(sml_mype)-istart2(sml_mype)+1,iend1(sml_mype)-istart1(sml_mype)+1] ! size of sub-region
         call MPI_Type_create_subarray(2, sizes2, subsizes2,starts2,MPI_ORDER_FORTRAN, MPI_double, newtype2, ierr)
              ! write(*,*) 'works 10',sml_mype,subsizes2
         call MPI_Type_commit(newtype2, ierr)
!         ndata=nvr*nvz*(iend2(sml_mype)-istart2(sml_mype)+1)*(iend1(sml_mype)-istart1(sml_mype)+1)

         do  isp = col_f_sp_s, col_f_sp_e
                  ! write(*,*) 'works 11',sml_mype,istart2(sml_mype)-1,istart1(sml_mype)
              call MPI_Send(int_n(istart2(sml_mype)-1,istart1(sml_mype),isp), &
                  1, newtype2, root, tag, MPI_COMM_WORLD, ierr)
              call MPI_Send(int_T(istart2(sml_mype)-1,istart1(sml_mype),isp), &
                  1, newtype2, root, tag, MPI_COMM_WORLD, ierr)
              call MPI_Send(int_P(istart2(sml_mype)-1,istart1(sml_mype),isp), &
                  1, newtype2, root, tag, MPI_COMM_WORLD, ierr)

              do jsp = tar_f_sp_s, tar_f_sp_e ! Caution :: colliding species [type var] has each target power
                  call MPI_Send(col_power_out(istart2(sml_mype)-1,istart1(sml_mype),isp, jsp), &  !JPL: reorder col_power_out  
                      1, newtype2, root, tag, MPI_COMM_WORLD, ierr)
                  call MPI_Send(col_dens_out(istart2(sml_mype)-1,istart1(sml_mype),isp, jsp), &  !JPL: reorder col_power_out  
                      1, newtype2, root, tag, MPI_COMM_WORLD, ierr)

                  do i_1 = 0, 2
                  do i_2 = 1, col_iteration
                      call MPI_Send(col_inner_info_out(istart2(sml_mype)-1,istart1(sml_mype),isp, jsp, i_1, i_2), &  !JPL: reorder col_power_out  
                          1, newtype2, root, tag, MPI_COMM_WORLD, ierr)
                  enddo
                  enddo
              enddo
         enddo
      else if (sml_mype .eq. root) then
         do ip1=2,nptot_mpi
           ipe=ip1_mpi(ip1)
              ! write(*,*) 'works 12',sml_mype,ipe
           !  ndata=nvr*nvz*(iend2(ipe)-istart2(ipe)+1)*(iend1(ipe)-istart1(ipe)+1)
           subsizes2 =[iend2(ipe)-istart2(ipe)+1,iend1(ipe)-istart1(ipe)+1]
            ! write(*,*) 'works 13',sml_mype,subsizes2
           call MPI_Type_create_subarray(2, sizes2, subsizes2,starts2,MPI_ORDER_FORTRAN,MPI_double, newtype2, ierr)
           call MPI_Type_commit(newtype2, ierr)
            ! write(*,*) 'works 14',sml_mype,subsizes2
           do  isp = col_f_sp_s, col_f_sp_e
                !  write(*,*) 'works 15',sml_mype,isp, ipe, istart2(ipe)
                call MPI_Recv(int_n(istart2(ipe)-1,istart1(ipe),isp), &
                     1, newtype2, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)
                call MPI_Recv(int_T(istart2(ipe)-1,istart1(ipe),isp), &
                     1, newtype2, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)
                call MPI_Recv(int_P(istart2(ipe)-1,istart1(ipe),isp), &
                     1, newtype2, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)
                !  write(*,*) 'works 16',sml_mype,isp

                do jsp = tar_f_sp_s, tar_f_sp_e ! Caution :: colliding species [type var] has each target power

                    call MPI_Recv(col_power_out(istart2(ipe)-1,istart1(ipe),isp,jsp), &
                              1, newtype2, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)           
                    call MPI_Recv(col_dens_out(istart2(ipe)-1,istart1(ipe),isp,jsp), &
                              1, newtype2, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)   

                    do i_1 = 0, 2
                    do i_2 = 1, col_iteration
                        call MPI_Recv(col_inner_info_out(istart2(ipe)-1,istart1(ipe),isp,jsp,i_1,i_2), &
                             1, newtype2, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)   
                    enddo
                    enddo
                enddo
           enddo  !isp
         enddo  !ip1
      endif

call MPI_Barrier(MPI_COMM_WORLD, ierr)
 do ipsi=1,isize    ! isize = outer thread num?????
    do isp = col_f_sp_s, col_f_sp_e

        ndata=nth
        call MPI_BCAST(int_n(0,ipsi,isp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(int_T(0,ipsi,isp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(int_p(0,ipsi,isp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
        do jsp = tar_f_sp_s, tar_f_sp_e 
            call MPI_BCAST(col_power_out(0,ipsi,isp,jsp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(col_dens_out(0,ipsi,isp,jsp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)

            do i_1 = 0, 2
            do i_2 = 1, col_iteration
                call MPI_BCAST(col_inner_info_out(0,ipsi,isp,jsp,i_1,i_2), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
            enddo
            enddo

        enddo
    enddo
enddo


      call MPI_Barrier(MPI_COMM_WORLD, ierr)
     ! ndata=nth*npsi*(col_f_sp_e-col_f_sp_s+1)
     ! call MPI_BCAST(int_n, ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)

      !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !  write(*,*) 'int_n1',sml_mype, int_n(0,1,0)
      ! write(*,*) 'int_n22',sml_mype, int_n(10,1,0)
      !call MPI_BCAST(int_p, ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)

      !call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !call MPI_BCAST(int_T, ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
      ! write(*,*) 'works 93'

      !call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !ndata=nth*npsi*(col_f_sp_e-col_f_sp_s+1)*(tar_f_sp_e-tar_f_sp_s+1)
      !call MPI_BCAST(col_power_out, ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
      !write(*,*) 'works 94'
     !call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! write(*,*) 'dfc_out1',sml_mype, dfc_out(-col_f_nvzm1/2,0,0,1,0)
      ! write(*,*) 'dfc_out2',sml_mype, dfc_out(-col_f_nvzm1/2,0,10,1,0)
!       post process


      ! do ipsi=1,npsi
      ! do ith=0,col_f_nthm1
      ! do isp = col_f_sp_s,col_f_sp_e
      ! do jsp = tar_f_sp_s,tar_f_sp_e
      !       write(*,*) ''
      !       write(*,*) 'STEP 3'
      !       write(*,*) 'ith,isp,jsp',ith,isp,jsp,'col_dens ',col_dens_out (ith,ipsi,isp,jsp)
      !       write(*,*) 'ith,isp,jsp',ith,isp,jsp,'col_power',col_power_out(ith,ipsi,isp,jsp)
      !       write(*,*) 'ith,isp,jsp',ith,isp,jsp,'col_inner_info_out 0',col_inner_info_out(ith,ipsi,isp,jsp,0,:)
      !       write(*,*) 'ith,isp,jsp',ith,isp,jsp,'col_inner_info_out 1',col_inner_info_out(ith,ipsi,isp,jsp,1,:)
      !       write(*,*) 'ith,isp,jsp',ith,isp,jsp,'col_inner_info_out 2',col_inner_info_out(ith,ipsi,isp,jsp,2,:)
      !       write(*,*) ''

      ! enddo
      ! enddo
      ! enddo
      ! enddo
      

       do ipsi=1,npsi
        do   isp = col_f_sp_s,col_f_sp_e
          do ith=0,col_f_nthm1
            FSAn(ipsi,isp)=FSAn(ipsi,isp)+int_n(ith,ipsi,isp)*eq_J(ith,ipsi)*col_f_dth
            FSAp(ipsi,isp)=FSAp(ipsi,isp)+int_p(ith,ipsi,isp)*eq_J(ith,ipsi)*col_f_dth
            FSAT(ipsi,isp)=FSAT(ipsi,isp)+int_T(ith,ipsi,isp)*eq_J(ith,ipsi)*col_f_dth
          enddo
          FSAn(ipsi,isp)=FSAn(ipsi,isp)/eq_intJ(ipsi)
          FSAT(ipsi,isp)=FSAT(ipsi,isp)/eq_intJ(ipsi)
          FSAp(ipsi,isp)=FSAp(ipsi,isp)/eq_intJ(ipsi)
        enddo
      enddo


      call finalize_col_f_module_dynamic( isize )
      if (sml_mype.eq.0) then
        do isp=col_f_sp_s,col_f_sp_e
          write(100+isp,'(i8,6e20.12)') sml_istep,FSAn(:,isp)
          write(200+isp,'(i8,6e20.12)') sml_istep,FSAT(:,isp)
          write(300+isp,'(i8,6e20.12)') sml_istep,FSAp(:,isp)
        enddo
      endif 
!rh unnecessary
!#ifdef _OPENMP
!      !  ---------------------
!      !  reset omp_num_threads
!      !  ---------------------
!      nthreads = omp_get_max_threads()
!      call omp_set_num_threads( nthreads )
!#endif

    end subroutine fm_collision

    ! fm_QLD from fm_collision
    subroutine fm_QLD(f_in,g1_in,df_QLD_out,nthreads,&                 
                 qld_B,qld_C,qld_F,f0_M_tar_in)

      use FP4D_globals, only: sml_nstep  !JPL

      use FP4D_globals, only: sml_totalpe

      use FP4D_globals, only: nptot_mpi, ip1_mpi,  npsi, nth, nvr, nvz, root
      ! youknow
      use FP4D_globals, only:sml_mass,sml_charge,RF_kpar
      use FP4D_globals, only: eq_den,eq_t_ev

      use FP4D_globals, only:is_s, is_e, tar_s, tar_e, ip_s, ip_e
      use FP4D_globals, only:ith_s, ith_e, ir_s, ir_e, iz_s, iz_e
#ifdef _OPENMP
      use omp_lib
#endif
      use mpi
      implicit none
      !
      ! youknow_TORIC : caution for order of dimension
      real (kind=8), intent(in) :: qld_B(1:col_f_npsi, 0:col_f_nthm1, 1:col_f_nvrm1, 1:col_f_nvzm1)
      real (kind=8), intent(in) :: qld_C(1:col_f_npsi, 0:col_f_nthm1, 1:col_f_nvrm1, 1:col_f_nvzm1)
      real (kind=8), intent(in) :: qld_F(1:col_f_npsi, 0:col_f_nthm1, 1:col_f_nvrm1, 1:col_f_nvzm1)
      
      integer, intent(in) :: nthreads
      real(kind=8), intent(in) :: f_in       (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)

      real(kind=8), intent(in) :: g1_in      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)
      real(kind=8), intent(in) :: f0_M_tar_in(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, tar_s:tar_e) ! 0:isp
                                                                    
      real(kind=8), intent(out) :: df_QLD_out   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)


      real(kind=8), dimension(iz_s:iz_e, ith_s:ith_e, ir_s:ir_e, is_s:is_e) :: f,g1 !!! youknow :: order is different with f_in, g1_in
      real(kind=8), dimension(iz_s:iz_e, ith_s:ith_e, ir_s:ir_e, tar_s:tar_e) :: f0_M_tar_4d  ! 0:isp   !!! youknow :: order is different with f_in, g1_in

      ! Hmm... adiabatic FSAn ? no!
      integer :: alloc_stat
      integer :: inode1, inode2, stride, isize, i_stride
      integer :: i_beg(size(col_f_mat_list)), i_end(size(col_f_mat_list))
      integer :: ipsi, ith, isp, imu, jsp
      integer :: ierr, ipe, ndata, tag,count,elements, i_1, i_2  !JPL
      integer :: ntemp_p, ip1
      integer :: newtype, newtype1, newtype2, newtype3
      integer :: newtype_send, newtype_recv

      integer, dimension(4) :: sizes, subsizes, starts
      integer, dimension(2) :: sizes2, subsizes2, starts2
      integer rstatus(MPI_STATUS_SIZE)
      logical :: update_ok


      df_QLD_out(:,:,:,:,:) = 0D0;

      inode1=1  !JPL
      inode2=1  !JPL
      stride=1
      call get_mesh_node_range_for_threads_in_this_mpirank ( inode1, inode2,stride,  &
                                        nthreads, isize, i_beg, i_end, i_stride)

      isize=col_f_npsi   !JPL:: isize is defined by npsi?

       call init_col_f_module_dynamic( isize )   !! allocating memory for M_s (and M_ss)
      
      !$omp parallel default (none)                                          &
      !$omp shared(isize, sml_mype, col_f_npsi, col_f_dth) &
      !$omp shared(col_sp_cell_all)&
      !$omp shared(col_f_sp_num, gammac_spall)&
      !$omp shared(col_f_sp_s, col_f_sp_e, col_f_nvrm1,  f_in)&
      !$omp shared(tar_f_sp_s, tar_f_sp_e) &
      !$omp shared(sml_mass, sml_charge, eq_den, eq_t_ev) &
      !$omp shared(df_QLD_out, col_f_nthm1, f, g1, opt_QLD, g1_in)&
      !$omp shared(f0_M_tar_in, f0_M_tar_4d)             &
      !$omp shared(qld_B, qld_C, qld_F,rf_pwrscale) &
      !$omp shared(istart1,iend1,istart2,iend2) &  !JPL
      !$omp shared(col_iteration, i_1, i_2) &  !youknow
      !$omp private(ith, ipsi, update_ok, isp, imu, jsp) num_threads(col_f_o_nthreads)
      !$omp do
      do ipsi=istart1(sml_mype),iend1(sml_mype)    !  JPL: mpi for psi and th
        do ith=istart2(sml_mype)-1,iend2(sml_mype)-1

          do imu=0, col_f_nvrm1
            do isp = tar_f_sp_s, tar_f_sp_e ! for define the f0_M_tar

              if (isp.lt.col_f_sp_s) then ! Array size is [0:isp] but [0:isp-1] 
                f0_M_tar_4d(:,ith,imu,isp) = f0_M_tar_in(:,imu,ith,ipsi,isp) ! caution for order
              else
                f (:,ith,imu,isp) = f_in (:,imu,ith,ipsi,isp)  ! caution for order
                g1(:,ith,imu,isp) = g1_in(:,imu,ith,ipsi,isp)  ! caution for order
              endif
            enddo  !isp
          enddo   !imu

        !do ith=0,col_f_nthm1  !JPL
          call setup_col_sp(f,g1,f0_M_tar_4d,ith, sml_mass, eq_t_ev(ipsi,:), eq_den(ipsi,:),sml_charge,&
                            col_sp_cell_all(:,ipsi), gammac_spall(:,:,ipsi),&
                            ith,ipsi)

           !call MPI_Barrier(MPI_COMM_WORLD, ierr)  !JPL
        !!##### TIME ADVANCE
          select case (opt_QLD)
            case(1)
                ! NEED TO UPDATE
            case(2)
              ! if (col_f_sp_num.eq.1) then !!## Single Species Collisions
                  call QLD_fm_core(ith,ipsi,col_sp_cell_all(:,ipsi),update_ok,&
                  qld_B(ipsi,ith,:,:), qld_C(ipsi,ith,:,:), qld_F(ipsi,ith,:,:))! qld_B (1:npsi,0:nthm1,0:nmu,-nvp:nvp)
          end select

          do isp = col_f_sp_s, col_f_sp_e
            if (update_ok.and.(opt_QLD.ne.0)) then
              do imu=0, col_f_nvrm1
                df_QLD_out(:,imu,ith,ipsi,isp) = col_sp_cell_all(isp,ipsi)%pdf_np1(imu+1,:)
              enddo
            endif
          enddo ! isp

        enddo   ! ith
      end do !JPL: ipsi
 
!$omp end do
!$acc wait
!$omp end parallel

!$acc wait



!       JPL: MPI transfer
!       send to root

      tag    = 3001
      sizes    = [col_f_nvz,col_f_nvr,nth,npsi]
      starts   = [0,0,0,0]     ! let's say we're looking at region "0"  

      ntemp_p=1
      do ip1=1,nptot_mpi
        ntemp_p=ntemp_p*(sml_mype-ip1_mpi(ip1))
      enddo
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! write(*,*) 'works 0',sizes
      if ((sml_mype .ne. root).and.(ntemp_p.eq.0)) then
         subsizes = [col_f_nvz,col_f_nvr,iend2(sml_mype)-istart2(sml_mype)+1,iend1(sml_mype)-istart1(sml_mype)+1]
! size of sub-region
         call MPI_Type_create_subarray(4, sizes, subsizes, starts,MPI_ORDER_FORTRAN, MPI_double, newtype, ierr)
              ! write(*,*) 'works 00',sml_mype,subsizes
         call MPI_Type_commit(newtype, ierr)
!         ndata=nvr*nvz*(iend2(sml_mype)-istart2(sml_mype)+1)*(iend1(sml_mype)-istart1(sml_mype)+1)

         do  isp = col_f_sp_s, col_f_sp_e
                ! write(*,*) 'works 01',sml_mype,istart2(sml_mype)-1,istart1(sml_mype)
            call MPI_Send(df_QLD_out(-col_f_nvzm1/2,0,istart2(sml_mype)-1,istart1(sml_mype),isp), &
                 1, newtype, root, tag, MPI_COMM_WORLD, ierr)
                ! write(*,*) 'works 011',sml_mype
         enddo
      else if (sml_mype .eq. root) then
         do ip1=2,nptot_mpi
           ipe=ip1_mpi(ip1)
              ! write(*,*) 'works 02',sml_mype,ipe
          !  ndata=nvr*nvz*(iend2(ipe)-istart2(ipe)+1)*(iend1(ipe)-istart1(ipe)+1)
           subsizes = [col_f_nvz,col_f_nvr,iend2(ipe)-istart2(ipe)+1,iend1(ipe)-istart1(ipe)+1]
            ! write(*,*) 'works 03',sml_mype,subsizes
           call MPI_Type_create_subarray(4, sizes, subsizes,starts,MPI_ORDER_FORTRAN,MPI_double, newtype, ierr)
           call MPI_Type_commit(newtype, ierr)
            ! write(*,*) 'works 04',sml_mype,subsizes
           do  isp = col_f_sp_s, col_f_sp_e
                !  write(*,*) 'works 50',sml_mype,isp, ipe, istart2(ipe),df_QLD_out(-col_f_nvzm1/2,0,10,1,0)
                call MPI_Recv(df_QLD_out(-col_f_nvzm1/2,0,istart2(ipe)-1,istart1(ipe),isp), &
                     1, newtype, ipe, tag, MPI_COMM_WORLD, rstatus,ierr)

                !  write(*,*) 'works51',sml_mype,isp,df_QLD_out(-col_f_nvzm1/2,0,10,1,0)
           enddo
         enddo
      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

 do ipsi=1,isize    ! isize = outer thread num?????
        do ith=0,col_f_nthm1
             do isp = col_f_sp_s, col_f_sp_e
     
      ndata=col_f_nvz*col_f_nvr
      !write(*,*) 'works 91',sml_mype, ndata, ipsi,ith,isp

      call MPI_BCAST(df_QLD_out(-col_f_nvzm1/2,0,ith,ipsi,isp), ndata, MPI_double, ROOT, MPI_COMM_WORLD, ierr)
      !write(*,*) 'works 92',sml_mype, ndata, ipsi,ith,isp
enddo
enddo
enddo


       call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! write(*,*) 'works 1',sizes
         

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call finalize_col_f_module_dynamic( isize )

!rh unnecessary
!#ifdef _OPENMP
!      !  ---------------------
!      !  reset omp_num_threads
!      !  ---------------------
!      nthreads = omp_get_max_threads()
!      call omp_set_num_threads( nthreads )
!#endif

    end subroutine fm_QLD

    subroutine col_fm_core(ith,ipsi,col_spall,gammac,update_ok,sml_istep,&
                           qld_B34,qld_C34,qld_F34,&
                           col_f_mat,col_f_ksp,col_f_vecb,col_f_vecx,M_s, M_ab)
      use FP4D_globals, only: neg_frac, RF_n_harmonic ! youknow
      use FP4D_globals, only: opt_QLD

      use FP4D_globals, only : opt_time_solver, opt_time_species
      use FP4D_globals, only : eq_tau

      implicit none
      integer,intent(in) :: ith,ipsi
      integer :: col_f_mat, col_f_ksp, col_f_vecb, col_f_vecx
      !
      real (8), dimension(1:col_f_nvrm1,1:col_f_nvzm1) :: qld_B34, qld_C34, qld_F34 ! youknow

      ! type(col_f_sp_type), dimension(col_f_sp_s:col_f_sp_e) :: col_spall
      type(col_f_sp_type), dimension(tar_f_sp_s:tar_f_sp_e) :: col_spall
      ! real(kind=8), dimension(col_f_sp_s:col_f_sp_e, col_f_sp_s:col_f_sp_e) ::gammac
      real(kind=8), dimension(tar_f_sp_s:tar_f_sp_e, col_f_sp_s:col_f_sp_e) ::gammac
      logical :: update_ok
      real (kind=8) :: M_s( col_f_nvrm1, 5, col_f_nvrm1_nvzm1,col_f_sp_s:col_f_sp_e )
      ! real (kind=8), optional :: M_ab( col_f_nvrm1_nvzm1, 3, col_f_nvrm1_nvzm1, &
      !                                   col_f_sp_s:col_f_sp_e-1, col_f_sp_s:col_f_sp_e )
      
      ! youknow M_ab called M_ab(:,:,:,spj,spi)

      real (kind=8), optional :: M_ab( col_f_nvrm1_nvzm1, 3, col_f_nvrm1_nvzm1, &
                                        tar_f_sp_s:tar_f_sp_e-1, col_f_sp_s:col_f_sp_e )
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),2,(col_f_nvr-1)*(col_f_nvz-1)) :: M_ra_za
      real (kind=8), dimension(col_f_nvr, col_f_nvz) :: df  ! local
      real (kind=8), dimension(col_f_nvr, col_f_nvz) :: df_corr   ! local
      real (kind=8), dimension(col_f_ntotal_v, col_f_sp_s:col_f_sp_e) :: dist_col
      real (kind=8) :: negative_count
      real (kind=8),dimension(col_f_nvrm1,col_f_nvzm1) :: fhalf,dfdr,dfdz
      real (kind=8), dimension(5,col_f_nvrm1,col_f_nvzm1) :: ED
      integer :: iter_inter, spi, spj, spj_mod
      real (kind=8) :: col_dw_sum, col_dp_sum, col_w_sum, col_p_sum
      real (kind=8) :: col_dw, col_dp, col_dn_n, col_dn_n_max
      real (kind=8), dimension(LU_nnz, col_f_sp_s:col_f_sp_e) :: LU_values ! Super LU Variables
      real (kind=8), dimension(LU_nnz, tar_f_sp_s:tar_f_sp_e, col_f_sp_s:col_f_sp_e) :: LU_values_save ! Super LU Variables ! JPL
      real (kind=8), dimension(LU_nnz) :: LU_values_tmp ! Super LU Variables
      integer :: index_1,index_2,index_0, ir
      real (kind=8) :: smu_n
      logical :: dens_exit_ok, en_exit_ok, mom_exit_ok
      real (kind=8) :: min_p_thres
      integer :: index_I, index_J, index_c
      logical, external :: is_nan
      integer, intent(in) :: sml_istep

      real (kind=8) :: dt_QLD

#ifdef _OPENACC
      logical :: pM_ab
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num

      pM_ab = present(M_ab)
      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif

      update_ok = .false.   !!just for the case where sudden returning during calculation
      !! ## ----- INITIALIZATION
      col_w_sum = 0D0   !! sum of energy over species
      col_p_sum = 0D0   !! sum of momentum over species
      min_p_thres=1D99  !! momentum error criterion; just large value to be reinitialized with a small value

      ! do spi=col_f_sp_s, col_f_sp_e
      do spi=tar_f_sp_s, tar_f_sp_e ! youknow : IMPORTANT this must be define as target (all species)
          call col_fm_core_init(col_spall(spi))  !! evaluate col_f_sp_type%numeric_vth2, numeric_T, dens, ens, mom from " col_f_sp_type%pdf_n "
        col_w_sum = col_w_sum + col_spall(spi)%ens
        col_p_sum = col_p_sum + col_spall(spi)%mom
        min_p_thres=min(col_spall(spi)%mass*col_spall(spi)%dens  &
                   * sqrt(col_spall(spi)%numeric_vth2) , min_p_thres)
      enddo !! spi
      !rh Calculate finite difference parameters
      if (opt_Col.ne.0) then
      ! do spi=col_f_sp_s, col_f_sp_e
      do spi=tar_f_sp_s, tar_f_sp_e ! Actually, delta_init conflict with col_fm_f_df.
                                    ! delta_init(colliding,target) 
                                    ! -> define the colliding%delta for target 
                                    ! col_fm_f_df(target(target),colliding)
                                    ! -> call the target%delta for colliding
                                    ! But, the 'col_fm_f_df' maybe correct...?
                                    ! And.. it doesn't care expand the array size 
                                    ! from col&delta[tar] to tar%delta[tar] ==> let's go ~
          ! do spj=col_f_sp_s, col_f_sp_e
        do spj=tar_f_sp_s, tar_f_sp_e
          call col_fm_core_delta_init(col_spall(spi), col_spall(spj))  !! quasi-equilibrium parameter
        enddo
      enddo
      do spi=col_f_sp_s, col_f_sp_e
        ! do spj=col_f_sp_s, col_f_sp_e
        do spj=tar_f_sp_s, tar_f_sp_e

          select case (opt_col)
           case (0,1) 
              if(spi .eq. spj) then
                call col_fm_angle_avg_s(col_spall(spi), M_s(:,:,:,spi))  !! C(spa,spa)
              else
                if(spj .gt. spi) then
                  spj_mod = spj-1
                else
                  spj_mod = spj
                endif
                call col_fm_angle_avg_ab(col_spall(spi), col_spall(spj), M_ab(:,:,:,spj_mod,spi)) !! C(spi->spj)
              endif
            case (2)
              if(spi .eq. spj) then
                spj_mod= spj
                call col_fm_angle_avg_rab(col_spall(spi), col_spall(spj), M_ab(:,:,:,spj_mod,spi),M_ra_za,sml_istep)  !! C(spa,spa)
              else
                if(spj .gt. spi) then
                  spj_mod = spj
                else
                  spj_mod = spj+1
                endif
                call col_fm_angle_avg_rab(col_spall(spi), col_spall(spj), M_ab(:,:,:,spj_mod,spi),M_ra_za,sml_istep) !! C(spi->spj)
              endif
          end select
        enddo
      enddo
      !! "dist_iter" initializtion (making it zero) has been moved to initialization at make_initial_pdf_bimaxwell
      !! dist_iter <-implicitly updated distribution  (2D)

      !! ## ----- IMPLICIT TIME ADVANCE LOOP
      !!    NOTE THAT "dist_iter" is the PDF iterated and being updated in this LOOP
      !!    dist_iter -> col_f_sp_type%pdf_np1
      do iter_inter=1, vpic_inner_iter_max
        !! ##  b (=dist_col) = pdf_n; Ax=b
        !$omp  parallel do default(none)                                      &
        !$omp& shared(dist_col,col_spall,iter_inter,opt_QLD)                  &
        !$omp& shared(den_converge, ene_converge, col_iteration, RF_n_harmonic)  &
        !$omp& shared(col_f_sp_s, col_f_sp_e, col_f_nvr,col_f_nvz,opt_col)    &
        !$omp& private(index_1,index_2,index_0,spi)                           &
        !$omp& num_threads(col_f_i_nthreads)
        do spi=col_f_sp_s, col_f_sp_e ! youknow : only for colliding
          do index_2=1,col_f_nvz
            do index_1=1,col_f_nvr
              index_0 = index_1 + (index_2-1)*col_f_nvr
              dist_col(index_0,spi) = col_spall(spi)%pdf_n(index_1,index_2)   !! keep pdf_n for convergence calc. dist_col = pdf_n for matrix calc.
            enddo
          enddo
        enddo
        !! ### ----- LU_values (which should be zero at the beginning of this sum)
        LU_values = 0D0;
        LU_values_save = 0D0; ! JPL
        do spi=col_f_sp_s, col_f_sp_e
          !  do spj=col_f_sp_s, col_f_sp_e
           do spj=tar_f_sp_s, tar_f_sp_e
              call col_fm_f_df(col_spall(spj), spi, fhalf, dfdr, dfdz)  !! fhalf, dfdr, dfdz for col_spall(target sp) based on "pdf_np1"
                                    !! For the 1st step, pdf_np1 is already initialized to have same values with pdf_n at [setup_col_sp]
              ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
              select case (opt_col)
                case (0,1)
                  if(spi .eq. spj) then
                    call col_fm_E_and_D_s(col_spall(spi),fhalf,dfdr,dfdz,M_s(:,:,:,spi),ED) ! i-i
                  else
                    if(spj .gt. spi) then
                      spj_mod = spj-1
                    else
                      spj_mod = spj
                    endif
                    ! write(*,*) 'spi', spi, 'spj', spj
                    ! write(*,*) 'spi mass', col_spall(spi)%mass
                    ! write(*,*) 'spj mass', col_spall(spj)%mass

                    call col_fm_E_and_D_ab(col_spall(spi), col_spall(spj), fhalf, dfdr, dfdz, M_ab(:,:,:,spj_mod,spi), ED)
                  endif
                case (2)
                  if(spi .eq. spj) then
                    call col_fm_E_and_D_rab(col_spall(spi), col_spall(spj), fhalf, dfdr, dfdz, M_ab(:,:,:,spj_mod,spi),M_ra_za, ED) ! i-i
                    spj_mod= spj
                  else
                    if(spj .gt. spi) then
                      spj_mod = spj
                    else
                      spj_mod = spj+1
                    endif
                    call col_fm_E_and_D_rab(col_spall(spi), col_spall(spj), fhalf, dfdr, dfdz, M_ab(:,:,:,spj_mod,spi),M_ra_za, ED)
                  endif
              end select
              call col_f_LU_matrix(spj, col_spall(spi), ED, LU_values_tmp)  !! use original
              LU_values_save(:,spj,spi) = LU_values_tmp*gammac(spj,spi) !JPL: save for collisional power partition
              LU_values(:,spi) = LU_values(:,spi) + LU_values_tmp*gammac(spj,spi)
           enddo !spj
        enddo !spi

        ! write(*,*) 'LU_values_save'
        ! write(*,*) LU_values_save(10,0,0)
        ! write(*,*) LU_values_save(10,1,0)
        ! write(*,*) LU_values(10,0)


! youknow recommend to add the RF term as follows
        do spi=col_f_sp_s, col_f_sp_e
          ! in here, spi=spj
            ! select case(opt_TORIC)

            ! Set the dt_QLD in here, using col_f_dt
            if (opt_time_solver.eq.-1) then
              ! colliding particle have to have the same collision time. --> fix the 'b' when define the tau_ab
              dt_QLD = col_f_dt * eq_tau(ipsi,opt_time_species,spi) ! eq_tau(psi, target, colliding)
            else
              dt_QLD = col_f_dt
            endif
            
            if (opt_QLD.eq.2) then
              if (RF_n_harmonic(spi).eq.1) then ! col_f_sp_e is tmp
                call col_fm_E_and_D_s_QLD(col_spall(spi), ED, qld_B34, qld_C34, qld_F34) ! youknow
                call col_f_LU_matrix_QLD (spi, col_spall(spi), ED, LU_values_tmp)         !! use original
                LU_values(:,spi) = LU_values(:,spi) + LU_values_tmp * dt_QLD             ! col_f_dt not gammac
              endif
            endif
        enddo
! end youknow
  


        !! ### ----- Time advance & Exit condition check
        col_dw_sum = 0D0  !! sum of energy change over species (i.e. error in energy)
        col_dp_sum = 0D0  !! sum of momentum change over species (i.e. error in momentum)
        col_dn_n_max = 0D0 !! maximum fraction of density changes among species (i.e. maximum relative error in density in all species )
        dens_exit_ok = .true.
        do spi=col_f_sp_s, col_f_sp_e
          call cal_power_partition(col_spall(spi), LU_values_save(:,:,spi), dist_col(:,spi), iter_inter) !JPL: find collisional power partition



          call col_fm_picard_step(iter_inter, LU_values(:,spi), dist_col(:,spi), col_spall(spi)%pdf_np1)
          call col_fm_convergence_eval(col_spall(spi), col_dw, col_dp, col_dn_n)  !! moments evalulation from difference betw. pdf_n and pdf_np1
          col_dw_sum = col_dw_sum + col_dw
          ! write 'C(f)*dt'
          ! write(*,*) 'spi', spi, 'col_dw', col_dw
          ! do spj=tar_f_sp_s, tar_f_sp_e
          !   write(*,*) 'spj', spj, 'power ', col_spall(spi)%power(spj)
          ! enddo


          ! youknow :: convergence criterion for maxwellian species
          ! if ( col_f_sp_s .ne. tar_f_sp_s ) then
          !   do spj=tar_f_sp_s, col_f_sp_s-1
          !     col_dw_sum = col_dw_sum - col_spall(spi)%power(spj)
          !   enddo
          ! endif

          col_dp_sum = col_dp_sum + col_dp
          if(col_dn_n_max .lt. col_dn_n) col_dn_n_max = col_dn_n
          ! dens_exit_ok = dens_exit_ok .and. (col_dn_n .lt. 1D-8)
          dens_exit_ok = dens_exit_ok .and. (col_dn_n .lt. den_converge)
          
        enddo !spi

        en_exit_ok = dabs(col_dw_sum/col_w_sum) .le. ene_converge
        mom_exit_ok = dabs(col_dp_sum)/(max(dabs(col_p_sum),min_p_thres)) .le. (col_f_vp_max/col_f_nvz)   !! v_par = max_v_par of simulation domain / vth
        
        update_ok = dens_exit_ok .and. en_exit_ok .and. mom_exit_ok
          
        if (( update_ok .and. iter_inter .ne. 1 ).or.(iter_inter.ge.col_iteration)) then
            write(*,*) 'exit @ith:,',ith, ' @iter_inter:', iter_inter
            exit
        else
            if( iter_inter .gt. vpic_inner_iter_max-2) then
              write(*,'(a,i0,a)') 'Not converged @ith: ', ith, ' by collisions because'
              write(*,'(a,l,e22.12,a,e10.2)') 'exit condition: den', dens_exit_ok, col_dn_n, '/', den_converge
              write(*,'(a,l2,1x,e22.12,a,1x,e10.2,1x,2(e22.12,1x))') 'exit condition: ene',   &
                  (dabs(col_dw_sum/col_w_sum) .le. ene_converge ),   &
                    col_dw_sum/col_w_sum, '/', ene_converge, col_dw_sum, col_w_sum
              write(*,'(a,l,3(e22.12,1x),a,e22.12)') 'exit condition: mom',  &
                  (dabs(col_dp_sum)/(max(dabs(col_p_sum),min_p_thres)) .le. (col_f_vp_max/col_f_nvz)),  &
                    col_dp_sum, col_p_sum, min_p_thres, '/',(col_f_vp_max/col_f_nvz)
              write(*,'(a,i4)') 'exit condition: iteration number', iter_inter
            endif
        endif

      enddo !iter_inter   !! ### ----- implicit iteration loop ends here

      !! ### ----- DISTRIBUTION UPDATE (start)
      update_ok = .true.  !! Now check all fail cases. Other than that, it shold be updated
      if (iter_inter .ge. vpic_inner_iter_max) then
        write(*,'(a,2(i8,1x))') 'inner iteration went to maximum number. All pdfs are NOT updated:   &
                 mype, node = ', sml_mype, col_spall(col_f_sp_s)%nodeid
        update_ok= .false.
      else
        do spi=col_f_sp_s, col_f_sp_e
        ! Prevent negative values in new distribution function
        ! Heaviside function
        !$omp parallel do private(index_1,index_2) num_threads(col_f_i_nthreads)
          do index_2=1,size(df_corr,2)
            do index_1=1,size(df_corr,1)
              df_corr(index_1,index_2)=0.5D0*(sign(1D0,col_spall(spi)%pdf_np1(index_1,index_2))+1D0)
            enddo
          enddo
          negative_count=real(col_f_nvr*col_f_nvz,8) - sum(df_corr)

          if (negative_count .gt. neg_frac *real(col_f_nvr*col_f_nvz,8)) then
            write(*,'(a,2(i8,1x))') 'Too many negative values at :', spi, floor(negative_count)
            write(*,'(a)') 'All pdfs are NOT updated by collisions'
            update_ok = .false.
            exit
          else
          ! df = dist_iter - dist_n  ! df captures variation of distribution function
          !$omp parallel do private(index_1,index_2) num_threads(col_f_i_nthreads)
            do index_2=1,size(df,2)
              do index_1=1,size(df,1)
                df(index_1,index_2) = df_corr(index_1,index_2) * &
                                     (col_spall(spi)%pdf_np1(index_1,index_2) - col_spall(spi)%pdf_n(index_1,index_2))
              enddo
            enddo
          endif

          ! temporary Update of the distribution because next species may have many negative and pdf will not be updated for that case
          !$omp parallel do private(index_1,index_2) num_threads(col_f_i_nthreads)
          do index_2=1,size(df,2)
            do index_1=1,size(df,1)
              col_spall(spi)%pdf_np1(index_1,index_2) = df(index_1,index_2)  !!+  col_spall(spi)%pdf_n(index_1,index_2) !! this is from kernel (diff. logic)
            enddo
          enddo
        enddo !spi
      endif
      !! following legacy
      if (update_ok) then
        do spi=col_f_sp_s, col_f_sp_e
          do ir=1,col_f_nvr
            col_spall(spi)%pdf_np1(ir,:) = col_spall(spi)%pdf_np1(ir,:)/col_spall(spi)%conv_factor
          enddo
        enddo
      endif
     
!$acc wait(istream)
      endif
    end subroutine col_fm_core

    ! QLD_fm_core from col_fm_core
    subroutine QLD_fm_core(ith,ipsi,QLD_spall,update_ok,&
                           qld_B34,qld_C34,qld_F34)
      use FP4D_globals, only: neg_frac, RF_n_harmonic ! youknow
      use FP4D_globals, only: opt_QLD

      use FP4D_globals, only : opt_time_solver, opt_time_species
      use FP4D_globals, only : eq_tau

      implicit none
      integer,intent(in) :: ith,ipsi
      !
      real (8), dimension(1:col_f_nvrm1,1:col_f_nvzm1) :: qld_B34, qld_C34, qld_F34 ! youknow

      ! type(col_f_sp_type), dimension(col_f_sp_s:col_f_sp_e) :: QLD_spall
      type(col_f_sp_type), dimension(tar_f_sp_s:tar_f_sp_e) :: QLD_spall
      logical :: update_ok
      
      real (kind=8), dimension(col_f_nvr, col_f_nvz) :: df  ! local
      real (kind=8), dimension(col_f_nvr, col_f_nvz) :: df_corr   ! local
      real (kind=8), dimension(col_f_ntotal_v, col_f_sp_s:col_f_sp_e) :: dist_QLD
      real (kind=8) :: negative_count
      real (kind=8),dimension(col_f_nvrm1,col_f_nvzm1) :: fhalf,dfdr,dfdz
      real (kind=8), dimension(5,col_f_nvrm1,col_f_nvzm1) :: ED
      integer :: iter_inter, spi, spj, spj_mod
      real (kind=8) :: col_dw_sum, col_dp_sum, col_w_sum, col_p_sum
      real (kind=8) :: col_dw, col_dp, col_dn_n, col_dn_n_max
      real (kind=8), dimension(LU_nnz, col_f_sp_s:col_f_sp_e) :: LU_values ! Super LU Variables
      real (kind=8), dimension(LU_nnz) :: LU_values_tmp ! Super LU Variables
      integer :: index_1,index_2,index_0, ir
      real (kind=8) :: smu_n
      logical :: dens_exit_ok, en_exit_ok, mom_exit_ok
      real (kind=8) :: min_p_thres
      integer :: index_I, index_J, index_c
      logical, external :: is_nan

      real (kind=8) :: dt_QLD

#ifdef _OPENACC
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num

      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif

      update_ok = .false.   !!just for the case where sudden returning during calculation
      !! ## ----- INITIALIZATION
      col_w_sum = 0D0   !! sum of energy over species
      col_p_sum = 0D0   !! sum of momentum over species
      min_p_thres=1D99  !! momentum error criterion; just large value to be reinitialized with a small value

      ! do spi=col_f_sp_s, col_f_sp_e
      do spi=tar_f_sp_s, tar_f_sp_e ! youknow : IMPORTANT this must be define as target (all species)
          call col_fm_core_init(QLD_spall(spi))  !! evaluate col_f_sp_type%numeric_vth2, numeric_T, dens, ens, mom from " col_f_sp_type%pdf_n "
        col_w_sum = col_w_sum + QLD_spall(spi)%ens
        col_p_sum = col_p_sum + QLD_spall(spi)%mom
        min_p_thres=min(QLD_spall(spi)%mass*QLD_spall(spi)%dens  &
                   * sqrt(QLD_spall(spi)%numeric_vth2) , min_p_thres)
      enddo !! spi

      !rh Calculate finite difference parameters

      do spi=tar_f_sp_s, tar_f_sp_e 
        do spj=tar_f_sp_s, tar_f_sp_e
          call col_fm_core_delta_init(QLD_spall(spi), QLD_spall(spj))  !! quasi-equilibrium parameter
        enddo
      enddo

      !! "dist_iter" initializtion (making it zero) has been moved to initialization at make_initial_pdf_bimaxwell
      !! dist_iter <-implicitly updated distribution  (2D)

      !! ## ----- IMPLICIT TIME ADVANCE LOOP
      !!    NOTE THAT "dist_iter" is the PDF iterated and being updated in this LOOP
      !!    dist_iter -> col_f_sp_type%pdf_np1
      do iter_inter=1, vpic_inner_iter_max
        !! ##  b (=dist_QLD) = pdf_n; Ax=b
        !$omp  parallel do default(none)                                      &
        !$omp& shared(dist_QLD,QLD_spall,iter_inter,opt_QLD)                  &
        !$omp& shared(den_converge, ene_converge, col_iteration, RF_n_harmonic)  &
        !$omp& shared(col_f_sp_s, col_f_sp_e, col_f_nvr,col_f_nvz,opt_col)    &
        !$omp& private(index_1,index_2,index_0,spi)                           &
        !$omp& num_threads(col_f_i_nthreads)
        do spi=col_f_sp_s, col_f_sp_e ! youknow : only for colliding
          do index_2=1,col_f_nvz
            do index_1=1,col_f_nvr
              index_0 = index_1 + (index_2-1)*col_f_nvr
              dist_QLD(index_0,spi) = QLD_spall(spi)%pdf_n(index_1,index_2)   !! keep pdf_n for convergence calc. dist_QLD = pdf_n for matrix calc.
            enddo
          enddo
        enddo



        !! ### ----- LU_values (which should be zero at the beginning of this sum)
        LU_values = 0D0;

        ! youknow recommend to add the RF term as follows
        do spi=col_f_sp_s, col_f_sp_e
            ! in here, spi=spj
            ! select case(opt_TORIC)

            ! Set the dt_QLD
            if (opt_time_solver.eq.-1) then
              ! colliding particle have to have the same collision time. --> fix the 'b' when define the tau_ab
              dt_QLD = QLD_f_dt * eq_tau(ipsi,opt_time_species,spi) ! eq_tau(psi, target, colliding)
            else
              dt_QLD = QLD_f_dt
            endif

            if (opt_QLD.eq.2) then
              if (RF_n_harmonic(spi).eq.1) then ! col_f_sp_e is tmp
                call col_fm_E_and_D_s_QLD(QLD_spall(spi), ED, qld_B34, qld_C34, qld_F34) ! youknow
                call col_f_LU_matrix_QLD (spi, QLD_spall(spi), ED, LU_values_tmp)         !! use original
                LU_values(:,spi) = LU_values(:,spi) + LU_values_tmp * dt_QLD             ! dt_QLD
              endif
            endif
        enddo
  
        !! ### ----- Time advance & Exit condition check
        col_dw_sum = 0D0  !! sum of energy change over species (i.e. error in energy)
        col_dp_sum = 0D0  !! sum of momentum change over species (i.e. error in momentum)
        col_dn_n_max = 0D0 !! maximum fraction of density changes among species (i.e. maximum relative error in density in all species )
        dens_exit_ok = .true.
        do spi=col_f_sp_s, col_f_sp_e

          call col_fm_picard_step(iter_inter, LU_values(:,spi), dist_QLD(:,spi), QLD_spall(spi)%pdf_np1)
          call col_fm_convergence_eval(QLD_spall(spi), col_dw, col_dp, col_dn_n)  !! moments evalulation from difference betw. pdf_n and pdf_np1
          col_dw_sum = col_dw_sum + col_dw
          ! write 'C(f)*dt'
          ! write(*,*) 'spi', spi, 'col_dw', col_dw
          ! do spj=tar_f_sp_s, tar_f_sp_e
          !   write(*,*) 'spj', spj, 'power ', QLD_spall(spi)%power(spj)
          ! enddo


          ! youknow :: convergence criterion for maxwellian species
          ! if ( col_f_sp_s .ne. tar_f_sp_s ) then
          !   do spj=tar_f_sp_s, col_f_sp_s-1
          !     col_dw_sum = col_dw_sum - QLD_spall(spi)%power(spj)
          !   enddo
          ! endif

          col_dp_sum = col_dp_sum + col_dp
          if(col_dn_n_max .lt. col_dn_n) col_dn_n_max = col_dn_n
          ! dens_exit_ok = dens_exit_ok .and. (col_dn_n .lt. 1D-8)
          dens_exit_ok = dens_exit_ok .and. (col_dn_n .lt. den_converge)
          
        enddo !spi

        en_exit_ok = dabs(col_dw_sum/col_w_sum) .le. ene_converge
        mom_exit_ok = dabs(col_dp_sum)/(max(dabs(col_p_sum),min_p_thres)) .le. (col_f_vp_max/col_f_nvz)   !! v_par = max_v_par of simulation domain / vth
        
        update_ok = dens_exit_ok .and. en_exit_ok .and. mom_exit_ok
          
        if (( update_ok .and. iter_inter .ne. 1 ).or.(iter_inter.ge.col_iteration)) then
            write(*,*) 'exit @ith:,',ith, ' @iter_inter:', iter_inter
            exit
        else
            if( iter_inter .gt. vpic_inner_iter_max-2) then
              write(*,'(a,i0,a)') 'Not converged @ith: ', ith, ' by collisions because'
              write(*,'(a,l,e22.12,a,e10.2)') 'exit condition: den', dens_exit_ok, col_dn_n, '/', den_converge
              write(*,'(a,l2,1x,e22.12,a,1x,e10.2,1x,2(e22.12,1x))') 'exit condition: ene',   &
                  (dabs(col_dw_sum/col_w_sum) .le. ene_converge ),   &
                    col_dw_sum/col_w_sum, '/', ene_converge, col_dw_sum, col_w_sum
              write(*,'(a,l,3(e22.12,1x),a,e22.12)') 'exit condition: mom',  &
                  (dabs(col_dp_sum)/(max(dabs(col_p_sum),min_p_thres)) .le. (col_f_vp_max/col_f_nvz)),  &
                    col_dp_sum, col_p_sum, min_p_thres, '/',(col_f_vp_max/col_f_nvz)
              write(*,'(a,i4)') 'exit condition: iteration number', iter_inter
            endif
        endif

      enddo !iter_inter   !! ### ----- implicit iteration loop ends here

      !! ### ----- DISTRIBUTION UPDATE (start)
      update_ok = .true.  !! Now check all fail cases. Other than that, it shold be updated
      if (iter_inter .ge. vpic_inner_iter_max) then
        write(*,'(a,2(i8,1x))') 'inner iteration went to maximum number. All pdfs are NOT updated:   &
                 mype, node = ', sml_mype, QLD_spall(col_f_sp_s)%nodeid
        update_ok= .false.
      else
        do spi=col_f_sp_s, col_f_sp_e
        ! Prevent negative values in new distribution function
        ! Heaviside function
        !$omp parallel do private(index_1,index_2) num_threads(col_f_i_nthreads)
          do index_2=1,size(df_corr,2)
            do index_1=1,size(df_corr,1)
              df_corr(index_1,index_2)=0.5D0*(sign(1D0,QLD_spall(spi)%pdf_np1(index_1,index_2))+1D0)
            enddo
          enddo
          negative_count=real(col_f_nvr*col_f_nvz,8) - sum(df_corr)

          if (negative_count .gt. neg_frac *real(col_f_nvr*col_f_nvz,8)) then
            write(*,'(a,2(i8,1x))') 'Too many negative values at :', spi, floor(negative_count)
            write(*,'(a)') 'All pdfs are NOT updated by collisions'
            update_ok = .false.
            exit
          else
          ! df = dist_iter - dist_n  ! df captures variation of distribution function
          !$omp parallel do private(index_1,index_2) num_threads(col_f_i_nthreads)
            do index_2=1,size(df,2)
              do index_1=1,size(df,1)
                df(index_1,index_2) = df_corr(index_1,index_2) * &
                                     (QLD_spall(spi)%pdf_np1(index_1,index_2) - QLD_spall(spi)%pdf_n(index_1,index_2))
              enddo
            enddo
          endif

          ! temporary Update of the distribution because next species may have many negative and pdf will not be updated for that case
          !$omp parallel do private(index_1,index_2) num_threads(col_f_i_nthreads)
          do index_2=1,size(df,2)
            do index_1=1,size(df,1)
              QLD_spall(spi)%pdf_np1(index_1,index_2) = df(index_1,index_2)  !!+  QLD_spall(spi)%pdf_n(index_1,index_2) !! this is from kernel (diff. logic)
            enddo
          enddo
        enddo !spi
      endif
      !! following legacy
      if (update_ok) then
        do spi=col_f_sp_s, col_f_sp_e
          do ir=1,col_f_nvr
            QLD_spall(spi)%pdf_np1(ir,:) = QLD_spall(spi)%pdf_np1(ir,:)/QLD_spall(spi)%conv_factor
          enddo
        enddo
      endif
     
!$acc wait(istream)
    end subroutine QLD_fm_core

    !> This routine computes the current moments of the distribution
    !! function: density, parallel and perpendicular temperature,
    !! mean parallel flow and entropy.
    !!
    subroutine col_fm_core_init(spa)
      use FP4D_globals, only:sml_istep
      use FP4D_globals, only : sml_e_charge, sml_ev2j
#ifdef _OPENMP
      use omp_lib, only: omp_get_thread_num
#endif
      implicit none
      type(col_f_sp_type) :: spa   !! species 'a'
      !
      real (8) :: densi, rdensi, numeric_Ti, eni, momi,entropy_i, eni_perp, eni_par, curi
      real (8) :: reni
      integer :: index_i, index_j, ithread,z,r

      ! Get Eq. Temperature
      densi = 0D0
      momi = 0D0
      eni_perp = 0D0
      eni_par = 0D0
      entropy_i = 0D0
      !$omp parallel do private(index_i) reduction(+:densi,eni_perp,eni_par, &
      !$omp momi,curi,entropy_i) shared(col_f_vp_max) num_threads(col_f_i_nthreads)
      do index_i=1, col_f_nvz
        z=col_f_dvp*(real(index_i)-real(col_f_nvz,8)/2D0)  
        densi = densi + sum(spa%pdf_n(:,index_i)*spa%vol)
        momi = momi+sum(spa%pdf_n(:,index_i)*spa%mesh_z(index_i)*spa%vol)
        eni_perp = eni_perp+sum(spa%pdf_n(:,index_i)*spa%mesh_r**2*spa%vol)
        eni_par  = eni_par+sum(spa%pdf_n(:,index_i)*((spa%vpar_beg+spa%mesh_dz*(index_i-1))**2)*spa%vol)
        entropy_i = entropy_i - sum(spa%pdf_n(:,index_i)*log(spa%pdf_n(:,index_i))*spa%vol)
        if (opt_Col.eq.2) then
          do index_j=1,col_f_nvr
            r=col_f_dsmu*real(index_j-1,8)	
	    reni = eni + 2*(sqrt(1+(((spa%mesh_z(index_i))**2+(spa%mesh_r(index_j))**2))/(299792458D0**2))-1)*(299792458D0**2)*(spa%pdf_n(index_j,index_i))*spa%vol(index_j)
            if ((z*z+r*r)**0.5D0.LT.col_f_vp_max) then
              rdensi = rdensi + sum(spa%pdf_n(index_j,index_i)*spa%vol)
            else
            end if
          end do
        endif
      end do
      curi = momi*spa%Ze
      momi=spa%mass*momi
      eni_perp = eni_perp*spa%mass
      eni_par = eni_par*spa%mass
      eni = eni_par+eni_perp
      reni = reni*spa%mass
      numeric_Ti = eni/(3D0*densi*sml_e_charge)
      !rh We cannot take the flow into account to get the equilibrium temperature in the
      !rh multi-species case. In equilibrium, ions and electrons would flow together. However
      !rh the equilibrium flow is not known. Part of the ion and electron flows will dissipate to
      !rh heat. Therefore, we use the mean kinetic energy to estimate the equilibrium temperature
      spa%numeric_t      = numeric_Ti
      spa%numeric_vth2   = numeric_Ti*sml_e_charge/spa%mass !  sp%mass --> ptl_mass
      !!spa%numeric_dens   = densi

      select case (opt_Col)
        case (0,1)
          spa%dens           = densi
          spa%ens            = eni
          spa%mom            = momi
        case (2)
          spa%dens           = rdensi
          spa%ens            = reni
          spa%mom            = momi
      end select
      return

    end subroutine col_fm_core_init


    !> This routine evaluates the gyroaveraged Landau
    !! interaction tensor for collisions between different
    !! particle species --> Sec. 2.2 in Hager et al.
    !!
    subroutine col_fm_angle_avg_ab(c1, c2, M_ab)
      use FP4D_globals, only : sml_pi     !! XGC1_3
      use elliptics_mod
      implicit none
      type(col_f_sp_type) :: c1, c2
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),3,(col_f_nvr-1)*(col_f_nvz-1)) :: M_ab
      integer :: index_I, index_ip, index_J, index_jp
      integer :: index_rz, index_ac
      real (8) :: r, z, a, c, dz, dz2, r2, a2, lambda, k, kp1_sqrt, km1
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1)) :: EK, EE, k_eff
      integer, dimension((col_f_nvr-1)*(col_f_nvz-1)) :: vpic_ierr
      real (8) :: EE_k, EK_k
      real (8) :: I1, I2, I3, temp_cons
      integer :: mesh_Nrm1, mesh_Nzm1
      real(8)::c1_mesh_r_half(lbound(c1%mesh_r_half,1):ubound(c1%mesh_r_half,1))
      real(8)::c1_mesh_z_half(lbound(c1%mesh_z_half,1):ubound(c1%mesh_z_half,1))
      real(8)::c2_mesh_r_half(lbound(c2%mesh_r_half,1):ubound(c2%mesh_r_half,1))
      real(8)::c2_mesh_z_half(lbound(c2%mesh_z_half,1):ubound(c2%mesh_z_half,1))
      logical, parameter :: use_ellip_subroutine = .true.
      !!============================== use_ellip_subroutine=.false.
      integer, parameter ::n_order = 10
      integer :: i,order
      real(kind=8), parameter :: e_tol = 1.0d-8
      integer, parameter :: vl = 64
      real(kind=8), dimension(vl) :: x_n,y_n,Mx_n,My_n,Mz_n
      real(kind=8) :: x_np1,Mx_np1,My_np1,rd
      integer :: iibeg,iiend,iisize,ii
      !!===============================
      logical, external :: is_nan
#ifdef _OPENACC
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num

      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
!!!!!$acc enter data pcreate(M_ab) ASYNC(istream)
#endif

      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1

      c1_mesh_r_half = c1%mesh_r_half
      c1_mesh_z_half = c1%mesh_z_half
      c2_mesh_r_half = c2%mesh_r_half
      c2_mesh_z_half = c2%mesh_z_half
      ! new routine
      ! (r,z) : colliding paritcle "c1" - capital index, (a,c) : target particle "c2" - small index
#ifdef _OPENACC
!$acc  kernels  present(M_ab)  ASYNC(istream) &
!$acc& pcopyin(mesh_Nzm1,mesh_Nrm1) &
!$acc& pcopyin(c1_mesh_r_half,c1_mesh_z_half) &
!$acc& pcopyin(c2_mesh_r_half,c2_mesh_z_half)

!$acc  loop independent collapse(2)  gang &
!$acc& private( index_I, z, index_J, index_rz, r, r2) &
!$acc& private(index_ip, c, dz, dz2, index_jp, a, a2, index_ac, lambda) &
!$acc& private(k, k_eff, kp1_sqrt, km1, EK, EE, vpic_ierr, I1, I2, I3)  &
!$acc& private(temp_cons) &
!$acc& private(EE_k, EK_k)    &
!$acc& private(x_n,y_n,Mx_n,My_n,Mz_n,i,x_np1,Mx_np1,My_np1,rd) &
!$acc& private(ii,iibeg,iiend,iisize)
#else
!$OMP PARALLEL DO &
!$OMP& default(none) &
!$OMP& shared(mesh_Nzm1,mesh_Nrm1,c1,c2,M_ab, &
!$OMP&        c1_mesh_z_half, c1_mesh_r_half, c2_mesh_z_half, c2_mesh_r_half) &
!$OMP& PRIVATE( index_I, z, index_J, index_rz, r, r2, index_ip, c, dz, dz2, index_jp, a, a2, index_ac, lambda, &
!$OMP&           k, k_eff, kp1_sqrt, km1, EK, EE, vpic_ierr, I1, I2, I3, temp_cons, &
!$OMP&           EE_k, EK_k, x_n, y_n, Mx_n, My_n, Mz_n, i, x_np1, Mx_np1, My_np1, rd, &
!$OMP&           ii,iibeg,iiend,iisize) &
!$OMP&           num_threads(col_f_i_nthreads)
#endif
      do index_I=1,mesh_Nzm1
          do index_J=1,mesh_Nrm1
              z = c1_mesh_z_half(index_I)
              index_rz = index_J +mesh_Nrm1*(index_I-1)
              r = c1_mesh_r_half(index_J)
              r2 = r*r

#ifdef _OPENACC
!$acc loop independent collapse(2) worker vector
#endif
              do index_ip=1, mesh_Nzm1
                  do index_jp=1, mesh_Nrm1
                      c = c2_mesh_z_half(index_ip)
                      dz = z-c
                      dz2 = dz*dz
                      a = c2_mesh_r_half(index_jp)
                      a2 = a*a
                      index_ac = index_jp+mesh_Nrm1*(index_ip-1)

                      !lambda
                      lambda = r2+a2+dz2     !symmetric (r<->a, z<->c)
                                             ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z) for SAME grid
                                             !      BUT!!, (r,z,a,c) = (a,c,r,z) only due to limitation of our domain for each species
      
                      k = 2D0*r*a/lambda     !symmetric
                      k_eff(index_ac) = 2D0*k/(1D0+k)
                  enddo
              enddo
#ifdef USE_VECTOR
              call ellip_agm_v(k_eff, EK, EE, vpic_ierr(1), size(k_eff))
#else
              if (use_ellip_subroutine) then
#ifdef _OPENACC
!$acc loop independent collapse(2) worker vector
#endif
                  do index_ip=1, mesh_Nzm1
                      do index_jp=1, mesh_Nrm1
                          index_ac = index_jp+mesh_Nrm1*(index_ip-1)
                          call ellip_agm( k_eff(index_ac),  &
                              EK(index_ac),EE(index_ac),vpic_ierr(index_ac))
                      enddo
                  enddo
              else
#ifdef _OPENACC
!$acc loop independent  worker
#endif
                  do iibeg=1,mesh_Nzm1*mesh_Nrm1,vl
                      iiend = min(mesh_Nzm1*mesh_Nrm1, iibeg+vl-1)
                      iisize = iiend - iibeg + 1
#ifdef _OPENACC
!$acc loop independent  vector
#endif
                      do i=1,iisize
                          ii = iibeg + (i-1)
                          Mx_n(i) = 1D0
                          My_n(i) = 1D0 - k_eff(ii)    !beta^2
                          Mz_n(i) = 0D0
                          x_n(i) = 1D0
                          y_n(i) = sqrt(My_n(i))  !beta
                      enddo

                      do order=1, n_order
#ifdef _OPENACC
!$acc loop independent vector
#endif
                          do i=1,iisize
                              x_np1 = (x_n(i)+y_n(i))*0.5D0
                              y_n(i) = sqrt(x_n(i)*y_n(i))
                              x_n(i) = x_np1 !update results
                              !magm
                              Mx_np1 = (Mx_n(i)+My_n(i))*0.5D0
                              rd = sqrt( (Mx_n(i)-Mz_n(i))*(My_n(i)-Mz_n(i)) )
                              My_np1 = Mz_n(i)+rd
                              Mz_n(i) = Mz_n(i)-rd
                              Mx_n(i) = Mx_np1
                              My_n(i) = My_np1
                          enddo
                      enddo
#ifdef _OPENACC
!$acc loop independent  vector
#endif
                      do i=1,iisize
                          ii = iibeg + (i-1)
                          EK(ii) = (0.5D0*sml_pi)/x_n(i)
                          EE(ii) = Mx_n(i)*EK(ii)
                      enddo
                  enddo
              endif
#endif

#ifdef _OPENACC
!$acc loop independent collapse(2) worker vector
#endif
              do index_ip=1, mesh_Nzm1
                  do index_jp=1, mesh_Nrm1
                      c = c2_mesh_z_half(index_ip)
                      dz = z-c
                      dz2 = dz*dz
                      a = c2_mesh_r_half(index_jp)
                      a2 = a*a
                      index_ac = index_jp+mesh_Nrm1*(index_ip-1)

                      if(index_ac .eq. index_rz) then
                          M_ab(index_ac,:,index_rz) = 0D0  !!! THIS LINE SHOULD BE INVESTIGATED. CONDITIONAL NEEDED?
                          cycle
                      endif

                      !lambda
                      lambda = r2+a2+dz2     !symmetric (r<->a, z<->c)
                                             ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z) for SAME grid
                                             !      BUT!!, (r,z,a,c) = (a,c,r,z) only due to limitation of our domain for each species
  
                      k = 2D0*r*a/lambda     !symmetric
                      kp1_sqrt = sqrt(1D0+k)
                      km1 = k-1D0

                      EE_k = EE(index_ac)
                      EK_k = EK(index_ac)
                      !Calulation of M coeff. (all symmetric)
                      I1 = -4D0*((1D0+k)*EE_k-EK_k)/(k*k*kp1_sqrt)
                      I2 = -2D0*EE_k/(km1*kp1_sqrt)
                      I3 = -2D0*(EE_k+km1*EK_k)/(km1*k*kp1_sqrt )
                      !I4 = I2-I1 !For exact Numerical Conservation, and It should be mathematically
                      temp_cons = 4D0*sml_pi/(lambda*sqrt(lambda))

                      M_ab(index_ac,1,index_rz) = temp_cons*(I1*a2+I2*dz2)
                      M_ab(index_ac,2,index_rz) = temp_cons*dz*(I3*a-I2*r)
                      M_ab(index_ac,3,index_rz) = temp_cons*(I2*(r2+a2) - I3*2D0*r*a)
                  enddo !index_jp
              enddo !index_ip
          enddo !index_J
      enddo ! index_I

!$acc end kernels
!$acc wait(istream)

    end subroutine col_fm_angle_avg_ab

    subroutine col_fm_angle_avg_rab(c1, c2, M_ab,M_ra_za,sml_istep)
      use FP4D_globals, only : sml_pi
      use elliptics_mod
      use spline_interpolation
      use H5_rdwt
      use HDF5

      implicit none
      type(col_f_sp_type) :: c1, c2
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),3,(col_f_nvr-1)*(col_f_nvz-1)) :: M_ab
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),5,(col_f_nvr-1)*(col_f_nvz-1)) :: M_Ucomp
      ! M_Ucomp saves all 5 U components for first time step, after first time step substitute M_Ucomp value into M_ab and M_ra_za 
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),2,(col_f_nvr-1)*(col_f_nvz-1)) :: M_ra_za
      ! NUperpperpprime and NUparlperpprime are saved in M_ra_za
      real (8) :: NUperpperp,NUperpparl,NUparlparl,NUperpperpprime,NUparlperpprime	
      integer :: index_I, index_ip, index_J, index_jp
      integer :: index_rz, index_ac
      real (8) :: r, z, a, c, dz, dz2, r2, a2, lambda, k, kp1_sqrt, km1
      real (8), dimension((col_f_nvr-1)*(col_f_nvz-1)) :: EK, EE, k_eff
      integer, dimension((col_f_nvr-1)*(col_f_nvz-1)) :: vpic_ierr
      real (8) :: EE_k, EK_k
      real (8) :: I1, I2, I3, temp_cons, result
      integer :: mesh_Nrm1, mesh_Nzm1
      real(8)::c1_mesh_r_half(lbound(c1%mesh_r_half,1):ubound(c1%mesh_r_half,1))
      real(8)::c1_mesh_z_half(lbound(c1%mesh_z_half,1):ubound(c1%mesh_z_half,1))
      real(8)::c2_mesh_r_half(lbound(c2%mesh_r_half,1):ubound(c2%mesh_r_half,1))
      real(8)::c2_mesh_z_half(lbound(c2%mesh_z_half,1):ubound(c2%mesh_z_half,1))
      real(8), Dimension(4)	:: utable			! utable = ur=1, urp=2, uz=3, uzp=4
      logical, parameter :: use_ellip_subroutine = .true.
      !!============================== use_ellip_subroutine=.false.
      integer, parameter ::n_order = 10
      integer :: i,order
      real(kind=8), parameter :: e_tol = 1.0d-8
      integer, parameter :: vl = 64
      real(kind=8), dimension(vl) :: x_n,y_n,Mx_n,My_n,Mz_n
      real(kind=8) :: x_np1,Mx_np1,My_np1,rd
      real(kind=8) :: k1value, tvalue
      integer :: iibeg,iiend,iisize,ii
      Double precision, parameter		:: cvel = 299792458D0		! light velocity
      integer :: sml_istep	
      !!===============================
      logical, external :: is_nan


#ifdef _OPENACC
	
      call cpu_time(tinit)
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num

      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
!!!!!$acc enter data pcreate(M_ab) ASYNC(istream)
#endif

      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1   
      c1_mesh_r_half = c1%mesh_r_half
      c1_mesh_z_half = c1%mesh_z_half
      c2_mesh_r_half = c2%mesh_r_half
      c2_mesh_z_half = c2%mesh_z_half
      ! new routine
      ! (r,z) : colliding paritcle "c1" - capital index, (a,c) : target particle "c2" - small index
#ifdef _OPENACC
!$acc  kernels  present(M_ab)  ASYNC(istream) &
!$acc& pcopyin(mesh_Nzm1,mesh_Nrm1) &
!$acc& pcopyin(c1_mesh_r_half,c1_mesh_z_half) &
!$acc& pcopyin(c2_mesh_r_half,c2_mesh_z_half)

!$acc  loop independent collapse(2)  gang &
!$acc& private( index_I, z, index_J, index_rz, r, r2) &
!$acc& private(index_ip, c, dz, dz2, index_jp, a, a2, index_ac, lambda) &
!$acc& private(k, k_eff, kp1_sqrt, km1, EK, EE, vpic_ierr, I1, I2, I3)  &
!$acc& private(temp_cons) &
!$acc& private(EE_k, EK_k)    &
!$acc& private(x_n,y_n,Mx_n,My_n,Mz_n,i,x_np1,Mx_np1,My_np1,rd) &
!$acc& private(ii,iibeg,iiend,iisize)
#else
!$OMP PARALLEL DO &
!$OMP& default(none) &
!$OMP& shared(mesh_Nzm1,mesh_Nrm1,c1,c2,M_ab,M_ra_za,M_Ucomp, &
!$OMP&		
!$OMP&          c1_mesh_z_half, c1_mesh_r_half, c2_mesh_z_half, c2_mesh_r_half,&
!$OMP&		 sml_istep) &
!$OMP& PRIVATE( index_I, z, index_J, index_rz, r, r2, index_ip, c, dz, dz2, index_jp, a, a2, index_ac, lambda, &
!$OMP&           k, k_eff, kp1_sqrt, km1, EK, EE, vpic_ierr, I1, I2, I3, temp_cons, &
!$OMP&		NUperpperp,NUperpparl,NUparlparl,NUperpperpprime,NUparlperpprime, &
!$OMP&           EE_k, EK_k, x_n, y_n, Mx_n, My_n, Mz_n, i, x_np1, Mx_np1, My_np1, rd, result, k1value, tvalue, &
!$OMP&           ii,iibeg,iiend,iisize) &
!$OMP&           num_threads(col_f_i_nthreads)
#endif

      do index_I=1,mesh_Nzm1
          do index_J=1,mesh_Nrm1
              z = c1_mesh_z_half(index_I)
	       	
              index_rz = index_J +mesh_Nrm1*(index_I-1)
              r = c1_mesh_r_half(index_J)
              r2 = r*r	      

#ifdef _OPENACC
!$acc loop independent collapse(2) worker vector
#endif
              do index_ip=1, mesh_Nzm1
                  do index_jp=1, mesh_Nrm1
                      c = c2_mesh_z_half(index_ip)
		      
                      dz = z-c
                      dz2 = dz*dz
                      a = c2_mesh_r_half(index_jp)
                      a2 = a*a
                      index_ac = index_jp+mesh_Nrm1*(index_ip-1)

                      lambda = r2+a2+dz2     !symmetric (r<->a, z<->c)
                                             ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z) for SAME grid
                                             !      BUT!!, (r,z,a,c) = (a,c,r,z) only due to limitation of our domain for each species
      
                      k = 2D0*r*a/lambda     !symmetric
                      k_eff(index_ac) = 2D0*k/(1D0+k)
                  enddo
              enddo
#ifdef USE_VECTOR
              call ellip_agm_v(k_eff, EK, EE, vpic_ierr(1), size(k_eff))
#else
              if (use_ellip_subroutine) then
#ifdef _OPENACC
!$acc loop independent collapse(2) worker vector
#endif
                  do index_ip=1, mesh_Nzm1
                      do index_jp=1, mesh_Nrm1
                          index_ac = index_jp+mesh_Nrm1*(index_ip-1)
                          call ellip_agm( k_eff(index_ac),  &
                              EK(index_ac),EE(index_ac),vpic_ierr(index_ac))
                      enddo
                  enddo
              else
#ifdef _OPENACC
!$acc loop independent  worker
#endif
                  do iibeg=1,mesh_Nzm1*mesh_Nrm1,vl
                      iiend = min(mesh_Nzm1*mesh_Nrm1, iibeg+vl-1)
                      iisize = iiend - iibeg + 1
#ifdef _OPENACC
!$acc loop independent  vector
#endif
                      do i=1,iisize
                          ii = iibeg + (i-1)
                          Mx_n(i) = 1D0
                          My_n(i) = 1D0 - k_eff(ii)    !beta^2
                          Mz_n(i) = 0D0
                          x_n(i) = 1D0
                          y_n(i) = sqrt(My_n(i))  !beta
                      enddo

                      do order=1, n_order
#ifdef _OPENACC
!$acc loop independent vector
#endif
                          do i=1,iisize
                              x_np1 = (x_n(i)+y_n(i))*0.5D0
                              y_n(i) = sqrt(x_n(i)*y_n(i))
                              x_n(i) = x_np1 !update results
                              !magm
                              Mx_np1 = (Mx_n(i)+My_n(i))*0.5D0
                              rd = sqrt( (Mx_n(i)-Mz_n(i))*(My_n(i)-Mz_n(i)) )
                              My_np1 = Mz_n(i)+rd
                              Mz_n(i) = Mz_n(i)-rd
                              Mx_n(i) = Mx_np1
                              My_n(i) = My_np1
                          enddo
                      enddo
#ifdef _OPENACC
!$acc loop independent  vector
#endif
                      do i=1,iisize
                          ii = iibeg + (i-1)
                          EK(ii) = (0.5D0*sml_pi)/x_n(i)
                          EE(ii) = Mx_n(i)*EK(ii)
                      enddo
                  enddo
              endif
#endif

#ifdef _OPENACC
!$acc loop independent collapse(2) worker vector
#endif
              do index_ip=1, mesh_Nzm1
                  do index_jp=1, mesh_Nrm1
                      c = c2_mesh_z_half(index_ip)
		      
                      dz = z-c
                      dz2 = dz*dz
                      a = c2_mesh_r_half(index_jp)
		     
                      a2 = a*a
                      index_ac = index_jp+mesh_Nrm1*(index_ip-1)

                      if(index_ac .eq. index_rz) then
                          M_ab(index_ac,:,index_rz) = 0D0  !!! THIS LINE SHOULD BE INVESTIGATED. CONDITIONAL NEEDED?
                          cycle
                      endif

                      !lambda
                      lambda = r2+a2+dz2     !symmetric (r<->a, z<->c)
                                             ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z) for SAME grid
                                             !      BUT!!, (r,z,a,c) = (a,c,r,z) only due to limitation of our domain for each species
  
                      k = 2D0*r*a/lambda     !symmetric
                      kp1_sqrt = sqrt(1D0+k)
                      km1 = k-1D0
		      
                      EE_k = EE(index_ac)
                      EK_k = EK(index_ac)
                      !Calulation of M coeff. (all symmetric)
                      I1 = -4D0*((1D0+k)*EE_k-EK_k)/(k*k*kp1_sqrt)
                      I2 = -2D0*EE_k/(km1*kp1_sqrt)
                      I3 = -2D0*(EE_k+km1*EK_k)/(km1*k*kp1_sqrt )
                      !I4 = I2-I1 !For exact Numerical Conservation, and It should be mathematically
                      temp_cons = 4D0*sml_pi/(lambda*sqrt(lambda))

		      if (sml_istep .eq. 1) then	
		      !======================================================== Numerical Integration================================================!
		      
		      call gauss_chebyshev_quad(r,a,z,c,32, NUperpperp,NUperpparl,NUparlparl,NUperpperpprime,NUparlperpprime,Tvalue)

		      !======================================================== Numerical Integration================================================!
		      
		      M_ab(index_ac,1,index_rz) = NUperpperp
		      M_ab(index_ac,2,index_rz) = NUperpparl
                      M_ab(index_ac,3,index_rz) = NUparlparl
		      M_ra_za(index_ac,1,index_rz) = NUperpperpprime
		      M_ra_za(index_ac,2,index_rz) = NUparlperpprime

		      else
		 	! if sml_istep != 1, substitute saved value to M_ab, M_ra_za
		      	M_ab(index_ac,1,index_rz)= M_Ucomp(index_ac,1,index_rz)
		      	M_ab(index_ac,2,index_rz)= M_Ucomp(index_ac,2,index_rz)
		      	M_ab(index_ac,3,index_rz)= M_Ucomp(index_ac,3,index_rz)
		      	M_ra_za(index_ac,1,index_rz)= M_Ucomp(index_ac,4,index_rz)
		      	M_ra_za(index_ac,2,index_rz)= M_Ucomp(index_ac,5,index_rz)
		      	cycle	

		      endif
		      M_Ucomp(index_ac,1,index_rz)=M_ab(index_ac,1,index_rz)
		      M_Ucomp(index_ac,2,index_rz)=M_ab(index_ac,2,index_rz)
		      M_Ucomp(index_ac,3,index_rz)=M_ab(index_ac,3,index_rz)
		      M_Ucomp(index_ac,4,index_rz)=M_ra_za(index_ac,1,index_rz)
		      M_Ucomp(index_ac,5,index_rz)=M_ra_za(index_ac,2,index_rz)



                  enddo !index_jp
              enddo !index_ip
          enddo !index_J
      enddo ! index_I

!$acc end kernels
!$acc wait(istream)


    end subroutine col_fm_angle_avg_rab


    !> This routine evaluates the gyroaveraged Landau
    !! interaction tensor for collisions between particesl of
    !! the same species --> Sec. 2.2 in Hager et al.
    !! This is cheaper because some symmetries can be exploited.
    !!
    !!M_zr => Ms(2,:,:)
    !!M_rc => Ms(2,:,:)
    !!M_zc => Ms(3,:,:)
    !! cs1 : colliding, cs2 : target
    !! f_half, dfdr, dfdz => cs2
    subroutine col_fm_E_and_D_ab(cs1, cs2, f_half, dfdr, dfdz, Ms, EDs)
      implicit none
      type(col_f_sp_type) :: cs1, cs2
      real (8), dimension(col_f_nvr-1,col_f_nvz-1) :: f_half, dfdr, dfdz
      !!real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),3,(col_f_nvr-1)*(col_f_nvz-1)) :: Ms
      real (8), dimension(col_f_nvrm1_nvzm1,3,col_f_nvrm1_nvzm1) :: Ms
      real (8), dimension(5,col_f_nvr-1, col_f_nvz-1) :: EDs
      integer :: index_I, index_J, index_ip, index_jp, mesh_Nrm1, mesh_Nzm1, index_2dp, index_2D
      real (8) :: tmp_vol, tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v, mass1_inv, mass2_inv
      real (8) :: tmpr1, tmpr2, tmpr3, tmpr4, tmpr5  !Dmitry
      real (8) :: r, a, z, c, dz, M_ra, M_za ! ES: array reduction
      real (8) :: cs1_mesh_r_half(lbound(cs1%mesh_r_half,1):ubound(cs1%mesh_r_half,1))
      real (8) :: cs1_mesh_z_half(lbound(cs1%mesh_z_half,1):ubound(cs1%mesh_z_half,1))
      real (8) :: cs2_mesh_r_half(lbound(cs2%mesh_r_half,1):ubound(cs2%mesh_r_half,1))
      real (8) :: cs2_mesh_z_half(lbound(cs2%mesh_z_half,1):ubound(cs2%mesh_z_half,1))
      real (8) :: cs2_local_center_volume(lbound(cs2%local_center_volume,1):ubound(cs2%local_center_volume,1))
#ifdef _OPENACC
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num
      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif

      cs1_mesh_r_half = cs1%mesh_r_half
      cs1_mesh_z_half = cs1%mesh_z_half
      cs2_mesh_r_half = cs2%mesh_r_half
      cs2_mesh_z_half = cs2%mesh_z_half
      cs2_local_center_volume = cs2%local_center_volume

      mass1_inv = 1D0/cs1%mass
      mass2_inv = 1D0/cs2%mass

      ! write(*,*) 'mass1_inv', mass1_inv, 'mass2_inv', mass2_inv

      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1

#ifdef _OPENACC
!$acc enter data pcreate(EDs) ASYNC(istream)
!$acc kernels present(Ms,EDs) ASYNC(istream)
!$acc loop independent  collapse(2) gang
#else
!$OMP  PARALLEL DO &
!$OMP& default(none) &
!$OMP& shared(mesh_Nzm1,mesh_Nrm1,f_half,dfdr,dfdz,Ms) &
!$OMP& shared(cs1,cs2,EDs,mass1_inv,mass2_inv, &
!$OMP&        cs1_mesh_r_half, cs1_mesh_z_half, cs2_mesh_r_half, cs2_mesh_z_half,cs2_local_center_volume ) &
!$OMP& PRIVATE(index_I,index_J, index_2D, index_ip, index_jp, index_2dp, &
!$OMP&          tmp_vol, tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v,  &
!$OMP&          r,a,z,c,dz,M_ra,M_za, &
!$OMP&          tmpr1, tmpr2, tmpr3, tmpr4, tmpr5 ) &
!$OMP& num_threads(col_f_i_nthreads)
#endif
      do index_I=1, mesh_Nzm1
          do index_J=1, mesh_Nrm1
              z = cs1_mesh_z_half(index_I)
              r = cs1_mesh_r_half(index_J)
              index_2D = index_J+mesh_Nrm1*(index_I-1)
              tmpr1=0D0
              tmpr2=0D0
              tmpr3=0D0
              tmpr4=0D0
              tmpr5=0D0
!$acc         loop independent collapse(2) vector
              do index_ip = 1, mesh_Nzm1
                  do index_jp = 1, mesh_Nrm1
                      c = cs2_mesh_z_half(index_ip)
                      dz = z-c
                      a = cs2_mesh_r_half(index_jp)
                      index_2dp = index_jp+mesh_Nrm1*(index_ip-1)
                      tmp_vol = cs2_local_center_volume(index_jp)
                      tmp_f_half_v = f_half(index_jp, index_ip) * tmp_vol
                      tmp_dfdr_v = dfdr(index_jp, index_ip) * tmp_vol
                      tmp_dfdz_v = dfdz(index_jp, index_ip) * tmp_vol
                      M_ra = (r*Ms(index_2dp,1,index_2D) + dz*Ms(index_2dp,2,index_2D))/a  !ES: array reduction
                      M_za = (r*Ms(index_2dp,2,index_2D) + dz*Ms(index_2dp,3,index_2D))/a  !ES: array reduction

                      ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
                      !  Ms = 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za
                      tmpr1 = tmpr1 + Ms(index_2dp,1,index_2D)*tmp_f_half_v
                      tmpr2 = tmpr2 + Ms(index_2dp,2,index_2D)*tmp_f_half_v
                      tmpr3 = tmpr3 + Ms(index_2dp,3,index_2D)*tmp_f_half_v
                      tmpr4 = tmpr4 + M_ra*tmp_dfdr_v + Ms(index_2dp,2,index_2D)*tmp_dfdz_v
                      tmpr5 = tmpr5 + M_za*tmp_dfdr_v + Ms(index_2dp,3,index_2D)*tmp_dfdz_v  !Multi-species
                  enddo !index_jp
              enddo !index_ip

              ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
              !mass correction
              EDs(1,index_J,index_I) = tmpr1*mass1_inv
              EDs(2,index_J,index_I) = tmpr2*mass1_inv
              EDs(3,index_J,index_I) = tmpr3*mass1_inv
              EDs(4,index_J,index_I) = tmpr4*mass2_inv
              EDs(5,index_J,index_I) = tmpr5*mass2_inv
          enddo
      enddo

!$acc end kernels
!$acc wait(istream)
!$acc exit data copyout(EDs) ASYNC(istream)
!$acc wait(istream)
    end subroutine col_fm_E_and_D_ab

    subroutine col_fm_E_and_D_rab(cs1, cs2, f_half, dfdr, dfdz, Ms, M_ra_za, EDs)
    ! What is Different? : M_ra_za is added
      implicit none
      type(col_f_sp_type) :: cs1, cs2
      real (8), dimension(col_f_nvr-1,col_f_nvz-1) :: f_half, dfdr, dfdz
      !!real (8), dimension((col_f_nvr-1)*(col_f_nvz-1),3,(col_f_nvr-1)*(col_f_nvz-1)) :: Ms
      real (8), dimension(col_f_nvrm1_nvzm1,3,col_f_nvrm1_nvzm1) :: Ms
      real (8), dimension(col_f_nvrm1_nvzm1,2,col_f_nvrm1_nvzm1) :: M_ra_za
    ! What is Different? : M_ra_za is added
      real (8), dimension(5,col_f_nvr-1, col_f_nvz-1) :: EDs
      integer :: index_I, index_J, index_ip, index_jp, mesh_Nrm1, mesh_Nzm1, index_2dp, index_2D
      real (8) :: tmp_vol, tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v, mass1_inv, mass2_inv
      real (8) :: tmpr1, tmpr2, tmpr3, tmpr4, tmpr5 !Dmitry
      real (8) :: r, a, z, c, dz, M_ra, M_za ! ES: array reduction
      real (8) :: cs1_mesh_r_half(lbound(cs1%mesh_r_half,1):ubound(cs1%mesh_r_half,1))
      real (8) :: cs1_mesh_z_half(lbound(cs1%mesh_z_half,1):ubound(cs1%mesh_z_half,1))
      real (8) :: cs2_mesh_r_half(lbound(cs2%mesh_r_half,1):ubound(cs2%mesh_r_half,1))
      real (8) :: cs2_mesh_z_half(lbound(cs2%mesh_z_half,1):ubound(cs2%mesh_z_half,1))
      real (8) :: cs2_local_center_volume(lbound(cs2%local_center_volume,1):ubound(cs2%local_center_volume,1))
#ifdef _OPENACC
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num
      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif
     
      cs1_mesh_r_half = cs1%mesh_r_half
      cs1_mesh_z_half = cs1%mesh_z_half
      cs2_mesh_r_half = cs2%mesh_r_half
      cs2_mesh_z_half = cs2%mesh_z_half
      cs2_local_center_volume = cs2%local_center_volume

      mass1_inv = 1D0/cs1%mass
      mass2_inv = 1D0/cs2%mass

      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1

#ifdef _OPENACC
!$acc enter data pcreate(EDs) ASYNC(istream)
!$acc kernels present(Ms,EDs) ASYNC(istream)
!$acc loop independent  collapse(2) gang
#else
!$OMP  PARALLEL DO &
!$OMP& default(none) &
!$OMP& shared(mesh_Nzm1,mesh_Nrm1,f_half,dfdr,dfdz,Ms) &
!$OMP& shared(cs1,cs2,EDs,mass1_inv,mass2_inv, &
!$OMP&        cs1_mesh_r_half, cs1_mesh_z_half, cs2_mesh_r_half, cs2_mesh_z_half,cs2_local_center_volume,M_ra_za) &
!$OMP& PRIVATE(index_I,index_J, index_2D, index_ip, index_jp, index_2dp, &
!$OMP&          tmp_vol, tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v,  &
!$OMP&          r,a,z,c,dz,M_ra,M_za, &
!$OMP&          tmpr1, tmpr2, tmpr3, tmpr4, tmpr5 ) &
!$OMP& num_threads(col_f_i_nthreads)
#endif
      do index_I=1, mesh_Nzm1
          do index_J=1, mesh_Nrm1
              z = cs1_mesh_z_half(index_I)
              r = cs1_mesh_r_half(index_J)
              index_2D = index_J+mesh_Nrm1*(index_I-1)
              tmpr1=0D0
              tmpr2=0D0
              tmpr3=0D0
              tmpr4=0D0
              tmpr5=0D0
!$acc         loop independent collapse(2) vector
              do index_ip = 1, mesh_Nzm1
                  do index_jp = 1, mesh_Nrm1
                      c = cs2_mesh_z_half(index_ip)
                      dz = z-c
                      a = cs2_mesh_r_half(index_jp)
                      index_2dp = index_jp+mesh_Nrm1*(index_ip-1)
                      tmp_vol = cs2_local_center_volume(index_jp)
                      tmp_f_half_v = f_half(index_jp, index_ip) * tmp_vol
                      tmp_dfdr_v = dfdr(index_jp, index_ip) * tmp_vol
                      tmp_dfdz_v = dfdz(index_jp, index_ip) * tmp_vol
		      ! Substitute M_ra_za data into M_ra, M_za
		      M_ra  = M_ra_za(index_2dp,1,index_2D)
		      M_za  = M_ra_za(index_2dp,2,index_2D)

                      ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
                      !  Ms = 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za
                      tmpr1 = tmpr1 + Ms(index_2dp,1,index_2D)*tmp_f_half_v
                      tmpr2 = tmpr2 + Ms(index_2dp,2,index_2D)*tmp_f_half_v
                      tmpr3 = tmpr3 + Ms(index_2dp,3,index_2D)*tmp_f_half_v
                      tmpr4 = tmpr4 + M_ra*tmp_dfdr_v + Ms(index_2dp,2,index_2D)*tmp_dfdz_v
                      tmpr5 = tmpr5 + M_za*tmp_dfdr_v + Ms(index_2dp,3,index_2D)*tmp_dfdz_v  !Multi-species

                  enddo !index_jp
              enddo !index_ip
	      
              ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
              !mass correction
              EDs(1,index_J,index_I) = tmpr1*mass1_inv
              EDs(2,index_J,index_I) = tmpr2*mass1_inv
              EDs(3,index_J,index_I) = tmpr3*mass1_inv
              EDs(4,index_J,index_I) = tmpr4*mass2_inv
              EDs(5,index_J,index_I) = tmpr5*mass2_inv
          enddo
      enddo
!$acc end kernels
!$acc wait(istream)
!$acc exit data copyout(EDs) ASYNC(istream)
!$acc wait(istream)
    end subroutine col_fm_E_and_D_rab

   !! same with col_f_angle_avg_s. But need to check perofrmance issue for simplified code script
    subroutine col_fm_angle_avg_s(cs, Ms)
      use FP4D_globals, only : sml_pi     !! XGC1_3
      use elliptics_mod
      implicit none
      type(col_f_sp_type) :: cs   !! DIFF: type
      !!real (8), dimension(col_f_nvr-1,5,(col_f_nvr-1)*(col_f_nvz-1)) :: Ms
      real (8), dimension(col_f_nvrm1,5,col_f_nvrm1_nvzm1) :: Ms
      integer :: index_dz, index_J, index_jp
      integer :: index_rz
      real (8) :: r, a, dz, dz2, r2, a2, lambda, k, k_eff, kp1_sqrt, km1, EK, EE, mesh_dz
      integer :: vpic_ierr
      real (8) :: I1, I2, I3, temp_cons
      integer :: mesh_Nrm1, mesh_Nzm1
      !!real (8) :: M_rr_v, M_rz_v, M_zz_v, M_ra_v, M_za_v
      real (8) :: tmp_vol(col_f_nvr-1), tmp_volr
      integer :: vpic_ierr0
      real*8::cs_mesh_r_half(lbound(cs%mesh_r_half,1):ubound(cs%mesh_r_half,1))
#ifdef _OPENACC
      integer :: ithread,istream
      !integer, external :: omp_get_thread_num

      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif
      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1

      vpic_ierr = 0
      vpic_ierr0 = 0
      cs_mesh_r_half = cs%mesh_r_half
      mesh_dz = cs%mesh_dz
      tmp_vol = cs%local_center_volume  ! It is volume for prime coordinate. cs2. For same species, both cs1 and cs 2 has same
#ifdef _OPENACC
!$acc kernels                                              &
!$acc& pcopyin(mesh_Nrm1,mesh_Nzm1,cs_mesh_r_half,tmp_vol,mesh_dz)       &
!$acc& present(Ms)

!$acc loop independent collapse(2) gang  &
!$acc& private(EE,EK, index_J, index_dz, index_jp, r, r2, index_rz, a, a2, tmp_volr, dz, dz2) &
!$acc& private(lambda, k, k_eff, kp1_sqrt, km1, I1, I2, I3, temp_cons,vpic_ierr0)
#else
!$omp parallel do default(none) collapse(3)                             &
!$omp& shared(mesh_Nrm1,mesh_Nzm1,cs_mesh_r_half,tmp_vol,mesh_dz)       &
!$omp& shared(Ms)                                                       &
!$omp& private(index_J,r,r2,index_dz,index_rz,index_jp)                 &
!$omp& private(a,a2,tmp_volr,dz,dz2,lambda,k,k_eff,kp1_sqrt,km1)        &
!$omp& private(vpic_ierr0) reduction(max:vpic_ierr)                     &
!$omp& private(EK,EE,I1,I2,I3,temp_cons)                                &
!$omp& num_threads(col_f_i_nthreads)
#endif
     do index_J=1,mesh_Nrm1
        do index_dz=0,mesh_Nzm1-1
#ifdef _OPENACC
!$acc loop independent vector
!!$acc& private(index_jp, r, r2, index_rz, a, a2, tmp_volr, dz, dz2) &
!!$acc& private(lambda, k, k_eff, kp1_sqrt, km1, I1, I2, I3, temp_cons)&
!!!!$acc& private(M_rr_v, M_rz_v, M_zz_v, M_ra_v, M_za_v)
#endif
           do index_jp=1, mesh_Nrm1
              r=cs_mesh_r_half(index_J)
              r2 = r*r

              index_rz = index_J +mesh_Nrm1*index_dz   ! (=index_2D in MATLAB)
              !!index_ac = index_jp ! (=index_2d in MATLAB)
              if(index_jp .eq. index_rz) then
                 Ms(index_jp,:,index_rz) = 0D0
                 cycle
              endif
              a=cs_mesh_r_half(index_jp)
              a2 = a*a
              tmp_volr = tmp_vol(index_jp)

              dz = mesh_dz*index_dz
              dz2 = dz*dz


              lambda = r2+a2+dz2     !symmetric (r<->a, z<->c)
                                   ! i.e. (r,z,a,c) = (a,z,r,c) = (r,c,a,z) = (a,c,r,z)
              k = 2D0*r*a/lambda     !symmetric
              k_eff = 2D0*k/(1D0+k)
              kp1_sqrt = sqrt(1D0+k)
              km1 = k-1D0

              !k_eff=min(k_eff,0.95D0)
               
              call ellip_agm(k_eff,EK,EE,vpic_ierr0)
              !vpic_ierr = max(vpic_ierr,vpic_ierr0)   ! ES : This line is to be considered again
              vpic_ierr = 0   ! ES : This line is to be considered again

              !Calulation of M coeff. (all symmetric)
              I1 = -4D0*((1D0+k)*EE-EK)/(k*k*kp1_sqrt)
              I2 = -2D0*EE/(km1*kp1_sqrt)
              I3 = -2D0*(EE+km1*EK)/(km1*k*kp1_sqrt )
              !I4 = I2-I1 !For exact Numerical Conservation, and It shoudl be mathematically
              temp_cons = 4D0*sml_pi/(lambda*sqrt(lambda)) * tmp_volr  ! tmp_volr is from E_and_D  - PATCH02-2

              !NOTE THAT CONSIDERING SYMMETRY, we already have
              ! I1, I2, I3 at (r,z,a,c), (a,z,r,c), (r,c,a,z), and (a,c,r,z) ->(jp,ip,J,I)
              ! Using this values, we calculate symmetric and asymmetric compoenents of matrix
              ! index_rz = index_J +mesh_Nrm1*index_I
              ! index_rc = index_J +mesh_Nrm1*index_ip
              ! index_ac = index_jp+mesh_Nrm1*index_ip
              ! index_az = index_jp+mesh_Nrm1*index_I
              Ms(index_jp,1,index_rz) = temp_cons*(I1*a2+I2*dz2)
              Ms(index_jp,2,index_rz) = temp_cons*dz*(I3*a-I2*r)
              Ms(index_jp,3,index_rz) = temp_cons*(I2*(r2+a2) - I3*2D0*r*a)
              Ms(index_jp,4,index_rz) = temp_cons*(I1*r*a + I3*dz2)
              Ms(index_jp,5,index_rz) = temp_cons*dz*(I2*a-I3*r)
              !(a,c,r,z)=(r,z,a,c)
              ! 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra
              !!Ms(index_ac,1,index_rz) = M_rr_v
              !!Ms(index_ac,2,index_rz) = M_rz_v
              !!Ms(index_ac,3,index_rz) = M_zz_v
              !!Ms(index_ac,4,index_rz) = M_ra_v
              !!Ms(index_ac,5,index_rz) = M_za_v
           enddo  !index_jp
        enddo ! index_dz
     enddo !index_J
#ifdef _OPENACC
!$acc end kernels
!$acc wait(istream)
#endif
    
    end subroutine col_fm_angle_avg_s


    !> Calculation of the distribution function and its finite difference
    !! derivative between grid points for conservation of the thermal equilibrium
    !! as discussed in section 2.5 in Hager et al.
    subroutine col_fm_f_df(cs, op_mode, f_half, dfdr, dfdz)
      implicit none
      type(col_f_sp_type), intent(in) :: cs
      integer, intent(in) :: op_mode
      real(kind=8), dimension(col_f_nvr-1, col_f_nvz-1) :: f_half, dfdr, dfdz
      real(kind=8) :: tmpr1, tmpr2, tmpr3, tmpr4
      integer :: index_I, index_J
      ! For finite difference parameter
      real (kind=8) :: delr, cdelr, delz, cdelz

      !$omp parallel do private(index_I,index_J,tmpr1,tmpr2,tmpr3,tmpr4,delr,cdelr,delz,cdelz) num_threads(col_f_i_nthreads)
      do index_I=1,col_f_nvz-1
          do index_J=1, col_f_nvr-1
              delr = cs%delta_r(index_J,op_mode)
              cdelr = 1D0-delr
              delz = cs%delta_z(index_I,op_mode)
              cdelz = 1D0-delz

              tmpr1 = cs%pdf_np1(index_J, index_I)
              tmpr2 = cs%pdf_np1(index_J+1, index_I)
              tmpr3 = cs%pdf_np1(index_J,index_I+1)
              tmpr4 = cs%pdf_np1(index_J+1,index_I+1)

              ! With delta_(r,z)=0.5
              !f_half(index_J, index_I) = (tmpr1 + tmpr3 + tmpr2 + tmpr4)*0.25D0
              !dfdr(index_J, index_I) = ((tmpr2 - tmpr1)+(tmpr4 - tmpr3))*0.5D0
              !dfdz(index_J, index_I) = ((tmpr3 - tmpr1)+(tmpr4 - tmpr2))*0.5D0

              ! With finite difference factor
              f_half(index_J, index_I) = tmpr1 * delr*delz &
                                        + tmpr3 * delr*cdelz &
                                        + tmpr2 * cdelr*delz &
                                        + tmpr4 * cdelr*cdelz
              dfdr(index_J, index_I) = (tmpr2 - tmpr1)*delz &
                                      +(tmpr4 - tmpr3)*cdelz
              dfdz(index_J, index_I) = (tmpr3 - tmpr1)*delr &
                                      +(tmpr4 - tmpr2)*cdelr
          enddo
      enddo

      ! dfdr = dfdr/mesh_dr
      !$omp parallel do private(index_I,index_J) num_threads(col_f_i_nthreads)
      do index_I=1,size(dfdr,2)
        do index_J=1,size(dfdr,1)
           dfdr(index_J,index_I) = dfdr(index_J,index_I)/cs%mesh_dr
        enddo
      enddo

      ! dfdz = dfdz/mesh_dz
      !$omp parallel do private(index_I,index_J) num_threads(col_f_i_nthreads)
      do index_I=1,size(dfdz,2)
        do index_J=1,size(dfdz,1)
          dfdz(index_J,index_I) = dfdz(index_J,index_I)/cs%mesh_dz
        enddo
      enddo

      return

    end subroutine col_fm_f_df


   subroutine col_fm_E_and_D_s(cs,f_half,dfdr,dfdz,Ms,EDs)
      !cs1 : colliding, cs2 : target
      !f_half, dfdr, dfdz => cs2
      implicit none
      type(col_f_sp_type) :: cs
      real (8) :: mass,Ze
      real (8), dimension(col_f_nvr-1,col_f_nvz-1) :: f_half, dfdr, dfdz
      real (8), dimension(col_f_nvr-1,5, col_f_nvrm1_nvzm1) :: Ms
      real (8), dimension(5,col_f_nvrm1, col_f_nvzm1) :: EDs
      integer :: index_I, index_J, index_ip, index_jp, mesh_Nrm1, mesh_Nzm1, index_2D, sav_i, index_dz
      real (8) :: tmp_f_half_v, tmp_dfdr_v, tmp_dfdz_v
      real (8) :: sum1,sum2,sum3,sum4,sum5,inv_mass
      real (8) :: cs_mesh_r_half(lbound(cs%mesh_r_half,1):ubound(cs%mesh_r_half,1))
      real (8) :: cs_mesh_z_half(lbound(cs%mesh_z_half,1):ubound(cs%mesh_z_half,1))
      real (8) :: csign
      integer :: lb1,ub1, lb2,ub2, lb3,ub3,spid
#ifdef _OPENACC
      integer :: istream, ithread
      !integer, external :: omp_get_thread_num

      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif

      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1
      mass      = cs%mass
      Ze        = cs%Ze
      spid      = cs%spid
      cs_mesh_r_half = cs%mesh_r_half
      cs_mesh_z_half = cs%mesh_z_half
#ifdef _OPENACC
!$acc  enter data pcreate(EDs) ASYNC(istream)

!$acc  parallel ASYNC(istream)                           &
!$acc& present(Ms,EDs)                                   &
!$acc& pcopyin(mesh_Nrm1,mesh_Nzm1,f_half,dfdr,dfdz)

!$acc loop independent  collapse(2) gang  &
!$acc& private(csign,index_jp,index_ip,index_2D,index_dz,sav_i,tmp_f_half_v,tmp_dfdr_v,tmp_dfdz_v)
#else
!$omp parallel do default(none)                                      &
!$omp& collapse(2)                                                   &
!$omp& shared(mesh_Nrm1,mesh_Nzm1,f_half,dfdr,dfdz,Ze,spid)       &
!$omp& shared(Ms,EDs,mass,cs_mesh_r_half,cs_mesh_z_half)&
!$omp& private(index_J,index_dz,index_2D,index_jp)                   &
!$omp& private(sav_i,index_I,index_ip)                               &
!$omp& private(tmp_f_half_v,tmp_dfdr_v,tmp_dfdz_v)                   &
!$omp& private(sum1,sum2,sum3,sum4,sum5)                     &
!$omp& private(csign)                                                &
!$omp& num_threads(col_f_i_nthreads)
#endif
      do index_J=1, mesh_Nrm1
       do index_I=1,mesh_Nzm1
          !
          ! (1) index_I = 1 + index_dz + sav_i
          ! (2) 0 <= sav_i <= mesh_Nzm1-1
          ! (3) 0 <= index_dz <= mesh_Nzm1-1
          ! from (1)    sav_i = (index_I-1) - index_dz
          ! from (3) (index_I-1) - (mesh_Nzm1-1) <=  sav_i <= (index_I-1) - 0
          ! from (2) max(0,(index_I-1)-(mesh_Nzm1-1)) <= sav_i <= min((index_I-1),mesh_Nzm1-1)
          ! index_dz = index_I-1-sav_i

          sum1 = 0
          sum2 = 0
          sum3 = 0
          sum4 = 0
          sum5 = 0

!!$acc loop independent worker reduction(+:sum1,sum2,sum3,sum4,sum5)
          do sav_i=0,(mesh_Nzm1-1)
            if (sav_i <= (index_I-1)) then
              ! index_I  = 1+index_dz+sav_i
              index_dz = index_I-1-sav_i
              csign = 1
            else
              !index_I = 1-index_dz+sav_i
              index_dz = 1 + sav_i - index_I
              csign = -1
            endif
            index_ip = 1+sav_i
            index_2D = index_J+mesh_Nrm1*index_dz  ! This is index for the Ms


            ! (original) EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
            ! (reduced)  EDs = 1: Drr, 2: Drz, 3:Dzz, 4: Er, 5:Ez   ! Drz = Dzr
            !  Ms = 1: M_rr, 2:M_rz, 3:M_zz, 4:M_ra, 5:M_za
!!$acc loop independent vector
            do index_jp=1, mesh_Nrm1
              !!index_2dp = index_jp  ! This is index for the Ms
              tmp_f_half_v = f_half(index_jp, index_ip)
              tmp_dfdr_v   = dfdr(index_jp, index_ip)
              tmp_dfdz_v   = dfdz(index_jp, index_ip)
              sum1 = sum1 + Ms(index_jp,1,index_2D)*tmp_f_half_v
              sum2 = sum2 + csign*Ms(index_jp,2,index_2D)*tmp_f_half_v
              sum3 = sum3 + Ms(index_jp,3,index_2D)*tmp_f_half_v
              sum4 = sum4 + Ms(index_jp,4,index_2D)*tmp_dfdr_v + &
                            csign*Ms(index_jp,2,index_2D)*tmp_dfdz_v
              sum5 = sum5 + Ms(index_jp,3,index_2D)*tmp_dfdz_v + &
                            csign*Ms(index_jp,5,index_2D)*tmp_dfdr_v
            enddo
          enddo
          EDs(1,index_J,index_I) =  sum1
          EDs(2,index_J,index_I) =  sum2
          EDs(3,index_J,index_I) =  sum3
          EDs(4,index_J,index_I) =  sum4
          EDs(5,index_J,index_I) =  sum5
        enddo !index_I
      enddo !index_J
!$acc end parallel
!$acc wait(istream)
      ! (original) EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
      ! (reduced)  EDs = 1: Drr, 2: Drz, 3:Dzz, 4: Er, 5:Ez
      !mass correction
      !EDs = EDs/mass
      inv_mass = 1D0/mass
      lb1 = lbound(EDs,1)
      lb2 = lbound(EDs,2)
      lb3 = lbound(EDs,3)
      ub1 = ubound(EDs,1)
      ub2 = ubound(EDs,2)
      ub3 = ubound(EDs,3)

#ifdef _OPENACC
!$acc  kernels  ASYNC(istream)                          &
!$acc& pcopyin(inv_mass)                                &
!$acc& pcopyin(lb1,ub1,lb2,ub2,lb3,ub3)                 &
!$acc& present(EDs)
#else
!$omp  parallel do default(none)                        &
!$omp& private(index_i,index_j,index_ip)                &
!$omp& shared(lb1,ub1,lb2,ub2,lb3,ub3)                  &
!$omp& shared(EDs,inv_mass)                             &
!$omp& num_threads(col_f_i_nthreads)
#endif
      do index_j=lb3,ub3
         do index_i=lb2,ub2
            do index_ip=lb1,ub1
               EDs(index_ip,index_i,index_j) = EDs(index_ip,index_i,index_j)*inv_mass
            enddo
         enddo
      enddo
#ifdef _OPENACC
!$acc end kernels
!$acc wait(istream)
#endif

#ifdef _OPENACC
!$acc  exit data   copyout(EDs) ASYNC(istream)
!$acc  wait(istream)
#endif

    end subroutine col_fm_E_and_D_s


    subroutine col_fm_E_and_D_s_QLD(cs,EDs,qld_B34, qld_C34, qld_F34)
      use FP4D_globals, only : sml_2pi
      implicit none
      type(col_f_sp_type) :: cs
      real (8) :: mass,Ze
      real (8), dimension(5,col_f_nvrm1, col_f_nvzm1) :: EDs
      real (8), dimension(1:col_f_nvrm1,1:col_f_nvzm1) :: qld_B34, qld_C34, qld_F34 ! youknow
      integer :: index_I, index_J, mesh_Nrm1, mesh_Nzm1
      real (8) :: cs_mesh_r_half(lbound(cs%mesh_r_half,1):ubound(cs%mesh_r_half,1))
      real (8) :: cs_mesh_z_half(lbound(cs%mesh_z_half,1):ubound(cs%mesh_z_half,1))

#ifdef _OPENACC
      integer :: istream, ithread
      !integer, external :: omp_get_thread_num
      ithread = 0
      !ithread = omp_get_thread_num()
      istream = ithread + 1
#endif
      mass      = cs%mass
      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1
      Ze        = cs%Ze
      cs_mesh_r_half = cs%mesh_r_half
      cs_mesh_z_half = cs%mesh_z_half
#ifdef _OPENACC
!$acc  enter data pcreate(EDs) ASYNC(istream)
!$acc  parallel ASYNC(istream)                           &
!$acc& present(EDs)                                   &
!$acc& pcopyin(mesh_Nrm1,mesh_Nzm1)
!$acc loop independent  collapse(2) gang  &
#else
!$omp parallel do default(none)                                      &
!$omp& collapse(2)                                                   &
!$omp& shared(mesh_Nrm1,mesh_Nzm1,Ze,mass)       &
!$omp& shared(EDs,cs_mesh_r_half,cs_mesh_z_half)&
!$omp& shared(qld_B34,qld_C34,qld_F34,rf_pwrscale)                    &
!$omp& private(index_J,index_I)                   &
!$omp& num_threads(col_f_i_nthreads)
#endif

      do index_J=1, mesh_Nrm1
        do index_I=1,mesh_Nzm1
          EDs(1,index_J,index_I) = +sml_2pi*(1D-10)*rf_pwrscale * qld_B34(index_J,index_I) ! D_rr => -B
          EDs(2,index_J,index_I) = +sml_2pi*(1D-12)*rf_pwrscale * qld_C34(index_J,index_I) ! D_rz => -C
          EDs(3,index_J,index_I) = +sml_2pi*(1D-14)*rf_pwrscale * qld_F34(index_J,index_I) ! D_zz => -F
          EDs(4,index_J,index_I) = 0 ! youknow
          EDs(5,index_J,index_I) = 0 ! youknow
        enddo !index_I
      enddo !index_J
!$acc end parallel
!$acc wait(istream)

#ifdef _OPENACC
!$acc  exit data   copyout(EDs) ASYNC(istream)
!$acc  wait(istream)
#endif

    end subroutine col_fm_E_and_D_s_QLD


    !> Single Picard step in backward Euler scheme.
    !!
    subroutine col_fm_picard_step(iter_inter, LU_values, dist_col, dist_iter)
      implicit none
      integer :: iter_inter
      real (kind=8), dimension(LU_nnz) :: LU_values
      real (kind=8), dimension (col_f_ntotal_v) :: dist_col
      real (kind=8), dimension(col_f_nvr,col_f_nvz) :: dist_iter
      integer :: LU_info
      integer :: index_I, index_J, mat_pos, mat_pos_d, local_ij, local_i, local_j

      if( iter_inter .eq. 1 ) then

        !$omp   parallel do private(index_I,index_J) num_threads(col_f_i_nthreads)
        do index_I=1,size(dist_iter,2)
          do index_J=1,size(dist_iter,1)
             dist_iter(index_J,index_I) = 0
          enddo
        enddo

        !$omp   parallel do                                                         &
        !$omp&  default(none)                                                       &
        !$omp&  shared(col_f_nvr,col_f_nvz)                                         &
        !$omp&  shared(dist_iter,LU_values,LU_Cvalues,index_map_LU,dist_col)        &
        !$omp&  private(index_I,index_J,mat_pos,local_ij,local_i,local_j,mat_pos_d) &
        !$omp&  num_threads(col_f_i_nthreads)
        do index_I=1,col_f_nvz
            do index_J=1, col_f_nvr
               mat_pos = index_J+(index_I-1)*col_f_nvr

               !explicit time marching for 1st guess
               do local_ij=1,9
                  local_i = (local_ij-1)/3 - 1
                  local_j = mod(local_ij-1,3) - 1
                  if( index_J+local_j .gt. 0 .and. index_I+local_i .gt. 0  &
                      .and. index_J+local_j .lt. col_f_nvr+1               &
                      .and. index_I+local_i .lt. col_f_nvz+1 ) then
                      !valid mesh

                      mat_pos_d = (index_J+local_j)+((index_I+local_i)-1)*col_f_nvr
                      if( local_ij .eq. 5 ) then
                          ! (I + L(f^n) dt) f^n  : diagonal part
                          dist_iter(index_J, index_I) = dist_iter(index_J, index_I)  &
                                    +(1D0 + LU_values(LU_Cvalues(mat_pos_d))) * dist_col(mat_pos_d)
                          !print *, LU_values(LU_Cvalues(mat_pos))*vpic_tstep*vpic_gamma, vpic_tstep
                      else
                          ! L(f^n) dt f^n : off-diagonal part
                          dist_iter(index_J, index_I) = dist_iter(index_J, index_I)  &
                                    +LU_values(index_map_LU(local_ij,mat_pos)) * dist_col(mat_pos_d)
                      endif
                  endif
                  
                  ! if ((index_i.eq.10).and.(index_j.eq.10)) then
                  !   write(*,*) 'pic',local_ij,dist_iter(index_J, index_I), LU_values(index_map_LU(local_ij,mat_pos)),dist_col(mat_pos_d)
                  ! endif

               enddo   !local_ij
            enddo  !index_J
        enddo  !index_I

        ! write(*,*)  'inner iter 1', dist_iter(10, 10)
        
      else
        !IMPLICIT TIME MARCHING
        ! NOTE!!!! THAT the below "-" sign has some reasons.
        ! There can be easy mistakes that I(identity) - M does not mean
        ! 1-diagonal part. Off-diagonal part will change to -M. Be careful!

        ! LU_values = -LU_values
        !$omp   parallel do private(index_I) num_threads(col_f_i_nthreads)
        do index_I=1,size(LU_values)
          LU_values(index_I) = -LU_values(index_I)
        enddo

        LU_values(LU_cvalues) = 1D0 + LU_values(LU_cvalues)

        !Then, finally solve! Super LU!! =)
        !SUPER LU!!!! - We can directly call superLU wihtout below wrapper function!
        call bsolver(LU_n,LU_nnz,LU_nrhs,LU_values, LU_rowindx, LU_colptr, &
                     dist_col, LU_ldb, LU_info )

        if(LU_info .ne. 0) then
            write(*,*) 'LU solver : Info = ',LU_info
            write(*,*) 'It is very possible that you got NAN or INF ', &
                       'since matrix component is NAN or INF'
        endif

        !$omp   parallel do private(index_I,index_J,mat_pos) num_threads(col_f_i_nthreads)
        do index_I=1,col_f_nvz
          do index_J=1,col_f_nvr
            mat_pos = index_J+(index_I-1)*col_f_nvr
            dist_iter(index_J,index_I) = dist_col( mat_pos )
          enddo
        enddo
      endif

      return

    end subroutine col_fm_picard_step

    subroutine cal_power_partition(cs, LU_values_save, dist_col, iter_inter)

      ! youknow :: the dimension of result is [mv^2 * C(f) * dt]

      implicit none
      type(col_f_sp_type) :: cs
      !real (8) :: pow_spj
      !! dist_n <- pdf_n, dist_iter <- pdf_np1
      !!real (kind=8), dimension(col_f_nvr, col_f_nvz) :: dist_iter   ! local
      real (8) :: vpic_dw !! were output
      real (8) :: vpic_dw_0 !! were output ! youknow : [C(f) * dt]

      integer :: index_SPJ
      real (8) :: tmpr1, tmpr12, tmpr2, tmpr3, dt
      real (kind=8), dimension(LU_nnz,tar_f_sp_s:tar_f_sp_e) :: LU_values_save
      real (kind=8), dimension (col_f_ntotal_v) :: dist_col
      real (kind=8), dimension(col_f_nvr,col_f_nvz) :: col_term
      integer :: LU_info
      integer :: index_I, index_J, mat_pos, mat_pos_d, local_ij, local_i, local_j
      integer :: iter_inter
      ! DIST_ITER has new updated PDF.


      do index_spj=tar_f_sp_s, tar_f_sp_e
        vpic_dw = 0D0
        vpic_dw_0 = 0D0

        do index_I=1,size(col_term,2)
          do index_J=1,size(col_term,1)
             col_term(index_J,index_I) = 0
          enddo
        enddo

       ! calculte each collision term for spj
        do index_I=1,col_f_nvz

            tmpr1= cs%mesh_z(index_I)
            tmpr12=tmpr1*tmpr1
    
            do index_J=1, col_f_nvr
              tmpr2 = cs%mesh_r(index_J)
              tmpr2 = tmpr2*tmpr2     ! mesh_r^2
              tmpr3 = tmpr12+tmpr2     ! mesh_z^2+mesh_r^2
             
                    ! col_term(index_J,index_I)=cs%pdf_np1(index_J,index_I) 
                    ! col_term(index_J,index_I)=0D0 !youknow
              mat_pos = index_J+(index_I-1)*col_f_nvr

               !explicit time marching for 1st guess
               do local_ij=1,9
                  local_i = (local_ij-1)/3 - 1
                  local_j = mod(local_ij-1,3) - 1
                  if( index_J+local_j .gt. 0 .and. index_I+local_i .gt. 0  &
                      .and. index_J+local_j .lt. col_f_nvr+1               &
                      .and. index_I+local_i .lt. col_f_nvz+1 ) then
                      !valid mesh

                      mat_pos_d =(index_J+local_j)+((index_I+local_i)-1)*col_f_nvr
                      if( local_ij .eq. 5 ) then
                          ! (I + L(f^n) dt) f^n  : diagonal part
                          col_term(index_J,index_I)=col_term(index_J,index_I)+ dist_col(mat_pos_d) &!cs%pdf_np1(index_J,index_I) &
                                    ! *(1D0 + LU_values_save(LU_Cvalues(mat_pos_d),index_spj)) 
                                    *(0D0 + LU_values_save(LU_Cvalues(mat_pos_d),index_spj)) 
                          !print *,
                          !LU_values(LU_Cvalues(mat_pos))*vpic_tstep*vpic_gamma,
                          !vpic_tstep
                      else
                          ! L(f^n) dt f^n : off-diagonal part
                          col_term(index_J,index_I)=col_term(index_J,index_I)+dist_col(mat_pos_d) &!cs%pdf_np1(index_J,index_I)  &
                                    *LU_values_save(index_map_LU(local_ij,mat_pos),index_spj)
                      endif
                  endif

                  ! if ((index_i.eq.10).and.(index_j.eq.10)) then
                  !   write(*,*) 'pow',local_ij,col_term(index_J, index_I), LU_values_save(index_map_LU(local_ij,mat_pos),index_spj),dist_col(mat_pos_d)

                  ! endif
               enddo   !local_ij

                ! vpic_dw=vpic_dw+ col_term(index_J,index_I)*tmpr3 !JPL
              vpic_dw = vpic_dw + col_term(index_J,index_I) * cs%vol(index_J) * tmpr3
              vpic_dw_0 = vpic_dw_0 + col_term(index_J,index_I) * cs%vol(index_J) ! youknow

            enddo  !index_J
        enddo  !index_I
        ! write(*,*) 'pdf_np1', cs%pdf_np1(10,10)

        ! write(*,*) 'col_term', col_term(10,10)
        ! cs%power(index_spj)=vpic_dw*cs%mass/dt !JPL
        cs%power(index_spj) = vpic_dw * cs%mass
        cs%col_dens(index_spj) = vpic_dw_0

        ! if (i_opt_col_inner_iter.eq.1) then
          cs%inner_info(index_spj, 0, iter_inter) = vpic_dw_0
          cs%inner_info(index_spj, 1, iter_inter) = 0D0 ! tmp.
          cs%inner_info(index_spj, 2, iter_inter) = vpic_dw * cs%mass
        ! endif

        ! if (sml_mype.eq.0) then
        !   write(*,*) 'iter_inter', iter_inter
        !   write(*,*) 'cs%power(index_spj)', cs%power(index_spj)
        !   write(*,*) 'cs%inner_info(index_spj, 2, iter_inter)', cs%inner_info(index_spj, 2, iter_inter) 
        ! endif



      enddo  !index_spj

    end subroutine

    !> Checks the convergence criterion (i.e. conservation of
    !! mass, energy and momentum.
    !!
    subroutine col_fm_convergence_eval(cs, col_dw, col_dp, col_dn_n)
      implicit none
      type(col_f_sp_type) :: cs
      real (8) :: col_dw, col_dp, col_dn_n
      !! dist_n <- pdf_n, dist_iter <- pdf_np1
      !!real (kind=8), dimension(col_f_nvr, col_f_nvz) :: dist_iter   ! local
      real (8) :: vpic_dn, vpic_dw, vpic_dp !! were output
      real (8) :: vpic_dfc
      integer :: index_I, index_J
      real (8) :: tmpr1, tmpr12, tmpr2, tmpr3
      ! DIST_ITER has new updated PDF.

      vpic_dn = 0D0
      vpic_dw = 0D0
      vpic_dp = 0D0

      !--------------------------------
      !efd add more reduction variables
      !--------------------------------
      !$omp  parallel do &
      !$omp& shared(cs)  &
      !$omp& private(index_I,index_J,tmpr1,tmpr12,tmpr2,tmpr3,vpic_dfc) &
      !$omp& reduction(+:vpic_dn,vpic_dw,vpic_dp)            &
      !$omp& num_threads(col_f_i_nthreads)
      do index_I=1,col_f_nvz
          tmpr1 = cs%mesh_z(index_I)
          tmpr12 = tmpr1*tmpr1     ! mesh_z^2
          do index_J=1,col_f_nvr
              tmpr2 = cs%mesh_r(index_J)
              tmpr2 = tmpr2*tmpr2     ! mesh_r^2
              tmpr3 = tmpr12+tmpr2     ! mesh_z^2+mesh_r^2
              ! youknow error debug
              vpic_dfc = (cs%pdf_np1(index_J,index_I) - cs%pdf_n(index_J, index_I))*cs%vol(index_J)
              ! write(*,*) 'debug :: need to fix vpic_dif after debug'
              ! vpic_dfc = (cs%pdf_np1(index_J,index_I) )*cs%vol(index_J)

              vpic_dn = vpic_dn + vpic_dfc
              vpic_dw = vpic_dw + vpic_dfc * tmpr3
              vpic_dp = vpic_dp + vpic_dfc * cs%mesh_z(index_I)
          enddo
      enddo

      ! write(*,*) 'converge', (cs%pdf_np1(10,10) - cs%pdf_n(10, 10))

      col_dw     = vpic_dw*cs%mass
      col_dp     = vpic_dp*cs%mass
      col_dn_n   = dabs(vpic_dn/cs%dens)

      return

    end subroutine col_fm_convergence_eval


    !> Computes the shift for inter-grid quantities (Sec. 2.5 in Hager et al.)
    !!
    subroutine col_fm_core_delta_init(spa, spb)
      use FP4D_globals, only : sml_e_charge     !! XGC1_3
      implicit none
      type(col_f_sp_type), intent(inout) :: spa, spb
      real (kind=8) :: maxw_fac !, en_t, maxw(4), alpha(4)
      real (kind=8) :: lambda_r, lambda_z, delta(2) !, a, b, c, d, beta
      real (kind=8) :: arg
      ! Limits argument to exp(:), exp(227.0) = 3.8E98
      real (kind=8), parameter :: arg_max = 227.0d0
      real (kind=8), parameter :: arg_min = -arg_max
      integer :: mesh_Nrm1, mesh_Nzm1, index_i, index_j
      real (kind=8) :: numeric_teq, numeric_vtheq2

      mesh_Nrm1 = col_f_nvr-1
      mesh_Nzm1 = col_f_nvz-1

      if(spa%spid .eq. spb%spid) then !! ### ----- Same-species collisions

        !! ### ----- R-direction
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh_Nrm1,spa) &
        !$omp private(index_j,lambda_r,arg,maxw_fac,delta) num_threads(col_f_i_nthreads)
        do index_j=1, mesh_Nrm1
           ! Same-species collisions
           lambda_r = -spa%mesh_r_half(index_j)/spa%numeric_vth2 * spa%mesh_dr
           arg = (spa%mesh_dr*(spa%mesh_dr+2D0*spa%mesh_r(index_j))/ &
                      (2D0*spa%numeric_vth2))
           arg = max(arg_min,min(arg_max, arg))
           if (abs(arg) <= 1.0d-6) then
             maxw_fac = -arg*(1.0d0+arg/2.0d0*(1.0d0+arg/3.0d0))
           else
             maxw_fac=1D0-exp(arg)
           endif
           delta(1)=1D0/maxw_fac-1D0/lambda_r
           delta(2)=1D0/maxw_fac

           if (0D0 .lt. delta(1) .and. delta(1) .lt. 1) then
              spa%delta_r(index_j,spa%spid)=delta(1)
           elseif (0D0 .lt. delta(2) .and. delta(2) .lt. 1) then
              spa%delta_r(index_j,spa%spid)=delta(2)
           else
              write(*,'(a)') 'COL_F_CORE_DELTA_INIT: Error (1a-fm) in calculation   &
                              of finite difference parameter!'
              write(*,'(i8,5(e22.12,1x))') index_j, maxw_fac, lambda_r, delta(1), delta(2), arg
              write(*,'(4(e22.12,1x))') spa%mesh_r(index_j), spa%mesh_r_half(index_j),  &
                                        spa%numeric_vth2, spa%mesh_dr
              write(*,'(2(e22.12,1x))')  spa%numeric_t, spa%dens
              stop
           endif
        enddo

        !! ### ----- Z-direction
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh_Nzm1,spa) &
        !$omp private(index_i,lambda_z,arg,maxw_fac,delta) num_threads(col_f_i_nthreads)
        do index_i=1, mesh_Nzm1
           lambda_z = -spa%mesh_z_half(index_i)/spa%numeric_vth2 * spa%mesh_dz
           arg = (spa%mesh_dz*(spa%mesh_dz+2D0*spa%mesh_z(index_i))/ &
                      (2D0*spa%numeric_vth2))
           arg = max(arg_min,min(arg_max, arg))
           if (abs(arg) <= 1.0d-6) then
             maxw_fac = -arg*(1.0d0+arg/2.0d0*(1.0d0+arg/3.0d0))
           else
             maxw_fac=1D0-exp(arg)
           endif
           if(abs(spa%mesh_z_half(index_i)) .lt. 0.25D0*spa%mesh_dz) then
              delta(1) = 0.5D0
              delta(2) = 0.5D0
           else
              delta(1) = (lambda_z - maxw_fac)/(maxw_fac*lambda_z)
              delta(2)=1D0/maxw_fac
           endif
           if (0D0 .lt. delta(1) .and. delta(1) .lt. 1) then
              spa%delta_z(index_i,spa%spid)=delta(1)
           elseif (0D0 .lt. delta(2) .and. delta(2) .lt. 1) then
              spa%delta_z(index_i,spa%spid)=delta(2)
           else
              write(*,'(a)') 'COL_F_DELTA_INIT: Error (2a-fm) in calculation  &
                              of finite difference parameter!'
              write(*,'(i8,8(e22.12,1x))') index_i, maxw_fac, lambda_z, delta(1), delta(2),  &
                            spa%mesh_z_half(index_i), spa%numeric_vth2, spa%mesh_dz, arg
              stop
           endif
        enddo

      else   !! ### ----- Inter-species collisions

        numeric_teq    = (spa%den_lambda_gamma*spa%numeric_T + spb%den_lambda_gamma*spb%numeric_T)  &
                         /(spa%den_lambda_gamma+spb%den_lambda_gamma)
        numeric_vtheq2 = numeric_teq*sml_e_charge/spa%mass !  sp%mass --> ptl_mass !! ## ----- this is also problematic correspondingly.

        !! ### ----- R-direction
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh_Nrm1,spa,spb,numeric_vtheq2) &
        !$omp private(index_j,lambda_r,arg,maxw_fac,delta) num_threads(col_f_i_nthreads)
        do index_j=1, mesh_Nrm1
           lambda_r = -spa%mesh_r_half(index_j)/numeric_vtheq2 * spa%mesh_dr
           arg = (spa%mesh_dr*(spa%mesh_dr+2D0*spa%mesh_r(index_j))/ &
                     (2D0*numeric_vtheq2))
           arg = max(arg_min,min(arg_max, arg))
           if (abs(arg) <= 1.0d-6) then
             maxw_fac = -arg*(1.0d0+arg/2.0d0*(1.0d0+arg/3.0d0))
           else
             maxw_fac=1D0-exp(arg)
           endif
           delta(1)=1D0/maxw_fac-1D0/lambda_r
           delta(2)=1D0/maxw_fac
           if (0D0 .lt. delta(1) .and. delta(1) .lt. 1) then
              spa%delta_r(index_j,spb%spid)=delta(1)
           elseif (0D0 .lt. delta(2) .and. delta(2) .lt. 1) then
              spa%delta_r(index_j,spb%spid)=delta(2)
           else
              write(*,'(a)') 'COL_F_DELTA_INIT: Error (1b-fm) in calculation  &
                              of finite difference parameter!'
              write(*,'(i8,4(e22.12,1x))') index_j, maxw_fac, lambda_r, delta(1), delta(2)
              stop
           endif
        enddo

        !! ### ----- Z-direction
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh_Nzm1,spa,spb,numeric_vtheq2) &
        !$omp private(index_i,lambda_z,arg,maxw_fac,delta) num_threads(col_f_i_nthreads)
        do index_i=1, mesh_Nzm1
          lambda_z = -spa%mesh_z_half(index_i)/numeric_vtheq2 * spa%mesh_dz
          arg = (spa%mesh_dz*(spa%mesh_dz+2D0*spa%mesh_z(index_i))/ &
                (2D0*numeric_vtheq2))
          arg = max(arg_min, min(arg_max, arg ))
          if (abs(arg) <= 1.0d-6) then
            maxw_fac = -arg*(1.0d0+arg/2.0d0*(1.0d0+arg/3.0d0))
          else
            maxw_fac = 1.0d0  - exp(arg)
          endif

          !if(abs(cs%mesh_z_half(index_i)) .lt. epsilon(0D0)) then
          if(abs(spa%mesh_z_half(index_i)) .lt. 0.25D0*spa%mesh_dz) then
             delta(1) = 0.5D0
             delta(2) = 0.5D0
          else
             delta(1) = (lambda_z - maxw_fac)/(maxw_fac*lambda_z)
             delta(2)=1D0/maxw_fac
          endif
          if (0D0 .lt. delta(1) .and. delta(1) .lt. 1) then
             spa%delta_z(index_i,spb%spid)=delta(1)
          elseif (0D0 .lt. delta(2) .and. delta(2) .lt. 1) then
             spa%delta_z(index_i,spb%spid)=delta(2)
          else
             write(*,'(a)') 'COL_F_DELTA_INIT: Error (2b-fm) in calculation  &
                             of finite difference parameter!'
             write(*,'(3(i8,1x),8(e22.12,1x))') spa%spid, spb%spid, index_i, maxw_fac,  &
                      lambda_z, delta(1), delta(2),  &
                      spa%mesh_z_half(index_i), numeric_vtheq2, spa%mesh_dz
             stop
          endif
        enddo
      endif
    end subroutine col_fm_core_delta_init


    ! op_mode==1 for same species, op_mode==2 for different species collision
    ! cell_I and cell_J seem to be cell index  (c.f. grid index)
    ! cs is required for cs%mesh_r_half , cs%mesh_dr/dz
    ! mat_pos_rel_indx
    subroutine col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs_comp, LU_values)
      implicit none
      integer :: op_mode
      type(col_f_sp_type) :: cs
      integer :: mat_pos_rel_indx(4)
      integer :: mat_pos, cell_I, cell_J
      real (8) :: coeff1, coeff2, coeff_loc1, coeff_loc2
      real (8) :: EDs_comp(5)
      real (kind=8), dimension(LU_nnz) :: LU_values
      real (8) :: M_rr_fOVERdr, M_rz_fOVERdz, M_zz_fOVERdz, M_zr_fOVERdr, tmp_ER_sum, tmp_EZ_sum, delr, cdelr, delz, cdelz

      coeff_loc1 = coeff1*cs%mesh_r_half(cell_J)
      coeff_loc2 = coeff2*cs%mesh_r_half(cell_J)
      ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
      M_rr_fOVERdr = EDs_comp(1)/cs%mesh_dr
      M_rz_fOVERdz = EDs_comp(2)/cs%mesh_dz
      M_zr_fOVERdr = EDs_comp(2)/cs%mesh_dr
      M_zz_fOVERdz = EDs_comp(3)/cs%mesh_dz
      tmp_ER_sum   = EDs_comp(4)
      tmp_EZ_sum   = EDs_comp(5)

      delr = cs%delta_r(cell_J,op_mode)
      ! delr = 0.5 ! youknow
      ! if (cell_J.eq.1) then
      !   delr = 0.0001 ! youknow
      ! elseif (cell_J.eq.2) then
      !   delr = cs%delta_r(cell_J,op_mode)+cs%delta_r(cell_J-1,op_mode)-0.0001
      ! endif

      cdelr = (1D0-delr)
      delz = cs%delta_z(cell_I,op_mode)
      cdelz = (1D0-delz)

      !for (I,J)
      LU_values(index_map_LU(mat_pos_rel_indx(1),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(1),mat_pos)) &
                            + coeff_loc1* (tmp_ER_sum*delr*delz - M_rr_fOVERdr*(-delz) - M_rz_fOVERdz*(-delr)) &
                            + coeff_loc2* (tmp_EZ_sum*delr*delz - M_zr_fOVERdr*(-delz) - M_zz_fOVERdz*(-delr))
      !for (I,J+1)
      LU_values(index_map_LU(mat_pos_rel_indx(2),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(2),mat_pos)) &
                            + coeff_loc1* (tmp_ER_sum*cdelr*delz - M_rr_fOVERdr*(delz) - M_rz_fOVERdz*(-cdelr)) &
                            + coeff_loc2* (tmp_EZ_sum*cdelr*delz - M_zr_fOVERdr*(delz) - M_zz_fOVERdz*(-cdelr))
      !for(I+1,J)
      LU_values(index_map_LU(mat_pos_rel_indx(3),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(3),mat_pos)) &
                            + coeff_loc1* (tmp_ER_sum*delr*cdelz - M_rr_fOVERdr*(-cdelz) - M_rz_fOVERdz*(delr)) &
                            + coeff_loc2* (tmp_EZ_sum*delr*cdelz - M_zr_fOVERdr*(-cdelz) - M_zz_fOVERdz*(delr))
      !for(I+1,J+1)
      LU_values(index_map_LU(mat_pos_rel_indx(4),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(4),mat_pos)) &
                            + coeff_loc1* (tmp_ER_sum*cdelr*cdelz - M_rr_fOVERdr*(cdelz) - M_rz_fOVERdz*(cdelr)) &
                            + coeff_loc2* (tmp_EZ_sum*cdelr*cdelz - M_zr_fOVERdr*(cdelz) - M_zz_fOVERdz*(cdelr))
      return

    end subroutine col_f_LU_matrix_ftn

    subroutine col_f_LU_matrix_ftn_youknow(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs_comp, LU_values)
      implicit none
      integer :: op_mode
      type(col_f_sp_type) :: cs
      integer :: mat_pos_rel_indx(4)
      integer :: mat_pos, cell_I, cell_J
      real (8) :: coeff1, coeff2, coeff_loc1, coeff_loc2
      real (8) :: EDs_comp(5)
      real (kind=8), dimension(LU_nnz) :: LU_values
      real (8) :: M_rr_fOVERdr, M_rz_fOVERdz, M_zz_fOVERdz, M_zr_fOVERdr, tmp_ER_sum, tmp_EZ_sum, delr, cdelr, delz, cdelz

      coeff_loc1 = coeff1*cs%mesh_r_half(cell_J)
      coeff_loc2 = coeff2*cs%mesh_r_half(cell_J)
      ! EDs = 1: Drr, 2: Drz, 3:Dzz, 4:Dzr, 5: Er, 6:Ez
      M_rr_fOVERdr = EDs_comp(1)/cs%mesh_dr
      M_rz_fOVERdz = EDs_comp(2)/cs%mesh_dz
      M_zr_fOVERdr = EDs_comp(2)/cs%mesh_dr
      M_zz_fOVERdz = EDs_comp(3)/cs%mesh_dz
      tmp_ER_sum   = EDs_comp(4)
      tmp_EZ_sum   = EDs_comp(5)

      delr = cs%delta_r(cell_J,op_mode)
      ! delr = 0.5 ! youknow

      cdelr = (1D0-delr)
      delz = cs%delta_z(cell_I,op_mode)
      cdelz = (1D0-delz)

      ! youknow
      coeff_loc1 = coeff1*cs%mesh_r_half(cell_J) / delr
      coeff_loc2 = coeff2*cs%mesh_r_half(cell_J) / delr

      !for (I,J)
      LU_values(index_map_LU(mat_pos_rel_indx(1),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(1),mat_pos)) &
                            + coeff_loc1* (tmp_ER_sum*delr*delz - M_rr_fOVERdr*(-delz) - M_rz_fOVERdz*(-delr)) &
                            + coeff_loc2* (tmp_EZ_sum*delr*delz - M_zr_fOVERdr*(-delz) - M_zz_fOVERdz*(-delr))
      ! !for (I,J+1)
      ! LU_values(index_map_LU(mat_pos_rel_indx(2),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(2),mat_pos)) &
      !                       + coeff_loc1* (tmp_ER_sum*cdelr*delz - M_rr_fOVERdr*(delz) - M_rz_fOVERdz*(-cdelr)) &
      !                       + coeff_loc2* (tmp_EZ_sum*cdelr*delz - M_zr_fOVERdr*(delz) - M_zz_fOVERdz*(-cdelr))
      !for(I+1,J)
      LU_values(index_map_LU(mat_pos_rel_indx(3),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(3),mat_pos)) &
                            + coeff_loc1* (tmp_ER_sum*delr*cdelz - M_rr_fOVERdr*(-cdelz) - M_rz_fOVERdz*(delr)) &
                            + coeff_loc2* (tmp_EZ_sum*delr*cdelz - M_zr_fOVERdr*(-cdelz) - M_zz_fOVERdz*(delr))
      ! !for(I+1,J+1)
      ! LU_values(index_map_LU(mat_pos_rel_indx(4),mat_pos)) = LU_values(index_map_LU(mat_pos_rel_indx(4),mat_pos)) &
      !                       + coeff_loc1* (tmp_ER_sum*cdelr*cdelz - M_rr_fOVERdr*(cdelz) - M_rz_fOVERdz*(cdelr)) &
      !                       + coeff_loc2* (tmp_EZ_sum*cdelr*cdelz - M_zr_fOVERdr*(cdelz) - M_zz_fOVERdz*(cdelr))
      return

    end subroutine col_f_LU_matrix_ftn_youknow

    subroutine col_f_LU_matrix(op_mode, cs, EDs, LU_values)
      ! EDs : target
      implicit none
      integer :: op_mode
      type(col_f_sp_type) :: cs
      real (8), dimension(5,col_f_nvrm1,col_f_nvzm1) :: EDs
      real (kind=8), dimension(LU_nnz) :: LU_values
      integer :: index_I, index_J, cell_I, cell_J, mesh_Nr,mesh_Nz
      integer :: mat_pos_rel_indx(4)
      real (8) :: coeff1_ab, coeff2_ab, coeff1, coeff2, local_vol
      real (8) :: M_rr_fOVERdr, M_rz_fOVERdz, M_zz_fOVERdz, M_zr_fOVERdr, tmp_ER_sum, tmp_EZ_sum, delr, cdelr, delz, cdelz
      integer :: local_ij, local_i, local_j, mat_pos, mat_pos_d

      mesh_Nr = col_f_nvr
      mesh_Nz = col_f_nvz

      LU_values(:) = 0D0
      do index_I=1,mesh_Nz
        do index_J=1, mesh_Nr
        !!! HJL Parallel Streaming B.C!!
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          mat_pos = index_J+(index_I-1)*col_f_nvr
             ! What we know : the number of cells related to cell (J,I) from LU_colptr
             !                index number of relevant cells from LU_rowindx
      
             ! NOTE that local_vol is for volume in LHS. Therefore, fixed values!
          local_vol = cs%vol(index_J)
          coeff1_ab    = 0.5D0*cs%mesh_dz/local_vol
          coeff2_ab    = 0.5D0*cs%mesh_dr/local_vol

          if((index_I .ne. mesh_Nz) .and. (index_J .ne. mesh_Nr)) then
          ! Existence of (I+1/2, J+1/2)
              cell_I = index_I
              cell_J = index_J
              mat_pos_rel_indx=(/5, 6, 8, 9/) !just reuse array
              coeff1 = -coeff1_ab
              coeff2 = -coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif
          if((index_I .ne. 1) .and. (index_J .ne. mesh_Nr)) then
          ! Existence of (I-1/2, J+1/2)
              cell_I = index_I-1
              cell_J = index_J
              mat_pos_rel_indx=(/2, 3, 5, 6/)
              coeff1 = -coeff1_ab
              coeff2 =  coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values) 
          endif

          if((index_I .ne. mesh_Nz) .and. (index_J .ne. 1)) then
          ! Existence of (I+1/2, J-1/2)
              cell_I = index_I
              cell_J = index_J-1
              mat_pos_rel_indx=(/4, 5, 7, 8/)
              coeff1 =  coeff1_ab
              coeff2 = -coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif

          if( (index_I .ne. 1) .and. (index_J .ne. 1) ) then
          ! Existence of (I-1/2, J-1/2)
              cell_I = index_I-1
              cell_J = index_J-1
              mat_pos_rel_indx=(/1, 2, 4, 5/)
              coeff1 = coeff1_ab
              coeff2 = coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif
        enddo !index_J
      enddo !index_I


    end subroutine col_f_LU_matrix


    subroutine col_f_LU_matrix_QLD(op_mode, cs, EDs, LU_values)
      ! EDs : target
      implicit none
      integer :: op_mode
      type(col_f_sp_type) :: cs
      real (8), dimension(5,col_f_nvrm1,col_f_nvzm1) :: EDs
      real (kind=8), dimension(LU_nnz) :: LU_values
      integer :: index_I, index_J, cell_I, cell_J, mesh_Nr,mesh_Nz
      integer :: mat_pos_rel_indx(4)
      real (8) :: coeff1_ab, coeff2_ab, coeff1, coeff2, local_vol
      real (8) :: M_rr_fOVERdr, M_rz_fOVERdz, M_zz_fOVERdz, M_zr_fOVERdr, tmp_ER_sum, tmp_EZ_sum, delr, cdelr, delz, cdelz
      integer :: local_ij, local_i, local_j, mat_pos, mat_pos_d

      mesh_Nr = col_f_nvr
      mesh_Nz = col_f_nvz

      LU_values(:) = 0D0
      do index_I=1,mesh_Nz
        do index_J=1, mesh_Nr
        !!! HJL Parallel Streaming B.C!!
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          mat_pos = index_J+(index_I-1)*col_f_nvr
             ! What we know : the number of cells related to cell (J,I) from LU_colptr
             !                index number of relevant cells from LU_rowindx
      
             ! NOTE that local_vol is for volume in LHS. Therefore, fixed values!
          local_vol = cs%vol(index_J)
          coeff1_ab    = 0.5D0*cs%mesh_dz/local_vol
          coeff2_ab    = 0.5D0*cs%mesh_dr/local_vol

          if((index_I .ne. mesh_Nz) .and. (index_J .ne. mesh_Nr)) then
          ! Existence of (I+1/2, J+1/2)
              cell_I = index_I
              cell_J = index_J
              mat_pos_rel_indx=(/5, 6, 8, 9/) !just reuse array
              coeff1 = -coeff1_ab
              coeff2 = -coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif
          if((index_I .ne. 1) .and. (index_J .ne. mesh_Nr)) then
          ! Existence of (I-1/2, J+1/2)
              cell_I = index_I-1
              cell_J = index_J
              mat_pos_rel_indx=(/2, 3, 5, 6/)
              coeff1 = -coeff1_ab
              coeff2 =  coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values) 
          endif

          if((index_I .ne. mesh_Nz) .and. (index_J .ne. 1)) then
          ! Existence of (I+1/2, J-1/2)
              cell_I = index_I
              cell_J = index_J-1
              mat_pos_rel_indx=(/4, 5, 7, 8/)
              coeff1 =  coeff1_ab
              coeff2 = -coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif

          if( (index_I .ne. 1) .and. (index_J .ne. 1) ) then
          ! Existence of (I-1/2, J-1/2)
              cell_I = index_I-1
              cell_J = index_J-1
              mat_pos_rel_indx=(/1, 2, 4, 5/)
              coeff1 = coeff1_ab
              coeff2 = coeff2_ab
              call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!    !!! @ J = 1. additional...

          if((index_I .ne. mesh_Nz) .and. (index_J .eq. 1)) then
            ! Existence of (I+1/2, J+1/2)
                cell_I = index_I
                cell_J = index_J
                mat_pos_rel_indx=(/5, 6, 8, 9/) !just reuse array
                coeff1 = (+1.0/4.0)*-coeff1_ab              ! -coeff1_ab
                coeff2 = 0                                  ! -coeff2_ab
                call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif
          if((index_I .ne. 1) .and. (index_J .eq. 1)) then
            ! Existence of (I-1/2, J+1/2)
                cell_I = index_I-1
                cell_J = index_J
                mat_pos_rel_indx=(/2, 3, 5, 6/)
                coeff1 = (+1.0/4.0)*-coeff1_ab              ! -coeff1_ab
                coeff2 = 0                                  ! coeff2_ab
                call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values) 
          endif

          if((index_I .ne. mesh_Nz) .and. (index_J .eq. 2)) then
              ! Existence of (I+1/2, J+1/2)
                cell_I = index_I
                cell_J = index_J
                mat_pos_rel_indx=(/5, 6, 8, 9/) !just reuse array
                coeff1 = (-1.0/12.0)*-coeff1_ab             ! -coeff1_ab
                coeff2 = 0                                  ! -coeff2_ab
                call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values)
          endif
          if((index_I .ne. 1) .and. (index_J .eq. 2)) then
              ! Existence of (I-1/2, J+1/2)
                cell_I = index_I-1
                cell_J = index_J
                mat_pos_rel_indx=(/2, 3, 5, 6/)
                coeff1 = (-1.0/12.0)*-coeff1_ab             ! -coeff1_ab
                coeff2 = 0                                  ! coeff2_ab
                call col_f_LU_matrix_ftn(op_mode,cell_I, cell_J, cs, mat_pos_rel_indx, mat_pos, coeff1, coeff2, EDs(:,cell_J,cell_I), LU_values) 
          endif

        enddo !index_J
      enddo !index_I

    end subroutine col_f_LU_matrix_QLD

    subroutine col_dt_gt(n)

      implicit none
      integer,intent(in) :: n
      col_f_dt=col_f_dt*real(n,8)

    end subroutine col_dt_gt


    subroutine QLD_dt_gt(n)

      implicit none
      integer,intent(in) :: n
      QLD_f_dt=QLD_f_dt*real(n,8)

    end subroutine QLD_dt_gt

end module col_f_module
