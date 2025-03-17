#ifdef USE_ASYNC
#define ASYNC(istream)  async(istream)
#else
#define ASYNC(istream)
#endif

module quasi_f_module
  use FP4D_globals
  use FP4D_math


  integer :: col_f_i_nthreads                     !! [col_f_setup] the number of inner OMP threads
  integer :: col_f_o_nthreads                     !! [col_f_setup] the number of outer OMP threads


  contains
    ! BLAS solver collision operator time-stepping
#include "bsolver.F90"

    !! @param[in]   mype   MPI rank ID
    !!
    subroutine quasi_f_setup(mype)
      use FP4D_globals, only: sml_charge_num
      implicit none
      integer, intent(in) :: mype

      real(kind=8) :: max_n0
      integer :: idx_max_n0
      
      integer :: is, ip, ith
   
      col_f_i_nthreads = col_f_nthreads
      call set_omp_outer_thread_num( col_f_o_nthreads )
      if(mype .eq. 0) then
          write(*,'(a,i4)') '[MULTI-PRL] inner thread = ', col_f_i_nthreads
          write(*,'(a,i4)') '[MULTI-PRL] outer thread = ', col_f_o_nthreads
      endif

      ! youknow - FP4D_v1.1
      opt_Potential = 2 ! basically adiabatic electron

      do is=tar_s,tar_e
        if ( sml_charge_num(is) .lt. 0D0 ) then
          opt_Potential = 1
        endif
      enddo

      allocate(n0_ae(ith_s:ith_e, ip_s:ip_e)) ! for rotation case
      allocate(T0_ae(ip_s:ip_e))

      n0_ae(:,:) = 0D0
      T0_ae(:) = 0D0

      if ( opt_Potential .eq. 2 ) then
        do ip=ip_s,ip_e
          do ith=ith_s,ith_e
              do is=tar_s,tar_e
                  n0_ae(ith,ip) = n0_ae(ith,ip) + sml_charge_num(is) * small_n0(ith,ip,is) ! rotation
              enddo
          enddo

          ! Assume that T_0e follows the T0 of most popular ion
              ! This part only searches for idx, so there is no issue even during rotation.
          do is=tar_s,tar_e
            if ( eq_den(ip,is) .gt. max_n0 ) then
              max_n0 = eq_den(ip,is)
              idx_max_n0 = is
            endif
          enddo
          T0_ae(ip) = eq_t_ev(ip,idx_max_n0)
        enddo
      endif

      write(*,*) 'opt_Potential', opt_Potential
      write(*,*) 'n0_ae', n0_ae
      write(*,*) 'T0_ae', T0_ae

      return

    end subroutine quasi_f_setup

  !> Deallocate quasi memory
    subroutine quasi_f_finalize
      implicit none

      if(allocated(n0_ae))          deallocate(n0_ae)
      if(allocated(T0_ae))          deallocate(T0_ae)

      return

    end subroutine quasi_f_finalize

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





    !! This routine performs the loop over configuration space grid points with
    !! OpenMP parallelization --> "Outer level" parallelism.
    !!
    subroutine fm_quasi(g1_in, g1_out, f1_in, f0_in,nthreads,&
                phi1_in, phi1_out, dphidth)
      use FP4D_globals, only: itmax
      use FP4D_globals, only: sml_charge_num
      use FP4D_globals, only: sml_2pi
      use FP4D_globals, only : g1_tar, f1_tar

#ifdef _OPENMP
      use omp_lib
#endif
      !use mpi
      implicit none
      
      integer, intent(in) :: nthreads
      real (kind=8), intent(in) :: g1_in (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e)
      real (kind=8), intent(out) :: g1_out (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e) !youknow
      real (kind=8), intent(in) :: f1_in   (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e)
      real (kind=8), intent(in) :: f0_in   (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e)

      real(kind=8), intent(in)  :: phi1_in (ith_s:ith_e,ip_s:ip_e)
      real(kind=8), intent(out) :: phi1_out(ith_s:ith_e,ip_s:ip_e)
      real(kind=8), intent(out) :: dphidth (ith_s:ith_e,ip_s:ip_e)

      ! let's assign as global vars in FP4D_globals and save...
      real(kind=8) :: result_abs (ith_s:ith_e, ip_s:ip_e)
      real(kind=8) :: result_rel (ith_s:ith_e, ip_s:ip_e)
      real(kind=8) :: result_iter(ith_s:ith_e, ip_s:ip_e)

      ! Parallelization for (ith,ip).
      ! vars for phi_iter 
      real(kind=8) :: g1_3D (iz_s:iz_e, ir_s:ir_e, is_s:is_e)
      real(kind=8) :: f0_3D (iz_s:iz_e, ir_s:ir_e, is_s:is_e)
      real(kind=8) :: f1_3D (iz_s:iz_e, ir_s:ir_e, is_s:is_e)
      real(kind=8) :: g1_3D_tar (iz_s:iz_e, ir_s:ir_e, tar_s:tar_e) ! opt_restart=3
      real(kind=8) :: f1_3D_tar (iz_s:iz_e, ir_s:ir_e, tar_s:tar_e) ! opt_restart=3

      real(kind=8) :: n1_g1(is_s:is_e)
      ! real(kind=8) :: n1_f1(is_s:is_e)
      real(kind=8) :: n1_g1_tar(tar_s:tar_e)
      real(kind=8) :: n1_f1_tar(tar_s:tar_e)
      

      real(kind=8) :: abs_error, rel_error
      real(kind=8) :: phi1_abs_tol, phi1_rel_tol ! Need to transfer to global var & input var

      real(kind=8) :: numer, denom    ! youknow
      real(kind=8) :: phi1 , phi1_old ! youknow

      logical :: absolute_ok, relative_ok, update_ok

      integer :: alloc_stat
      integer :: is, ip, ith, ir, iz, phi_iter


      real (kind=8) :: conv_factor



      opt_Potential_converge = 1  ! Need to transfer to global var & input var

      dphidth=0d0

      phi1_abs_tol = 1D-14  ! 나중에 read 하도록 바꾸기
      phi1_rel_tol = 1D-10   ! 나중에 read 하도록 바꾸기

      !$omp parallel default (none) num_threads(col_f_i_nthreads)&
      !$omp shared(g1_in, f0_in, f1_in, g1_out) &
      !$omp shared(g1_3D, f0_3D, f1_3D) &
      !$omp shared(opt_restart, g1_tar, g1_3D_tar, n1_g1_tar) &
      !$omp shared(f1_tar, f1_3D_tar, n1_f1_tar) &
      !$omp shared(n0_ae, T0_ae ) &
      !$omp shared(phi1_in, phi1_out) &
      !$omp shared(itmax, phi1_abs_tol, phi1_rel_tol) &
      !$omp shared(opt_Potential, opt_Potential_converge)&
      !$omp shared(absolute_ok, relative_ok, update_ok) &
      !$omp shared(abs_error, rel_error) &
      !$omp shared(result_abs, result_rel,result_iter) &
      !$omp shared(vth, vel_volume, sml_charge_num, small_n0, eq_t_ev) &
      !$omp shared(is_s, ip_s, ith_s, ir_s, iz_s ) &
      !$omp shared(is_e, ip_e, ith_e, ir_e, iz_e ) &
      !$omp shared(tar_s, tar_e) &
      !$omp private(is, ip, ith, ir, iz) &
      !$omp private(phi_iter, numer, denom, n1_g1, phi1, phi1_old, conv_factor) 
      !$omp do 
        
        do ip=ip_s,ip_e
        do ith=ith_s,ith_e

              ! initialization
              g1_3D = g1_in(:,:,ith,ip,:) ! 초기값 이렇게 주는지는 모르겠네.
              f0_3D = f0_in(:,:,ith,ip,:)
              f1_3D = f1_in(:,:,ith,ip,:)

              g1_3D_tar = g1_tar(:,:,ith,ip,:) ! opt_restart = 3
              f1_3D_tar = f1_tar(:,:,ith,ip,:) ! opt_restart = 3

              phi1      = 0D0
              phi1_old  = phi1_in(ith,ip)

              !! Need to know 'ip'  for 'small_n0' & eq_t_ev & vel_volume
              !! Need to know 'ith' for 'small_n0'
              ! If you want to make the subroutine below precess, Need to define as follows.
              ! n0(:)        = small_n0   (ith,ip,:)
              ! T0(:)        = eq_t_ev    (ip,:)
              ! dV(:,:)      = vel_volume (:,ip,:)

              ! opt_restart == 3
              if (opt_restart .eq. 3) then
                  do is=tar_s,is_s-1
                      conv_factor = vth(ip,is)*sqrt(sml_2pi)
                      conv_factor = 1D0/conv_factor**3

                      do ir=ir_s,ir_e
                      do iz=iz_s,iz_e
                          n1_g1_tar(is) = n1_g1_tar(is) + g1_3D_tar(iz,ir,is) &
                                                        * vel_volume(ir,ip,is) &!! Need to know 'ip'
                                                        * conv_factor
                      enddo
                      enddo
                  enddo
              elseif (opt_restart .eq. 4) then
                  do is=tar_s,is_s-1
                      conv_factor = vth(ip,is)*sqrt(sml_2pi)
                      conv_factor = 1D0/conv_factor**3

                      do ir=ir_s,ir_e
                      do iz=iz_s,iz_e
                          n1_f1_tar(is) = n1_f1_tar(is) + f1_3D_tar(iz,ir,is) &
                                                        * vel_volume(ir,ip,is) &!! Need to know 'ip'
                                                        * conv_factor
                      enddo
                      enddo
                  enddo
              endif

              ! Convergence criterion is different for each phi1(ith,ip)
              do phi_iter=0,itmax

                  n1_g1 = 0D0

                  numer = 0D0
                  denom = 0D0

                  do is=is_s,is_e

                      conv_factor = vth(ip,is)*sqrt(sml_2pi)
                      conv_factor = 1D0/conv_factor**3

                      do ir=ir_s,ir_e
                      do iz=iz_s,iz_e
                          n1_g1(is) = n1_g1(is) + g1_3D(iz,ir,is) &
                                                  * vel_volume(ir,ip,is) &!! Need to know 'ip'
                                                  * conv_factor

                      enddo
                      enddo
                  enddo



                  ! write(*,*) 'phi_iter', phi_iter, 'ith', ith, 'n1_g1', n1_g1

                  !youknow memo
                  ! e*Phi/T[J] -> Phi/T[eV]
                  ! sml_charge_num : Z
                  ! sml_charge     : Z*e

                  select case (opt_Potential)
                    case(1) ! default
                        do is=is_s,is_e
                          numer = numer + sml_charge_num(is) * n1_g1(is)
                        enddo

                        ! !! opt_restart = 3
                        if (opt_restart .eq. 3) then
                          do is=tar_s,is_s-1
                            numer = numer + sml_charge_num(is) * n1_g1_tar(is)
                          enddo
                        elseif (opt_restart .eq. 4) then
                          do is=tar_s,is_s-1
                            numer = numer + sml_charge_num(is) * n1_f1_tar(is)
                          enddo
                        endif

                        if (opt_restart .eq. 4) then
                        ! when opt_restart==4, don't include 0th order term in denom.
                          do is=is_s,is_e !! is_s, is_e
                            denom = denom + sml_charge_num(is)**2 * small_n0(ith,ip,is) / eq_t_ev(ip,is) !! Need to know 'ip'
                          enddo
                        else
                          do is=tar_s,tar_e !! tar_s, tar_e
                            denom = denom + sml_charge_num(is)**2 * small_n0(ith,ip,is) / eq_t_ev(ip,is) !! Need to know 'ip'
                          enddo
                        endif


                        phi1 = numer / denom
                    case(2) ! adiabatic electron
                        do is=is_s,is_e
                          numer = numer + sml_charge_num(is) * n1_g1(is)
                        enddo
                        do is=tar_s,tar_e !! tar_s, tar_e
                          denom = denom + sml_charge_num(is)**2 * small_n0(ith,ip,is) / eq_t_ev(ip,is) !! Need to know 'ip'
                        enddo

                        denom = denom + n0_ae (ith,ip) / T0_ae (ip)

                        phi1 = numer / denom
                        
                    case default
                      write(*,*) 'error :: opt_Potential'

                  end select

                  ! Upadte g1_out for iteration
                  do is=is_s,is_e
                  do ir=ir_s,ir_e
                  do iz=iz_s,iz_e
                    ! e*Phi/T[J] -> Phi/T[eV]
                    g1_3D(iz,ir,is) = f1_3D(iz,ir,is) &
                                              + sml_charge_num(is) * phi1 / eq_t_ev(ip,is)  & !! Need to know 'ip'
                                              * f0_3D(iz,ir,is) 
                  enddo
                  enddo
                  enddo

                  ! Calcualte the error
                  select case (opt_Potential_converge)
                    case(1) ! Potential
                        abs_error = phi1 - phi1_old ! 처음에는 phi1_old가 이전 phi1인데 괜찮겠지?
                        rel_error = abs_error / phi1
                    case(2) ! Species 별로 수렴의 속도가 다를 수도 있음. <-- 어떻게 판단?
                            ! 아니면 뭐 vperp, vpar에 대해 평가한다던가

                    case default
                      write(*,*) 'error :: opt_Potential_converge'
                  end select

                  ! Update phi1_old at the end of the process, not at this step.

                  absolute_ok = abs(abs_error) .lt. phi1_abs_tol
                  relative_ok = abs(rel_error) .lt. phi1_rel_tol   ! phi1=0일 때 오류 가능.

                  ! Is it fine?
                  update_ok = ( (phi_iter.gt.0) .and.(absolute_ok .or. relative_ok) ) .or. (phi_iter.eq.itmax )

                  ! Output
                  if (update_ok) then
                    if (phi_iter.eq.itmax) then
                      write(*,*) 'potentail is not converged. but updated @ ith=', ith
                      write(*,*) 'phi_iter', phi_iter, 'n1_g1', n1_g1
                      write(*,*) 'abs_error', abs_error, 'rel_error', rel_error
                    endif

                    g1_out(:,:,ith,ip,:) = g1_3D(:,:,:)
                    phi1_out(ith,ip) = phi1

                    result_abs (ith,ip) = abs_error
                    result_rel (ith,ip) = rel_error
                    result_iter(ith,ip) = phi_iter

                    exit
                  else
                  endif

                  ! Update phi1_old
                  phi1_old = phi1

              enddo ! phi_iter


        enddo
        enddo

      
      !$omp end do
      !$acc wait
      !$omp end parallel

      !$acc wait
    
      ! post-processing?

      ! This diffrential method is different with dfdth
      call cal_dphidth(phi1_out,dphidth)

      ! write(*,*) 'result_abs', result_abs
      ! write(*,*) 'result_rel', result_rel
      ! write(*,*) 'result_iter', result_iter

    end subroutine fm_quasi

   subroutine find_potential (g1_in,f0_in, f1_out, phi1_out, dphidth)
   ! this routine to calculate the f1 & potential from g1

      use FP4D_globals, only: itmax
      use FP4D_globals, only: sml_charge_num
      use FP4D_globals, only: sml_2pi
    real (kind=8), intent(in) :: g1_in (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e)

      real (kind=8), intent(in) :: f0_in   (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e)

      real (kind=8), intent(out) :: f1_out (iz_s:iz_e, ir_s:ir_e,  &
                                           ith_s:ith_e, ip_s:ip_e, is_s:is_e)
      real(kind=8), intent(out) :: phi1_out(ith_s:ith_e,ip_s:ip_e)
      real(kind=8), intent(out) :: dphidth (ith_s:ith_e,ip_s:ip_e)

      real (kind=8) :: conv_factor

      real(kind=8) :: g1_3D (iz_s:iz_e, ir_s:ir_e, is_s:is_e)
      real(kind=8) :: f1_3D (iz_s:iz_e, ir_s:ir_e, is_s:is_e)
      real(kind=8) :: f0_3D (iz_s:iz_e, ir_s:ir_e, is_s:is_e)

      real(kind=8) :: n1_g1(is_s:is_e)
      ! real(kind=8) :: n1_f1(is_s:is_e)

      real(kind=8) :: abs_error, rel_error
      real(kind=8) :: phi1_abs_tol, phi1_rel_tol ! Need to transfer to global

      real(kind=8) :: numer, denom    ! youknow
      real(kind=8) :: phi1 , phi1_old ! youknow

      logical :: absolute_ok, relative_ok, update_ok

      integer :: alloc_stat
      integer :: is, ip, ith, ir, iz

      dphidth=0d0

      !$omp parallel default (none) num_threads(col_f_i_nthreads)&
      !$omp shared(g1_in, f0_in, f1_out) &
      !$omp shared(g1_3D, f0_3D, f1_3D ) &
      !$omp shared(n0_ae, T0_ae ) &
      !$omp shared(phi1_out) &
      !$omp shared(opt_Potential, opt_Potential_converge)&
      !$omp shared(vth, vel_volume, sml_charge_num, small_n0, eq_t_ev) &
      !$omp shared(is_s, ip_s, ith_s, ir_s, iz_s ) &
      !$omp shared(is_e, ip_e, ith_e, ir_e, iz_e ) &
      !$omp shared(tar_s, tar_e) &
      !$omp private(is, ip, ith, ir, iz) &
      !$omp private(numer, denom, n1_g1, phi1, phi1_old, conv_factor) 
      !$omp do 

      do ip=ip_s,ip_e
      do ith=ith_s,ith_e

              ! initialization
              g1_3D = g1_in(:,:,ith,ip,:) ! ?~H기?~R ?~]??| ~G?~L

              f0_3D = f0_in(:,:,ith,ip,:)
              f1_3D = 0D0
                 n1_g1 = 0D0

                  numer = 0D0
                  denom = 0D0

                  do is=is_s,is_e

                      conv_factor = vth(ip,is)*sqrt(sml_2pi)
                      conv_factor = 1D0/conv_factor**3

                      do ir=ir_s,ir_e
                      do iz=iz_s,iz_e
                          n1_g1(is) = n1_g1(is) + g1_3D(iz,ir,is) &
                                                  * vel_volume(ir,ip,is) &!!Need to know 'ip'
                                                  * conv_factor
                      enddo
                      enddo
                  enddo


                  do is=is_s,is_e
                    numer = numer + sml_charge_num(is) * n1_g1(is)
                  enddo
                  do is=tar_s,tar_e !! tar_s, tar_e
                    denom = denom + sml_charge_num(is)**2 * small_n0(ith,ip,is) / eq_t_ev(ip,is) !! Need to know 'ith,ip'
                  enddo

                  phi1 = numer / denom

          
                  do is=is_s,is_e
                  do ir=ir_s,ir_e
                  do iz=iz_s,iz_e
                     f1_3D(iz,ir,is) = g1_3D(iz,ir,is) &
                                   - sml_charge_num(is) /eq_t_ev(ip,is)  & !!Need to know 'ip'
                                              * phi1 &
                                              * f0_3D(iz,ir,is)
                  enddo
                  enddo
                  enddo
                  f1_out(:,:,ith,ip,:) = f1_3D(:,:,:)
                  phi1_out(ith,ip) = phi1


      enddo
      enddo

      !$omp end do
      !$acc wait
      !$omp end parallel

      call cal_dphidth(phi1_out,dphidth)

  end subroutine find_potential     
end module quasi_f_module

