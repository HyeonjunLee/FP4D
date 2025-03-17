#ifdef USE_ASYNC
#define ASYNC(istream)  async(istream)
#else
#define ASYNC(istream)
#endif

module RHS_module
  use FP4D_globals
  use FP4D_math

  integer :: rhs_i_nthreads                     !! [rhs_setup] the number of inner OMP threads
  integer :: rhs_o_nthreads                     !! [rhs_setup] the number of outer OMP threads

  contains
    ! BLAS solver collision operator time-stepping
#include "bsolver.F90"

    !> Set up global variables of the collision operator
    !! @param[in]   mype   MPI rank ID
    !!
    subroutine RHS_f_setup(mype)
      implicit none
      integer, intent(in) ::  mype
    
      rhs_i_nthreads = col_f_nthreads
      call set_omp_outer_thread_num( rhs_o_nthreads )
      if(mype .eq. 0) then
          write(*,'(a,i4)') '[MULTI-PRL] RHS_module: inner thread = ', rhs_i_nthreads
          write(*,'(a,i4)') '[MULTI-PRL] RHS_module: outer thread = ', rhs_o_nthreads
      endif

      return

    end subroutine RHS_f_setup

    !> compute the number of threads used in the outer
    !! level OpenMP parallelization.
    !! The number of inner threads is set with an input
    !! parameter (rhs_nthreads)
    !! --> number of outer threads = OMP_NUM_THREADS/#(inner threads)
    !!
    subroutine set_omp_outer_thread_num( nthreads )
      use FP4D_globals, only : sml_mype
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



    subroutine cal_PS(f0_in,f1_in,g1_in,pdf_PS)
      use FP4D_globals, only: sml_pi, sml_e_charge
      use FP4D_globals, only: opt_PS, opt_PS_aux

      use FP4D_globals, only: sml_dt

      use FP4D_globals, only: sml_2pi

      use FP4D_globals, only: opt_time_solver, opt_time_species, opt_rotation

#ifdef _OPENMP
      use omp_lib
#endif
      !use mpi
      implicit none

      real (kind=8), intent(in) :: f0_in(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)
      real (kind=8), intent(in) :: f1_in(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)
      real (kind=8), intent(in) :: g1_in (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)

      real (kind=8), intent(out) :: pdf_PS(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)
                                         
      integer :: alloc_stat
      integer :: isize

      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dfdrho ! Taeyeon?
      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dfdth , d2fdth2
      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dfdr  , d2fdr2
      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dfdz  , d2fdz2

      real(kind=8), dimension(ith_s:ith_e, ip_s:ip_e) :: local_dPhi1dth, local_d2Phi1dth2

      integer :: stat
      
      real (kind=8) :: PS_term1, PS_term_rot
      real (kind=8) :: PS_term2, PS_term3, PS_term4, common_PS
      integer :: is, ip, ith, ir, iz

      pdf_PS = 0D0

      ! write(*,*) 'norm_dr', norm_dr
      ! write(*,*) 'norm_dz', norm_dz
      ! write(*,*) 'eq_dtheta', eq_dtheta

      select case (opt_PS)
        case (1)
          !-----------------------------------------------------------------------
          ! START DIFFERENTIATION for PS
          !   - Using collapse(4). 
          !   - Using only collapse(3) results in low thread utilization.
          !     - Because Nsp(=2)*Npsi(=1)*Nth(=16) = 32 for basic case.
          !
          dfdr  = 0D0; d2fdr2  = 0D0
          dfdz  = 0D0; d2fdz2  = 0D0
          dfdth = 0D0; d2fdth2 = 0D0

          ! dPhi1dth is calculated @ quasi_f_module.F90
          ! local_dPhi1dth   = 0D0
          ! local_d2Phi1dth2 = 0D0
          ! do ip = ip_s, ip_e
          !   call cal_dfdth6  (eq_phi1(:,ip), local_dPhi1dth  (:,ip))
          !   call cal_dfdths  (eq_phi1(:,ip), local_d2Phi1dth2(:,ip))
          ! enddo
          
          !$omp parallel do collapse(4) private(is, ip, ith, ir, iz) num_threads(rhs_i_nthreads)
          do is = is_s, is_e
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
          do ir = ir_s, ir_e
            call cal_dfdz8   (g1_in(:,ir,ith,ip,is), dfdz  (:,ir,ith,ip,is))
            ! call cal_d2fdz28 (g1_in(:,ir,ith,ip,is), d2fdz2(:,ir,ith,ip,is))
          enddo ! ir
          enddo ! ith
          enddo ! ip
          enddo ! is
          !$omp end parallel do

          !$omp parallel do collapse(4) private(is, ip, ith, ir, iz) num_threads(rhs_i_nthreads)
          do is = is_s, is_e
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
          do iz = iz_s, iz_e
            call cal_dfdr8   (g1_in(iz,:,ith,ip,is), dfdr  (iz,:,ith,ip,is))
            ! call cal_d2fdr28 (g1_in(iz,:,ith,ip,is), d2fdr2(iz,:,ith,ip,is))
          enddo ! iz
          enddo ! ith
          enddo ! ip
          enddo ! is
          !$omp end parallel do

          !$omp parallel do collapse(4) private(is, ip, ith, ir, iz) num_threads(rhs_i_nthreads)
          do is = is_s, is_e
          do ip = ip_s, ip_e
          do ir = ir_s, ir_e
          do iz = iz_s, iz_e
              call cal_dfdth6  (g1_in(iz,ir,:,ip,is), dfdth  (iz,ir,:,ip,is))
              call cal_dfdths  (g1_in(iz,ir,:,ip,is), d2fdth2(iz,ir,:,ip,is))
          enddo ! iz
          enddo ! ir
          enddo ! ip
          enddo ! is
          !$omp end parallel do

          ! !$omp parallel do collapse(4) private(is, ip, ith, ir, iz) num_threads(rhs_i_nthreads)
          ! do is = is_s, is_e
          ! do ith = ith_s, ith_e
          ! do ir = ir_s, ir_e
          ! do iz = iz_s, iz_e
          !     ! call cal_dfdrho  (g1_in(iz,ir,ith,:,is), dfdpsi  (iz,ir,ith,:,is))
          ! enddo ! iz
          ! enddo ! ith
          ! enddo ! ir
          ! enddo ! is
          ! !$omp end parallel do

          ! write(*,*) 'dfdr', dfdr(:,:,ith_s,ip_s,is_s)
          ! write(*,*) 'dfdz', dfdz(:,:,ith_s,ip_s,is_s)
          ! write(*,*) 'dfdth', dfdth(:,:,ith_s,ip_s,is_s)
          ! write(*,*) 'd2fdth2', d2fdth2(:,:,ith_s,ip_s,is_s)

          !
          ! END DIFFERENTIATION for PS
          !-----------------------------------------------------------------------

          !-----------------------------------------------------------------------
          ! START main process of PS
          !   - Using collapse(4) to preventing cashe overhead
          !

          !$omp parallel default(none) num_threads(rhs_i_nthreads) &
          !$omp shared(f0_in, f1_in  , g1_in, pdf_PS   ) &
          !$omp shared(dfdrho  , dfdth, d2fdth2, dfdr , dfdz     )&
          !$omp shared(is_s, ip_s, ith_s, ir_s, iz_s ) &
          !$omp shared(is_e, ip_e, ith_e, ir_e, iz_e ) &
          !$omp shared(sml_dt, eq_tau,&
          !$omp        norm_vz, norm_vr, vth, &
          !$omp        eq_dBdth, eq_B, eq_dphidth, eq_omega, eq_invJ, &
          !$omp        opt_PS_aux) &
          !$omp shared (opt_rotation, dlambda_rot_dth) &
          !$omp shared(opt_time_solver,opt_time_species) &
          !$omp private(is, ip, ith, ir, iz, &
          !$omp         PS_term1, PS_term_rot, &
          !$omp         PS_term2, PS_term3, PS_term4, common_PS)
          !$omp do collapse(4)

                ! collapse(4) enables the parallelization of (is, ip, ith, ir) with num_threads(rhs_i_nthreads).

            !! START PS
            do is = is_s, is_e
            do ip = ip_s, ip_e
            do ith = ith_s, ith_e
            do ir = ir_s, ir_e
            do iz = iz_s, iz_e
        
                PS_term1 = norm_vz(iz)*dfdr(iz,ir,ith,ip,is) - norm_vr(ir)*dfdz(iz,ir,ith,ip,is)
                PS_term1 = PS_term1 * 0.5 * norm_vr(ir) * eq_dBdth(ith,ip) / eq_B(ith,ip)

                if (opt_rotation.ne.0) then
                  PS_term_rot = - dlambda_rot_dth(ith,ip,is) * dfdz(iz,ir,ith,ip,is)
                  PS_term1 = PS_term1 + PS_term_rot
                endif

                PS_term2 = norm_vz(iz)*dfdth(iz,ir,ith,ip,is)

                PS_term3 = abs(norm_vz(iz))*d2fdth2(iz,ir,ith,ip,is)

                PS_term4 = - eq_dphidth(ith,ip) * eq_omega(ith,ip,is)/eq_B(ith,ip) * dfdz(iz,ir,ith,ip,is) / vth(ip,is)**2.0 ! 1/vth**2 is important@!!
                ! PS_term4 = - local_dPhi1dth(ith,ip) * eq_omega(ith,ip,is)/eq_B(ith,ip) * dfdz(iz,ir,ith,ip,is) / vth(ip,is)**2.0 ! 1/vth**2 is important@!!



                common_PS = -sml_dt *vth(ip,is) &         ! (-) sign due to PS from LHS
                            *eq_invJ(ith,ip)/eq_B(ith,ip)
                if (opt_time_solver.eq.-1) then ! steady
                  common_PS = common_PS * eq_tau(ip,opt_time_species,is) ! eq_tau(psi, target, colliding)
                endif
                ! if ( (is.eq.is_s) .and. (ip.eq.ip_s) .and. (ith.eq.ith_s)) then
                !   write(*,*) 'PS_term1', PS_term1
                !   write(*,*) 'PS_term2', PS_term2
                !   write(*,*) 'PS_term3', PS_term3
                !   write(*,*) 'PS_term4', PS_term4
                !   write(*,*) 'common_PS', common_PS
                ! endif

                select case (opt_PS_aux)
                  case (1) ! default : 1
                    pdf_PS(iz,ir,ith,ip,is) = common_PS * (PS_term1 + PS_term2 + PS_term3 + PS_term4)
                                            
                  case (100001)  ! Nth=1 case for NB ==> No dfdth
                    pdf_PS(iz,ir,ith,ip,is) = common_PS * (PS_term1 )

                  case (100002)  ! Nth=48 case for NB ==> dfdth O dissipation O but No dphidth
                    pdf_PS(iz,ir,ith,ip,is) = common_PS * (PS_term1 + PS_term2 + PS_term3 )

                  case (100003)  ! Nth=48 case for NB ==> dfdth O dissipation O but No dphidth No ..
                    pdf_PS(iz,ir,ith,ip,is) = common_PS * (PS_term2 +PS_term3)

                  case (100004)  ! Nth=48 case for NB ==> dfdth O dissipation O but No dphidth No ..
                    pdf_PS(iz,ir,ith,ip,is) = common_PS * (PS_term2)



                  case default
                      write(*,*) 'error :: opt_PS_aux'
                end select ! case (opt_PS_aux)
                  
            enddo ! iz
            enddo ! ir
            enddo ! ith
            enddo ! ip
            enddo ! is

            ! end PS term

                
            !$omp end do
            !$omp end parallel

          !$acc wait
          
          !
          ! END main process of PS
          !-----------------------------------------------------------------------

        case(0) ! case (opt_PS)
          pdf_PS = 0D0
        case default! case (opt_PS)
          write(*,*) 'error :: opt_PS'
      end select ! case (opt_PS)

    end subroutine cal_PS


    subroutine cal_RF(f0_in,f1_in,pdf_RF)
      use FP4D_globals, only: sml_pi, sml_e_charge
      use FP4D_globals, only: opt_QLD, opt_QLD_aux, opt_QLD_dist

      use FP4D_globals, only: sml_dt
      use FP4D_globals, only: sml_2pi

      use FP4D_globals, only: RF_n_harmonic, RF_omega, RF_kpar, RF_W0, RF_width

      use FP4D_globals, only: opt_time_solver, opt_time_species

      use f_module, only: f1_an, f1_neo

#ifdef _OPENMP
      use omp_lib
#endif
      !use mpi
      implicit none

      real (kind=8), intent(in) :: f0_in(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)
      real (kind=8), intent(in) :: f1_in(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)                                   

      real (kind=8), intent(out) :: pdf_RF(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e,&
                                         ip_s:ip_e, is_s:is_e)

      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dist_in
      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dfdr  , d2fdr2
      real(kind=8), dimension(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e) :: dfdz  , d2fdz2

      integer :: alloc_stat
      integer :: isize

      integer :: stat
      
      real (kind=8) :: RF_term_diff, common_RF
      real (kind=8) :: RF_term2_upper, RF_term2_lower, RF_term2_zero
      real (kind=8) :: RF_term3
      real (kind=8) :: RF_term23, RF_term23_upper, RF_term23_lower
      real (kind=8) :: RF_term4

      integer :: is, ip, ith, ir, iz

      logical :: cond_Karney

      pdf_RF = 0D0

      ! write(*,*) 'norm_dr', norm_dr
      ! write(*,*) 'norm_dz', norm_dz
      ! write(*,*) 'eq_dtheta', eq_dtheta


      select case (opt_QLD)
        case (1)
          !-----------------------------------------------------------------------
          ! START DIFFERENTIATION for QLD
          !   - Using collapse(4). 
          !   - Using only collapse(3) results in low thread utilization.
          !     - Because Nsp(=2)*Npsi(=1)*Nth(=16) = 32 for basic case.
          !
          dfdr  = 0D0; d2fdr2  = 0D0
          dfdz  = 0D0; d2fdz2  = 0D0

          select case (opt_QLD_dist)
            case(0)   ! default :: analytic differentiattion of f0
                ! Nothing
            case(1,2,3,4) ! numerical differentiattion

                if (opt_QLD_dist .eq. 1) then     ! f0
                  dist_in = f0_in
                elseif (opt_QLD_dist .eq. 2) then ! f0+f1
                  dist_in = f0_in + f1_in
                elseif (opt_QLD_dist .eq. 3) then ! f1_an [Helander Book] eq(8.16)
                  dist_in = f1_an
                elseif (opt_QLD_dist .eq. 4) then ! f1_an :: using restart f1 (considering collision + drift)
                  dist_in = f1_neo
                endif

                ! write(*,*) 'dist_in', dist_in

               !$omp parallel do collapse(4) private(is, ip, ith, ir) num_threads(rhs_i_nthreads)
                do is = is_s, is_e
                do ip = ip_s, ip_e
                do ith = ith_s, ith_e
                do ir = ir_s, ir_e
                  call cal_dfdz8   (dist_in(:,ir,ith,ip,is), dfdz  (:,ir,ith,ip,is))
                  call cal_d2fdz28 (dist_in(:,ir,ith,ip,is), d2fdz2(:,ir,ith,ip,is))
                enddo ! ir
                enddo ! ith
                enddo ! ip
                enddo ! is
               !$omp end parallel do

               !$omp parallel do collapse(4) private(is, ip, ith, iz) num_threads(rhs_i_nthreads)
                do is = is_s, is_e
                do ip = ip_s, ip_e
                do ith = ith_s, ith_e
                do iz = iz_s, iz_e
                  call cal_dfdr8   (dist_in(iz,:,ith,ip,is), dfdr  (iz,:,ith,ip,is))
                  call cal_d2fdr28 (dist_in(iz,:,ith,ip,is), d2fdr2(iz,:,ith,ip,is))
                enddo ! iz
                enddo ! ith
                enddo ! ip
                enddo ! is
               !$omp end parallel do

                ! write(*,*) 'dfdr', dfdr(:,:,ith_s,ip_s,is_s)
                ! write(*,*) 'dfdz', dfdz(:,:,ith_s,ip_s,is_s)
                ! write(*,*) 'd2fdr2', d2fdr2(:,:,ith_s,ip_s,is_s)
                ! write(*,*) 'd2fdz2', d2fdz2(:,:,ith_s,ip_s,is_s)

            case default
                write(*,*) 'error :: opt_QLD_dist'
          end select ! case(opt_QLD_dist)

          ! END DIFFERENTIATION for QLD
          !-----------------------------------------------------------------------


          !-----------------------------------------------------------------------
          ! START main process of QLD
          !   - Using collapse(4) to preventing cashe overhead
          !

          !$omp parallel default(none) num_threads(rhs_i_nthreads) &
          !$omp shared(f0_in, f1_in, pdf_RF   ) &
          !$omp shared( dfdr,d2fdr2, dfdz,d2fdz2    )&
          !$omp shared(is_s, ip_s, ith_s, ir_s, iz_s ) &
          !$omp shared(is_e, ip_e, ith_e, ir_e, iz_e ) &
          !$omp shared(sml_dt, eq_tau, nth, &
          !$omp        norm_vz, norm_vr, norm_vr0, vth, &
          !$omp        eq_omega,eq_dtheta, &
          !$omp        opt_QLD_dist, opt_QLD_aux, &
          !$omp        RF_n_harmonic, RF_kpar,RF_W0,RF_omega,RF_width) &
          !$omp shared(RF_vr1,RF_vr2,RF_vz1,RF_vz2,cond_Karney) &
          !$omp shared(f1_an,f1_neo) &
          !$omp shared(opt_time_solver,opt_time_species) &
          !$omp private(is, ip, ith, ir, iz, &
          !$omp          RF_term_diff, common_RF, &
          !$omp          RF_term2_upper, RF_term2_lower, RF_term2_zero, &
          !$omp          RF_term3, &
          !$omp          RF_term23, RF_term23_upper, RF_term23_lower, &
          !$omp          RF_term4 )
          !$omp do collapse(4)
                ! collapse(4) enables the parallelization of (is, ip, ith, ir) with num_threads(rhs_i_nthreads).

            !! START QLD
            do is = is_s, is_e
            do ip = ip_s, ip_e
            do ith = ith_s, ith_e
            do ir = ir_s, ir_e
            do iz = iz_s, iz_e
              
                select case (opt_QLD_dist)
                  case(0)
                    ! Analytic differentiattion of f0
                    RF_term_diff = (norm_vr(ir)**2.0-2.0*RF_n_harmonic(is)) &
                            *f0_in(iz,ir,ith,ip,is) 
                  case(1,2,3,4)
                    ! Numerical differentiattion of f0 | f0+f1 | f1_an
                    RF_term_diff = (2.0*RF_n_harmonic(is)-1.0)*dfdr(iz,ir,ith,ip,is)/norm_vr(ir) +d2fdr2(iz,ir,ith,ip,is)
                end select

                ! dirac-delta for theta
                RF_term2_upper = delta_rectangular(real(ith,kind=8),real(nth,kind=8)*1.0/4.0,eq_dtheta) ! dirac-delta(th-pi/2)
                RF_term2_lower = delta_rectangular(real(ith,kind=8),real(nth,kind=8)*3.0/4.0,eq_dtheta) ! dirac-delta(th+pi/2)
                ! RF_term2_zero  = delta_rectangular(ith,0      ,eq_dtheta) ! dirac-delta(th+0   )

                ! dirac-delta for vpar
                RF_term3 = delta_gaussian(norm_vz(iz), RF_kpar(is), RF_width(is))

                ! If RF_kpar = 0 --> error Nan
                RF_term23 = RF_omega(is) - RF_n_harmonic(is)*eq_omega(ith,ip,is)
                RF_term23 = RF_term23/vth(ip,is)/RF_kpar(is)
                RF_term23 = delta_gaussian(norm_vz(iz), RF_term23, RF_width(is))
                RF_term23 = RF_term23 / abs(RF_kpar(is)) / vth(ip,is)

                ! ith=0,nth-1 ! ex)Nth=16, upper=0~7, lower=8~15
                if (ith .lt. nth/2) then ! [0,pi)
                  RF_term23_upper = RF_term23
                  RF_term23_lower = 0D0
                elseif (ith .lt. nth/2) then ! [pi,2pi)
                  RF_term23_upper = 0D0
                  RF_term23_lower = RF_term23
                endif

                ! Karney-1981
                cond_Karney = (norm_vr0(ir) .ge. RF_vr1) .and. &
                                        (norm_vr0(ir) .le. RF_vr2) .and. &
                                        (norm_vz (iz) .ge. RF_vz1) .and. &
                                        (norm_vz (iz) .le. RF_vz2)

                if (cond_Karney) then
                  RF_term4 = 1D0
                else
                  RF_term4 = 0D0
                endif



                common_RF = +sml_dt * RF_W0(is) &
                        *norm_vr(ir)**(2.0*RF_n_harmonic(is)-2.0) &
                        / 2.0**(RF_n_harmonic(is)-1.0) &
                        /vth(ip,is)**2.0

                if (opt_time_solver.eq.-1) then ! steady
                  common_RF = common_RF * eq_tau(ip,opt_time_species,is) ! eq_tau(psi, target, colliding)
                endif

                ! if ( (is.eq.is_s) .and. (ip.eq.ip_s) .and. (ith.eq.ith_s)) then
                !   write(*,*) 'RF_term_diff', RF_term_diff
                !   write(*,*) 'RF_term3', RF_term3
                !   write(*,*) 'term4', term4
                !   write(*,*) 'common_RF', common_RF
                ! endif

                select case (opt_QLD_aux)
                  ! 1x : Q(f(+th)) :: upper mid-plane
                  ! 2x : Q(f(-th)) :: lower mid-plane
                  ! 3x : Q+ = [ Q(f(+th)) + Q(f(-th)) ] /2
                  ! 4x : Q- = [ Q(f(+th)) - Q(f(-th)) ] /2

                  ! x0 : delta(th)
                  ! x1 : delta(th) * delta(vpar) --> no-fundamental youknow form .. lol
                  ! x2 : delta(omega_RF-kpar*vpar-n_harmonic*Omega) --> general form
                  ! x3 : delta(th) * Kareny-1981

                  ! 1000 : same power Q for every theta w/o delta(vpar)
                  ! 1001 : same power Q for every theta w/  delta(vpar)

                  case (1000)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff
                  case (1001)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * RF_term3

                  case (10)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+1.0*RF_term2_upper)
                  case (11)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+1.0*RF_term2_upper) * RF_term3
                  case (12)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+1.0*RF_term23_upper)

                  case (20)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+1.0*RF_term2_lower)
                  case (21)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+1.0*RF_term2_lower) * RF_term3
                  case (22)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+1.0*RF_term23_lower)

                  case (30) ! Q+
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term2_upper+0.5*RF_term2_lower)
                  case (31)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term2_upper+0.5*RF_term2_lower) * RF_term3
                  case (32)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term23_upper+0.5*RF_term23_lower)
                  case (33)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term23_upper+0.5*RF_term23_lower) * RF_term4


                  case (40) ! Q- -- we need this term
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term2_upper-0.5*RF_term2_lower)
                  case (41)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term2_upper-0.5*RF_term2_lower) * RF_term3
                  case (42)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (+0.5*RF_term23_upper-0.5*RF_term23_lower)
 

                  case (50) 
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (-0.5*RF_term2_upper-0.5*RF_term2_lower)
                  case (51)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (-0.5*RF_term2_upper-0.5*RF_term2_lower) * RF_term3
                  case (52)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (-0.5*RF_term23_upper-0.5*RF_term23_lower)

                  case (60) 
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (-0.5*RF_term2_upper+0.5*RF_term2_lower)
                  case (61)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (-0.5*RF_term2_upper+0.5*RF_term2_lower) * RF_term3
                  case (62)
                      pdf_RF(iz,ir,ith,ip,is) = common_RF * RF_term_diff * (-0.5*RF_term23_upper+0.5*RF_term23_lower)

                  case default
                    write (*,*) 'error :: opt_QLD_aux'
                end select

            enddo ! iz
            enddo ! ir
            enddo ! ith
            enddo ! ip
            enddo ! is

            ! end RF term

            !$omp end do
            !$omp end parallel

          !$acc wait

          !
          ! END main process of QLD
          !-----------------------------------------------------------------------

        case(0,2) ! case (opt_QLD)
          pdf_RF = 0D0
        case default! case (opt_QLD)
          write(*,*) 'error :: opt_QLD'
      end select ! case (opt_QLD)

    end subroutine cal_RF



end module RHS_module
