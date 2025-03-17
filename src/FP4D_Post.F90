module FP4D_Post

  use FP4D_globals
  use FP4D_math

  integer :: post_i_nthreads
  integer :: post_o_nthreads

  contains


    subroutine Post_f_setup(mype)
      implicit none
      integer, intent(in) ::  mype
    
      post_i_nthreads = col_f_nthreads
      call set_omp_outer_thread_num( post_o_nthreads )
      if(mype .eq. 0) then
          write(*,'(a,i4)') '[MULTI-PRL] inner thread = ', post_i_nthreads
          write(*,'(a,i4)') '[MULTI-PRL] outer thread = ', post_o_nthreads
      endif

      return

    end subroutine Post_f_setup

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

  subroutine vel_integral()
    implicit none
    integer :: is, ip, ith, ir, iz
    real (kind=8) :: conv_factor, element, element_col

    call set_omp_outer_thread_num( post_o_nthreads )

    allocate(pflux    (ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(eflux    (ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(mflux    (ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(pflux_col(ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(eflux_col(ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(mflux_col(ith_s:ith_e, ip_s:ip_e, is_s:is_e))

    allocate(FSA_pflux(ip_s:ip_e, is_s:is_e))
    allocate(FSA_eflux(ip_s:ip_e, is_s:is_e))
    allocate(FSA_mflux(ip_s:ip_e, is_s:is_e))
    allocate(FSA_pflux_col(ip_s:ip_e, is_s:is_e))
    allocate(FSA_eflux_col(ip_s:ip_e, is_s:is_e))
    allocate(FSA_mflux_col(ip_s:ip_e, is_s:is_e))

    allocate(post_dens(ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(post_NUpar(ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(post_NUparB(ith_s:ith_e, ip_s:ip_e, is_s:is_e))

    allocate(FSA_post_dens(ip_s:ip_e, is_s:is_e))
    allocate(FSA_post_NUpar(ip_s:ip_e, is_s:is_e))
    allocate(FSA_post_NUparB(ip_s:ip_e, is_s:is_e))

    allocate(jpar     (ith_s:ith_e, ip_s:ip_e))
    allocate(jparB    (ith_s:ith_e, ip_s:ip_e))
    allocate(FSA_jpar (ip_s:ip_e))
    allocate(FSA_jparB(ip_s:ip_e))


    pflux = 0D0; eflux = 0D0; mflux = 0D0
    FSA_pflux = 0D0; FSA_eflux = 0D0; FSA_mflux = 0D0;

    pflux_col = 0D0; eflux_col = 0D0; mflux_col = 0D0;
    FSA_pflux_col = 0D0; FSA_eflux_col = 0D0; FSA_mflux_col = 0D0;

    post_dens = 0D0; post_NUpar = 0D0; post_NUparB = 0D0
    FSA_post_dens = 0D0; FSA_post_NUpar = 0D0; FSA_post_NUparB = 0D0

          !$omp parallel default(none) num_threads(post_o_nthreads) &
          !$omp shared(is_s, ip_s, ith_s, ir_s, iz_s ) &
          !$omp shared(is_e, ip_e, ith_e, ir_e, iz_e ) &
          !$omp shared(vth, vel_volume, norm_vz, norm_vr, norm_vr0  ) &
          !$omp shared(eq_B ,eq_J,eq_intJ,sml_charge_num,eq_dtheta ) &
          !$omp shared(drift_psi  ) &
          !$omp shared(sml_dt, step_Col, eq_tau, eq_I, eq_omega, sml_mass) &
          !$omp shared(g1, Ff, dfdt_c) &
          !$omp shared(pflux, eflux, mflux, pflux_col, eflux_col, mflux_col) &
          !$omp shared(FSA_pflux, FSA_eflux, FSA_mflux, FSA_pflux_col, FSA_eflux_col, FSA_mflux_col) &
          !$omp shared(post_dens, post_NUpar, post_NUparB) &
          !$omp shared(FSA_post_dens, FSA_post_NUpar, FSA_post_NUparB) &
          !$omp private(is, ip, ith, ir, iz, &
          !$omp         conv_factor, element, element_col)
          !$omp do

            do is = is_s, is_e
            do ip = ip_s, ip_e
              do ith = ith_s, ith_e
                  do ir = ir_s, ir_e
                  do iz = iz_s, iz_e

                      ! conv_factor
                      conv_factor = vth(ip,is)*sqrt(6.2831853071795862D0)
                      conv_factor = 1D0/conv_factor**3

                      element     = conv_factor * vel_volume(ir,ip,is) * g1(iz,ir,ith,ip,is)
                      element_col = conv_factor * vel_volume(ir,ip,is) * dfdt_c(iz,ir,ith,ip,is)

                      ! pflux for psi
                      pflux(ith,ip,is) = pflux(ith,ip,is)          &
                                          + element                &
                                          * drift_psi  (ith,ip,is) &
                                          * vth(ip,is)*vth(ip,is)  &
                                          * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir)/2D0)

                      ! eflux for psi
                      eflux(ith,ip,is) = eflux(ith,ip,is)          &
                                          + element                &
                                          * drift_psi  (ith,ip,is) &
                                          * vth(ip,is)*vth(ip,is)  &
                                          * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir)/2D0) &
                                          * 0.5 * sml_mass(is) * vth(ip,is)*vth(ip,is)  &
                                          * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir))

                      ! mflux for psi
                      mflux(ith,ip,is) = mflux(ith,ip,is)          &
                                          + element                &
                                          * drift_psi  (ith,ip,is) &
                                          * vth(ip,is)*vth(ip,is)  &
                                          * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir)/2D0) &
                                          * sml_mass(is) * vth(ip,is) * norm_vz(iz) &
                                          * eq_I(ith,ip) / eq_B(ith,ip) 


                      ! pflux for psi using collision
                      pflux_col(ith,ip,is) = pflux_col(ith,ip,is)  &
                                          - element_col            &
                                          / sml_dt / step_Col       &
                                          / eq_tau(ip,1,is)        & !! need to modify
                                          * (vth(ip,is) * norm_vz(iz)) &
                                          * eq_I(ith,ip) / eq_omega(ith,ip,is)
                                          
                      ! eflux for psi using collision
                      eflux_col(ith,ip,is) = eflux_col(ith,ip,is)  &
                                          - element_col            &
                                          / sml_dt / step_Col       &
                                          / eq_tau(ip,1,is)        & !! need to modify
                                          * (vth(ip,is) * norm_vz(iz)) &
                                          * eq_I(ith,ip) / eq_omega(ith,ip,is) &
                                          * 0.5 * sml_mass(is) * vth(ip,is)*vth(ip,is)  &
                                          * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir))

                      ! mflux for psi using collision
                      mflux_col(ith,ip,is) = mflux_col(ith,ip,is)  &
                                          - element_col            &
                                          / sml_dt / step_Col       &
                                          / eq_tau(ip,1,is)        & !! need to modify
                                          * (vth(ip,is) * norm_vz(iz)) &
                                          * eq_I(ith,ip) / eq_omega(ith,ip,is) &
                                          * sml_mass(is) * vth(ip,is) * norm_vz(iz) &
                                          * eq_I(ith,ip) / eq_B(ith,ip) 

                      ! post_***
                      post_dens(ith,ip,is) = post_dens(ith,ip,is)          &
                                          + element
                      post_NUpar(ith,ip,is) = post_NUpar(ith,ip,is)          &
                                          + element                &
                                          * (vth(ip,is) * norm_vz(iz))
                      post_NUparB(ith,ip,is) = post_NUparB(ith,ip,is)          &
                                          + element                &
                                          * (vth(ip,is) * norm_vz(iz)) * eq_B(ith,ip)



     
                  enddo ! iz
                  enddo ! ir

                  ! Caution :: integeral (nv) = NV
                  ! post_Upar(ith,ip,is) = post_NUpar(ith,ip,is)/post_dens(ith,ip,is)
                  ! post_UparB(ith,ip,is) = post_NUparB(ith,ip,is)/post_dens(ith,ip,is)


              enddo ! ith

              ! FSA (iterate for ip because iteration of ith conducted in compute_FSA)
              FSA_pflux(ip,is) = compute_FSA(pflux(:,ip,is),ip)
              FSA_eflux(ip,is) = compute_FSA(eflux(:,ip,is),ip)
              FSA_mflux(ip,is) = compute_FSA(mflux(:,ip,is),ip)

              FSA_pflux_col(ip,is) = compute_FSA(pflux_col(:,ip,is),ip)
              FSA_eflux_col(ip,is) = compute_FSA(eflux_col(:,ip,is),ip)
              FSA_mflux_col(ip,is) = compute_FSA(mflux_col(:,ip,is),ip)

              FSA_post_dens(ip,is) = compute_FSA(post_dens(:,ip,is),ip)
              FSA_post_NUpar(ip,is) = compute_FSA(post_NUpar(:,ip,is),ip)
              FSA_post_NUparB(ip,is) = compute_FSA(post_NUparB(:,ip,is),ip)

            enddo ! ip
            enddo ! is
          !$omp end do
          !$acc wait
          !$omp end parallel

          ! jpar(charge summation) must be calculated after intergration onf each species
          jpar  = 0D0; FSA_jpar  = 0D0;
          jparB = 0D0; FSA_jparB = 0D0

          do ip = ip_s, ip_e
            do ith = ith_s, ith_e
              do is = is_s, is_e
                jpar (ith,ip) = jpar (ith,ip) + post_NUpar (ith,ip,is) * sml_charge_num(is) * sml_e_charge
                jparB(ith,ip) = jparB(ith,ip) + post_NUparB(ith,ip,is) * sml_charge_num(is) * sml_e_charge
              enddo
            enddo
            FSA_jpar (ip) = compute_FSA(jpar (:,ip),ip)
            FSA_jparB(ip) = compute_FSA(jparB(:,ip),ip)
          enddo



  end subroutine vel_integral

  subroutine Post_record()
    implicit none
    integer :: is, ip, ith, ir, iz
    real :: sum_pflux(ip_s:ip_e)
    real :: inv_dpsidr(ip_s:ip_e)
    real(kind=8) :: vth_tmp

    ! print*,'----------------------------------------------------------------------'
    ! print*,'                               Post_record'
    ! print*,'----------------------------------------------------------------------'
    
    ! write(*,*) 'FSA_pflux for psi', FSA_pflux
    ! write(*,*) 'FSA_post_dens', FSA_post_dens
    ! write(*,*) 'FSA_post_NUpar', FSA_post_NUpar
    ! write(*,*) 'FSA_post_NUparB', FSA_post_NUparB

    ! write(*,*) 'FSA_jpar', FSA_jpar
    ! write(*,*) 'FSA_jparB', FSA_jparB


    ! ! refer NEO
    ! open(unit=io_FP4Dout, file=runfile_FP4Dout, status='replace')
    ! close(io_FP4Dout)

    open(unit=io_FP4Dout, file=runfile_FP4Dout, status='old', position='append')     
      
      sum_pflux=0.0
      do ip=ip_s,ip_e
        do is=is_s, is_e
          sum_pflux(ip) = sum_pflux(ip) + sml_charge_num(is) * FSA_pflux(ip,is)
        enddo
      enddo
      
      write(io_FP4Dout,*)

      write(io_FP4Dout,*) 'Caution :: eflux, mflux are not calcualted in present'
      write(io_FP4Dout,*) ''      

      write(io_FP4Dout,*) 'rho : sqrt(normalized pololidal flux)'
      write(io_FP4Dout,*) ''      

      do ip=ip_s,ip_e

        ! e14.5 type : 
        write(io_FP4Dout,'(t2,a,t21,e14.5)') 'rho = ', eq_psin(ip)
        ! es14.5 type
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'jpar = ', FSA_jpar(ip)
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'jparB = ', FSA_jparB(ip)
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'sum Z_s Gamma_s = ', sum_pflux(ip)
        write(io_FP4Dout,'(t3,a,t11,a,t25,a,t39,a)') &
            'Z', 'pflux', 'eflux', 'mflux'

        do is=is_s, is_e
          write (io_FP4Dout,'(f7.3,3(es14.5))') &
                sml_charge_num(is),FSA_pflux(ip,is),FSA_eflux(ip,is),FSA_mflux(ip,is)
        enddo

      enddo

      write(io_FP4Dout,*) '------------------'
      write(io_FP4Dout,*)

    close(io_FP4Dout)


!!!!!!! NEO unit

    ! considering :: inv_dpsidr = np.average(1.0/data['Bp']/data['geo_R'].T)

    open(unit=io_FP4Dout, file=runfile_FP4Dout, status='old', position='append')     
      write(io_FP4Dout,*) "-------------------------------"
      write(io_FP4Dout,*) 'To compare with NEO FLUX'
      write(io_FP4Dout,*) '  from psi to r'
      write(io_FP4Dout,*) '  normalizatino : refer to out.fp4d.neo_input'
      write(io_FP4Dout,*) "-------------------------------"

      inv_dpsidr = 0D0
      do ip=ip_s,ip_e
        do ith=ith_s,ith_e
          inv_dpsidr(ip) = inv_dpsidr(ip) + 1.0 / eq_Bp(ith,ip)/geo_R(ip,ith)
        enddo
        inv_dpsidr(ip) = inv_dpsidr(ip) / nth
      enddo

      do ip=ip_s,ip_e

        vth_tmp = sqrt(eq_t_ev(ip,0)*sml_ev2j/sml_prot_mass/2.0)

        ! e14.5 type : 
        write(io_FP4Dout,'(t2,a,t21,e14.5)') 'rho = ', eq_psin(ip)
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'inv_dpsidr = ', inv_dpsidr(ip)
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'normalized jparB = ', FSA_jparB(ip) / eq_den(ip,0) / sml_ev2j / vth_tmp / B0
        ! es14.5 type
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'sum Z_s Gamma_s = ', sum_pflux(ip)*inv_dpsidr(ip)/eq_den(ip,0)/vth_tmp
        write(io_FP4Dout,'(t3,a,t11,a,t25,a,t39,a)') &
            'Z', 'pflux', 'eflux', 'mflux'

        do is=is_s, is_e
          write (io_FP4Dout,'(f7.3,3(es14.5))') &
                sml_charge_num(is), &
                FSA_pflux(ip,is)*inv_dpsidr(ip) &
                  /eq_den(ip,0)/vth_tmp, &                          ! n_norm * v_norm
                FSA_eflux(ip,is)*inv_dpsidr(ip) &
                  /eq_den(ip,0)/vth_tmp/eq_t_ev(ip,0)/sml_ev2j, &   ! n_norm * v_norm * T_norm
                FSA_mflux(ip,is)*inv_dpsidr(ip) &
                  /eq_den(ip,0)/eq_t_ev(ip,0)/sml_ev2j/a0           ! n_norm * T_norm * a_norm
        enddo

      enddo

      write(io_FP4Dout,*) '------------------'
      write(io_FP4Dout,*)

    close(io_FP4Dout)



    open(unit=io_FP4Dout, file=runfile_FP4Dout, status='old', position='append')     

      sum_pflux=0.0
      do ip=ip_s,ip_e
        do is=is_s, is_e
          sum_pflux(ip) = sum_pflux(ip) + sml_charge_num(is) * FSA_pflux_col(ip,is)
        enddo
      enddo

      write(io_FP4Dout,*) "-------------------------------"
      write(io_FP4Dout,*) "----   COLLSIONAL FLUX ??  ----"
      write(io_FP4Dout,*) "-------------------------------"
      write(io_FP4Dout,*) 'To compare with NEO FLUX'
      write(io_FP4Dout,*) '  from psi to r'
      write(io_FP4Dout,*) '  normalizatino : refer to out.fp4d.neo_input'
      write(io_FP4Dout,*) "-------------------------------"

      inv_dpsidr = 0D0
      do ip=ip_s,ip_e
        do ith=ith_s,ith_e
          inv_dpsidr(ip) = inv_dpsidr(ip) + 1.0 / eq_Bp(ith,ip)/geo_R(ip,ith)
        enddo
        inv_dpsidr(ip) = inv_dpsidr(ip) / nth
      enddo

      do ip=ip_s,ip_e

        vth_tmp = sqrt(eq_t_ev(ip,0)*sml_ev2j/sml_prot_mass/2.0)

        ! e14.5 type : 
        write(io_FP4Dout,'(t2,a,t21,e14.5)') 'rho = ', eq_psin(ip)
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'inv_dpsidr = ', inv_dpsidr(ip)
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'normalized jparB = ', FSA_jparB(ip) / eq_den(ip,0) / sml_ev2j / vth_tmp / B0
        ! es14.5 type
        write(io_FP4Dout,'(t2,a,t21,es14.5)') 'sum Z_s Gamma_s = ', sum_pflux(ip)*inv_dpsidr(ip)/eq_den(ip,0)/vth_tmp
        write(io_FP4Dout,'(t3,a,t11,a,t25,a,t39,a)') &
            'Z', 'pflux', 'eflux', 'mflux'

        do is=is_s, is_e
          write (io_FP4Dout,'(f7.3,3(es14.5))') &
                sml_charge_num(is), &
                FSA_pflux_col(ip,is)*inv_dpsidr(ip) &
                  /eq_den(ip,0)/vth_tmp, &                              ! n_norm * v_norm
                FSA_eflux_col(ip,is)*inv_dpsidr(ip) &
                  /eq_den(ip,0)/vth_tmp/eq_t_ev(ip,0)/sml_ev2j, &       ! n_norm * v_norm * T_norm
                FSA_mflux_col(ip,is)*inv_dpsidr(ip) &
                  /eq_den(ip,0)/eq_t_ev(ip,0)/sml_ev2j/a0               ! n_norm * T_norm * a_norm
        enddo

      enddo

      write(io_FP4Dout,*) '------------------'
      write(io_FP4Dout,*)

    close(io_FP4Dout)

  end subroutine Post_record

  subroutine Post_finalize()
    implicit none




    deallocate(pflux)
    deallocate(eflux)
    deallocate(mflux)
    deallocate(pflux_col)
    deallocate(eflux_col)
    deallocate(mflux_col)

    deallocate(FSA_pflux)
    deallocate(FSA_eflux)
    deallocate(FSA_mflux)

    deallocate(FSA_pflux_col)
    deallocate(FSA_eflux_col)
    deallocate(FSA_mflux_col)

    deallocate(post_dens)
    deallocate(post_NUpar)
    deallocate(post_NUparB)

    deallocate(FSA_post_dens)
    deallocate(FSA_post_NUpar)
    deallocate(FSA_post_NUparB)

    deallocate(jpar)
    deallocate(jparB)
    deallocate(FSA_jpar)
    deallocate(FSA_jparB)


  endsubroutine Post_finalize

end module FP4D_Post