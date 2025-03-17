module FP4D_equilibrium

  use FP4D_globals
  use FP4D_math
  use eqmodule
  contains

  subroutine eq_init()

    use mpi
    implicit none
    ! integer, intent(in) :: isp, nsp, sml_eqb_opt, npsi, nth

    integer :: ip, ith, is, js, profile_npsi

    ! Allocate memory
    allocate(eq_psin   (ip_s:ip_e),            &
            ! eq_omega(ith_s:ith_e,ip_s:ip_e,isp:nsp), &
             eq_omega  (ith_s:ith_e, ip_s:ip_e, tar_s:tar_e), & ! target
             eq_theta  (ith_s:ith_e),                       &
             eq_invJ   (ith_s:ith_e, ip_s:ip_e),    eq_I     (ith_s:ith_e, ip_s:ip_e),      &
             eq_J      (ith_s:ith_e, ip_s:ip_e),    eq_B     (ith_s:ith_e, ip_s:ip_e),      &
             eq_Bt     (ith_s:ith_e, ip_s:ip_e),    eq_Bp    (ith_s:ith_e, ip_s:ip_e),      &
             eq_dBdth  (ith_s:ith_e, ip_s:ip_e),    eq_dBdpsi(ith_s:ith_e, ip_s:ip_e),      &
             eq_phi1   (ith_s:ith_e, ip_s:ip_e),    eq_dphidth   (ith_s:ith_e, ip_s:ip_e),  &
             eq_phi1_k1(ith_s:ith_e, ip_s:ip_e),    eq_dphidth_k1(ith_s:ith_e, ip_s:ip_e),  &
             eq_phi1_k2(ith_s:ith_e, ip_s:ip_e),    eq_dphidth_k2(ith_s:ith_e, ip_s:ip_e),  &
             eq_phi1_k3(ith_s:ith_e, ip_s:ip_e),    eq_dphidth_k3(ith_s:ith_e, ip_s:ip_e),  &
             eq_phi1_tmp(ith_s:ith_e, ip_s:ip_e), &
            !  eq_t_ev   (ip_s:ip_e,  isp:nsp),   eq_den    (ip_s:ip_e,  isp:nsp),    &
            !  eq_dlnT   (ip_s:ip_e,isp:nsp),     eq_dlnn   (ip_s:ip_e,  isp:nsp),    &
             eq_t_ev   (ip_s:ip_e,  tar_s:tar_e),     eq_den    (ip_s:ip_e,  tar_s:tar_e),    &
             eq_dlnT   (ip_s:ip_e,  tar_s:tar_e),     eq_dlnn   (ip_s:ip_e,  tar_s:tar_e),    &

             eq_lambda (ip_s:ip_e,  tar_s:tar_e, tar_s:tar_e), &
             eq_tau    (ip_s:ip_e,  tar_s:tar_e, tar_s:tar_e), &
             eq_nu     (ip_s:ip_e,  tar_s:tar_e, tar_s:tar_e), &

             geo_R    (ip_s:ip_e, ith_s:ith_e),     geo_Z    (ip_s:ip_e, ith_s:ith_e),      & !! caution :: (ip,ith) due to eq_util
             geo_dRdth(ip_s:ip_e, ith_s:ith_e), &
             geo_FSA_R2 (ip_s:ip_e), &

             geo_psi  (ip_s:ip_e, ith_s:ith_e),         &
             eq_intJ   (ip_s:ip_e), eq_dpsidrho(ip_s:ip_e),&
             eq_Volume (ith_s:ith_e, ip_s:ip_e))

    eq_psin  = 0D0; eq_theta = 0D0
    eq_Bt   = 0D0; eq_Bp    = 0D0; eq_B     = 0D0
    eq_dBdth= 0D0; eq_dBdpsi= 0D0; eq_J     = 0D0
    eq_invJ = 0D0; eq_omega = 0D0; eq_intJ  = 0D0; eq_Volume=0D0; eq_dpsidrho=0D0
    eq_I     = 0D0;
    eq_phi1 = 0D0; eq_dphidth =0D0;
    eq_phi1_k1 = 0D0; eq_dphidth_k1 =0D0;
    eq_phi1_k2 = 0D0; eq_dphidth_k2 =0D0;
    eq_phi1_k3 = 0D0; eq_dphidth_k3 =0D0;

    eq_lambda = 0D0; eq_tau = 0D0; eq_nu = 0D0; 


    geo_FSA_R2 = 0D0;

    !-----------------------------------------------------------------------
    ! Junhyuk Song added the 'eqmodule' from eqdsk
    ! 
    call set_eq
    !
    ! write(*,*) 'TEST eqmodule'
    ! write(*,*) 'shape of SHAPE(R_arr)', SHAPE(R_arr)
    ! write(*,*) 'R_arr', R_arr
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    ! define using input
    eq_dtheta = sml_2pi/real(nth,8)

    if (ip_e .eq. ip_s) then
      eq_drho = xx1-xx0             ! for Npsi=Nth=1 case...  NB case...
      eq_psin(ip_s) = xx0
    else
      eq_drho = (xx1-xx0)/(1D0*npsi-1D0)
      do ip = ip_s,ip_e
        eq_psin(ip) = xx0 + eq_drho * (ip*1D0-1D0)
      end do
    endif

    write(*,*) 'eq_psin', eq_psin
    write(*,*) 'deltapsi', deltapsi

    eq_dpsidrho = 2.0 * eq_psin * deltapsi   !  J = 1/(B cdot grad th)

    write(*,*) 'eq_dpsidrho', eq_dpsidrho

    do ith = ith_s,ith_e
      eq_theta(ith)=real(ith,8)*eq_dtheta
    end do
    !-----------------------------------------------------------------------



    !-----------------------------------------------------------------------
    ! extract the data from 'eqmodule' 
    !
    select case(sml_eqb_opt)
      case (0) ! using EqUtil
      case (1) ! using text-file

        do ip=ip_s,ip_e
          do ith=ith_s,ith_e
            eq_Bt    (ith,ip) = B_arr(ip,ith)
            eq_Bp    (ith,ip) = Bp_arr(ip,ith)
            eq_I     (ith,ip) = I_arr(ip,ith)
            eq_J     (ith,ip) = abs(J_arr(ip,ith))
            eq_dBdpsi(ith,ip) = dBdp_arr(ip,ith)
            eq_dBdth (ith,ip) = dBdt_arr(ip,ith)

            eq_B     (ith,ip) = sqrt(eq_Bt(ith,ip)**2.0+eq_Bp(ith,ip)**2.0)
            eq_invJ  (ith,ip) = 1D0/eq_J(ith,ip)
          enddo
        enddo


        !-----------------------------------------------------------------------
        ! write(*,*) 'B_arr', shape(B_arr)
        ! write(*,*) B_arr

        ! write(*,*) 'Bp_arr', shape(Bp_arr)
        ! write(*,*) Bp_arr

        ! write(*,*) 'I_arr', shape(I_arr)
        ! write(*,*) I_arr

        ! write(*,*) 'J_arr', shape(J_arr)
        ! write(*,*) J_arr

        ! write(*,*) 'dBdp_arr', shape(dBdp_arr)
        ! write(*,*) dBdp_arr


        ! write(*,*) 'dBdt_arr', shape(dBdt_arr)
        ! write(*,*) dBdt_arr
        !-----------------------------------------------------------------------


        !-----------------------------------------------------------------------
        ! write(*,*) 'eq_Bt', shape(eq_Bt)
        ! write(*,*) eq_Bt

        ! write(*,*) 'eq_Bp', shape(eq_Bp)
        ! write(*,*) eq_Bp

        ! write(*,*) 'eq_I', shape(eq_I)
        ! write(*,*) eq_I

        ! write(*,*) 'eq_J', shape(eq_J)
        ! write(*,*) eq_J

        ! write(*,*) 'eq_dBdpsi', shape(eq_dBdpsi)
        ! write(*,*) eq_dBdpsi

        ! write(*,*) 'eq_dBdth', shape(eq_dBdth)
        ! write(*,*) eq_dBdth
        !-----------------------------------------------------------------------


        do ip=ip_s,ip_e
          do ith=ith_s,ith_e
            eq_intJ(ip) = eq_intJ(ip) + eq_J(ith,ip) * eq_dtheta
          enddo
        enddo


        do ip=ip_s,ip_e
          do ith=ith_s,ith_e
            eq_Volume(ith,ip) = eq_dpsidrho(ip) * eq_J(ith,ip) & ! for rho not psi
                               * eq_dtheta * eq_drho * sml_2pi   ! include toroidal direction
          enddo
        enddo

        write(*,*) 'eq_drho', eq_drho

        ! write(*,*) 'eq_Volume', eq_Volume


        ! youknow. tmp code. for using TORIC Jacobian
        ! if (opt_QLD.eq.2) then
        !   open(1100,file='arr_J.txt',action='read')
        !   do ip=ip_s,ip_e
        !     read(1100,*) eq_J(:,ip)
        !   end do
        !   close(1100)
        !   write(*,*) 'eq_J', eq_J
        ! endif
        
        ! Correction of Numerical Difference for dBdth
        eq_dBdth(0,:)=0D0
        eq_dBdth(nth/2,:)=0D0
        ! do is = isp,nsp
        do is = tar_s,tar_e ! tar

          do ip=ip_s,ip_e
            eq_omega(:,ip,is)=sml_charge(is)*eq_B(:,ip)/sml_mass(is)
          end do
        end do

        !!! REPLACE 'arr_RZ.txt' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  geo_R [ip_s:ip_e  , ith_s:ith_e]
        !  R_arr [ip_s:ip_e+1, ith_s:ith_e]

        geo_R  (ip_s:ip_e,:) = R_arr(ip_s:ip_e,:)
        geo_Z  (ip_s:ip_e,:) = Z_arr(ip_s:ip_e,:)
        geo_psi  (ip_s:ip_e,:) = psi_arr(ip_s:ip_e,:)

        ! write(*,*) 'R_arr', shape(R_arr)
        ! write(*,*) R_arr

        ! write(*,*) 'geo_R', shape(geo_R)
        ! write(*,*) geo_R


        ! write(*,*) 'Z_arr', shape(Z_arr)
        ! write(*,*) Z_arr

        ! write(*,*) 'geo_Z', shape(geo_Z)
        ! write(*,*) geo_Z


        ! write(*,*) 'psi_arr', shape(psi_arr)
        ! write(*,*) psi_arr

        ! write(*,*) 'geo_psi', shape(geo_psi)
        ! write(*,*) geo_psi
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) 'filename_profile :', filename_profile

        geo_dRdth = 0D0
        do ip=ip_s,ip_e
          call cal_dfdth6(geo_R(ip,:), geo_dRdth(ip,:))
          do ith=ith_s,ith_e
            geo_FSA_R2(ip) = geo_FSA_R2(ip)     &
                           + geo_R(ip,ith)**2.0 &
                           * eq_J(ith,ip) * eq_dtheta &
                           / eq_intJ(ip)
          enddo
        enddo

        if (opt_global.eq.1) then
        !Read "profile.dat"
          open(io_profile,file=filename_profile,action='read')
          read(io_profile,*) profile_npsi
          write(*,*) "npsi of profile data :", profile_npsi

          ! allocate(eq_profile(1:profile_npsi,5),prf_t_ev(isp:nsp,1:profile_npsi),&
          !   prf_den(isp:nsp,1:profile_npsi),prf_dlnn_dr(isp:nsp,1:profile_npsi),&
          !   prf_dlnT_dr(isp:nsp,1:profile_npsi))
          ! youknow : maxwellian species..
          allocate(eq_profile(1:profile_npsi,5),prf_t_ev(tar_s:tar_e,1:profile_npsi),&
            prf_den(tar_s:tar_e,1:profile_npsi),prf_dlnn_dr(tar_s:tar_e,1:profile_npsi),&
            prf_dlnT_dr(tar_s:tar_e,1:profile_npsi))
          eq_profile=0D0

          do ip=1,profile_npsi
            read(io_profile,*) eq_profile (ip,:)
          enddo

          close(io_profile)
          prf_den(0,:)=eq_profile(:,2)*1D19 
          prf_den(1,:)=eq_profile(:,2)*1D19
          prf_t_ev(0,:)=eq_profile(:,3)*1D3
          prf_t_ev(1,:)=eq_profile(:,4)*1D3

          if (is_e.gt.1) then
            ! prf_den(2,:)=eq_profile(:,2)*1D19
            ! prf_den(2,:)=eq_profile(:,2)*1D19/5D4 ! tmp for NB (max. e, max. D, fast D) 
            
            ! prf_den(2,:)=eq_profile(:,2)*1D19/5D4/2D0 ! tmp for NB (max. e, max. D, fast D) 
            prf_den(2,:)=eq_profile(:,2)*1D19*NB_fast_den_factor ! NB_fast_den_factor :: (default) 1D-5
            
            ! write(*,*) 'flag1', prf_den (2,:)
            ! write(*,*) 'flag1', eq_profile (:,2)

            ! prf_den(2,:)=eq_profile(:,2)*1D19/5D4/2D0/1D0 ! tmp for NB (max. e, max. D, fast D) 

            ! prf_t_ev(2,:)=eq_profile(:,3)*1D3
            prf_t_ev(2,:)=eq_profile(:,4)*1D3

          endif

          do is=is_s,is_e
            do ip=1,profile_npsi
              if (ip.eq.1) then
                prf_dlnT_dr(is,ip)=(prf_t_ev(is,ip+1)-prf_t_ev(is,ip)) &
                                    /(eq_profile(ip+1,1)-eq_profile(ip,1))/prf_t_ev(is,ip)
                prf_dlnn_dr(is,ip)=(prf_den(is,ip+1)-prf_den(is,ip)) &
                                    /(eq_profile(ip+1,1)-eq_profile(ip,1))/prf_den(is,ip)
              elseif (ip.eq.profile_npsi) then
                prf_dlnT_dr(is,ip)=(prf_t_ev(is,ip)-prf_t_ev(is,ip-1)) &
                                    /(eq_profile(ip,1)-eq_profile(ip-1,1))/prf_t_ev(is,ip)
                prf_dlnn_dr(is,ip)=(prf_den(is,ip)-prf_den(is,ip-1)) &
                                    /(eq_profile(ip,1)-eq_profile(ip-1,1))/prf_den(is,ip) 
              else
                prf_dlnT_dr(is,ip)=(prf_t_ev(is,ip+1)-prf_t_ev(is,ip-1)) &
                                    /(eq_profile(ip+1,1)-eq_profile(ip-1,1))/prf_t_ev(is,ip)
                prf_dlnn_dr(is,ip)=(prf_den(is,ip+1)-prf_den(is,ip-1)) &
                                    /(eq_profile(ip+1,1)-eq_profile(ip-1,1))/prf_den(is,ip)
              endif
            enddo
          enddo

          ! Spline
          write(*,*) 'xx0', xx0
          write(*,*) 'xx1', xx1

          ! x = eq_profile(:,1)
          ! y = prf_den(:,is), ...
          ! n = profile_npsi
          ! xx = eq_psin(:)
          ! yy = ...
          ! nn = nsp
          do is = tar_s, tar_e ! considering maxweliian sp.
            ! profile.data :: rho definition is rho =  (Psi-Psi_0)/(Psi_b-Psi_0) 
            if (opt_profile_rho .eq. 1) then
              ! default
              ! if FP4D use 'rho = ( (Psi-Psi_0)/(Psi_b-Psi_0) )' then use this option.
              call spline_psi (eq_profile(:,1), prf_den    (is,:), profile_npsi, eq_psin, eq_den  (:,is), npsi)
              call spline_psi (eq_profile(:,1), prf_dlnn_dr(is,:), profile_npsi, eq_psin, eq_dlnn (:,is), npsi)
              call spline_psi (eq_profile(:,1), prf_t_ev   (is,:), profile_npsi, eq_psin, eq_t_ev (:,is), npsi)
              call spline_psi (eq_profile(:,1), prf_dlnT_dr(is,:), profile_npsi, eq_psin, eq_dlnT (:,is), npsi)
            elseif (opt_profile_rho .eq. 2) then
              ! rho =  (Psi-Psi_0)/(Psi_b-Psi_0)               
              ! if FP4D use 'rho = sqrt( (Psi-Psi_0)/(Psi_b-Psi_0) )' then use this option.
              call spline_psi (sqrt(eq_profile(:,1)), prf_den    (is,:), profile_npsi, eq_psin, eq_den  (:,is), npsi)
              call spline_psi (sqrt(eq_profile(:,1)), prf_dlnn_dr(is,:), profile_npsi, eq_psin, eq_dlnn (:,is), npsi)
              call spline_psi (sqrt(eq_profile(:,1)), prf_t_ev   (is,:), profile_npsi, eq_psin, eq_t_ev (:,is), npsi)
              call spline_psi (sqrt(eq_profile(:,1)), prf_dlnT_dr(is,:), profile_npsi, eq_psin, eq_dlnT (:,is), npsi)
            else
              write(*,*) 'error :: opt_profile_rho'
            endif
          end do
          ! write(*,*)  'eq_profile', eq_profile
          ! write(*,*)  'profile_npsi', profile_npsi
          write(*,*)  'eq_psin', eq_psin
          do is = tar_s, tar_e
            write(*,*)  'isp',is,'eq_den',eq_den (:,is)  
          enddo     

          ! if (npsi.ge.1) then
          !   ! do ip=ip_s,ip_e
          !     ! do is=isp,nsp
          !     do is=0,nsp

          !       eq_den (ip,is)=1D0*prf_den    (is,ip*(profile_npsi/(npsi+1)))
          !       write(*,*) 'flag2', 'ip', 'is', eq_den (ip,is)
          !       eq_dlnn(ip,is)=1D0*prf_dlnn_dr(is,ip*(profile_npsi/(npsi+1)))
          !       eq_t_ev(ip,is)=1D0*prf_t_ev   (is,ip*(profile_npsi/(npsi+1)))
          !       eq_dlnT(ip,is)=1D0*prf_dlnT_dr(is,ip*(profile_npsi/(npsi+1)))
          !     enddo
          !     eq_psin(ip)=1D0*eq_profile(ip*(profile_npsi/(npsi+1)),1)
          !   enddo
          ! else
          !   write(*,*) "Error:! sml_npsi must be greater than 1!"
          ! endif
        else if (opt_global.eq.0) then
          ! do is=isp,nsp
          do is = tar_s, tar_e
            eq_den (:,is)=1D0*sml_den (is)
            eq_t_ev(:,is)=1D0*sml_t_ev(is)
            eq_dlnn(:,is)=1D0*dlnn_dr (is)
            eq_dlnT(:,is)=1D0*dlnT_dr (is) 
          enddo
        endif
      case default
    end select     ! select case(sml_eqb_opt)

    do is = tar_s, tar_e  ! colliding
    do js = tar_s, tar_e  ! target
    do ip = ip_s, ip_e
        ! eq_lambda, eq_tau, eq_nu :: (ip, target, colliding)
          ! a : colliding --> is
          ! b : target    --> js
        call cal_tau(sml_charge_num(is), sml_charge_num(js), &
                      sml_mass(is), sml_mass(js),                     &
                      eq_den(ip,is), eq_den(ip,js),           &
                      eq_t_ev(ip,is), eq_t_ev(ip,js),         &
                      eq_lambda(ip,js,is),                   &
                      eq_tau(ip,js,is),                      &
                      eq_nu(ip,js,is))
    enddo
    enddo
    enddo

    !-----------------------------------------------------------------------
    ! record the profile to out.fp4d.run & out.fp4d.neo_input
    if (sml_mype==0) then
      ! out.fp4d.run
      call eq_record()

      ! out.fp4d.neo_input
      call eq_record_NEO()
    endif
    !-----------------------------------------------------------------------


    if (special_tau_W.eq.1) then
      write(*,*) ' '
      write(*,*) 'Caution :: special_tau_W = 1, tau definition is changed. Be careful when you post-process'
      write(*,*) ' '
      write(*,*) 'before tau', eq_tau(1,1,2) 
      eq_tau(1,1,2) = eq_tau(1,1,1) ! ip=1
      write(*,*) 'after  tau', eq_tau(1,1,2) 
    endif

    if (special_tau_He3.eq.1) then
      write(*,*) ' '
      write(*,*) 'Caution :: special_tau_He3 = 1, tau definition is changed. Be careful when you post-process'
      write(*,*) ' '
      write(*,*) 'before tau', eq_tau(1,1,3) 
      eq_tau(1,1,3) = eq_tau(1,1,1) ! ip=1
      write(*,*) 'after  tau', eq_tau(1,1,3) 
    endif


    call deallocate_eqmodule

    return

  end subroutine eq_init

  subroutine eq_record()
    implicit none
    integer :: ip, ith, is, js
    character(len=7) :: kind_str
        open(unit=io_FP4Dout, file=runfile_FP4Dout, status='old', position='append')     
          write(io_FP4Dout,*) ''      
          write(io_FP4Dout,*) 'rho : sqrt(normalized pololidal flux)'
          write(io_FP4Dout,*) ''      

          do ip = ip_s, ip_e
              write(io_FP4Dout, '(A, I2)') 'ip : ', ip
              write(io_FP4Dout, '(t2,a,t21,e14.5)') 'rho = ', eq_psin(ip)
              write(io_FP4Dout, '(A)') ""

              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "GEOMETRY PARAMETERS -> Future work"
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "rho, geometry option(efit,circular),"
              write(io_FP4Dout, '(A)') "R0,q0,roation,"
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "grid resolution..."
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "RF option... RF power..."
              write(io_FP4Dout, '(A)') "-------------------------------"

              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "EQUILIBRIUM PARAMETERS"
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "1/Ln, 1/Lt for r"

              write(io_FP4Dout,'(a)') &
                  'idx    kind     Z       n[m^-3]      T[eV]        m/m_H         1/Ln         1/Lt'
              do is = tar_s, tar_e
                if (is < is_s) then
                    kind_str = "Maxwell"
                else
                    kind_str = "General"
                end if

              write(io_FP4Dout,'(t2,i2,2x,a7,2x,f7.3,2x,6(1pe11.4,2x))') &
                      is,kind_str,sml_charge_num(is),eq_den(ip,is),eq_t_ev(ip,is),sml_mass(is)/sml_prot_mass,eq_dlnn(ip,is),eq_dlnT(ip,is)
              enddo



              ! eq_lambda
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "ln_Lambda_ab"
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A, *(I4,1X))') "       ", (js, js=tar_s,tar_e)

              do is = tar_s, tar_e
                  write(io_FP4Dout, '(I2, A, *(F10.4))') is, " |", (eq_lambda(ip, js, is), js=tar_s,tar_e)
              end do

              ! eq_tau

              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "tau_ab [same with Helander]"
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A, *(I4,1X))') "       ", (js, js=tar_s,tar_e)

              do is = tar_s, tar_e
                  write(io_FP4Dout, '(I2, A, *(ES15.6))') is, " |", (eq_tau(ip, js, is), js=tar_s,tar_e)
              end do

              ! eq_nu
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "nu_ab [same with Helander]"
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A, *(I4,1X))') "       ", (js, js=tar_s,tar_e)

              do is = tar_s, tar_e
                  write(io_FP4Dout, '(I2, A, *(ES15.6))') is, " |", (eq_nu(ip, js, is), js=tar_s,tar_e)
              end do
              write(io_FP4Dout, '(A)') "-------------------------------"
              write(io_FP4Dout, '(A)') "*****equilibrium end*******"
          end do

        close(io_FP4Dout)


  end subroutine eq_record

  subroutine eq_record_NEO()
    implicit none
    integer :: ip, ith, is, js
    real(kind=8) :: vth_tmp, rho_star_tmp

        open(unit=io_NEO, file=runfile_NEO_input, status='replace')
        close(io_NEO)
        
        open(unit=io_NEO, file=runfile_NEO_input, status='old', position='append')     
          write(io_NEO,*) ''      
          write(io_NEO,*) 'rho : sqrt(normalized pololidal flux)'
          write(io_NEO,*) ''      

          do ip = ip_s, ip_e
              write(io_NEO, '(A, I2)') 'ip : ', ip
              write(io_NEO, '(t2,a,t21,e14.5)') 'rho = ', eq_psin(ip)
              write(io_NEO, '(A)') ""

              write(io_NEO, '(A)') "-------------------------------"
              write(io_NEO, '(A)') "GEOMETRY PARAMETERS -> Future work"
              write(io_NEO, '(A)') "-------------------------------"
              write(io_NEO, '(A)') "rho, geometry option(efit,circular),"
              write(io_NEO, '(A)') "R0,q0,roation,"
              write(io_NEO, '(A)') "-------------------------------"
              write(io_NEO, '(A)') "grid resolution..."
              write(io_NEO, '(A)') "-------------------------------"

              write(io_NEO, '(A)') "-------------------------------"
              write(io_NEO, '(A)') "EQUILIBRIUM PARAMETERS"
              write(io_NEO, '(A)') "-------------------------------"
              write(io_NEO, '(A)') "NEO normalization :"
              write(io_NEO, '(A)') "  n_norm = n(0), T_norm = T(0), m_norm = m_Deuterium"
              write(io_NEO, '(A)') 'vth_norm = sqrt(T_norm/m_norm)'

              write(io_NEO, '(A12, ES15.6)') 'n_norm:', eq_den(ip,0)
              write(io_NEO, '(A12, ES15.6)') 'T_norm:', eq_t_ev(ip,0)
              write(io_NEO, '(A12, ES15.6)') 'm_norm:', sml_prot_mass*2.0
              write(io_NEO, '(A12, ES15.6)') 'B_norm:', B0 ! eq_B(nth/4,ip) ::need to modified...
              write(io_NEO, '(A12, ES15.6)') 'a_norm:', a0

              ! NEO_vth :: sqrt(T(0)/m_D)
              vth_tmp = sqrt(eq_t_ev(ip,0)*sml_ev2j/sml_prot_mass/2.0)
              rho_star_tmp = (sml_prot_mass*2.0)*vth_tmp/sml_ev2j/B0 ! eq_B(nth/4,ip) rather tahn B0 ...
              rho_star_tmp = rho_star_tmp/a0

              write(io_NEO, '(A, ES15.6)') 'vth_norm:', vth_tmp
              write(io_NEO, '(A)') ""
              write(io_NEO, '(A)') ""
              write(io_NEO, '(A)') "NU_1 define well if Lambda_ei=ee & First species must be electron"
              write(io_NEO, '(A)') "omega_rot define in f_module... so, need to modify this ..."
              write(io_NEO, '(A12, ES15.6)') 'NU_1:',      eq_nu(ip, 0, 0) /vth_tmp*a0
              write(io_NEO, '(A12, ES15.6)') 'RHO_STAR:',  rho_star_tmp
              write(io_NEO, '(A12, ES15.6)') 'OMEGA_ROT:', omega_rot_in/(vth_tmp/a0)
              write(io_NEO, '(A)') ""
              write(io_NEO, '(A)') ""

              write(io_NEO,'(a)') &
                  'idx    Z       n/n_norm     T/T_norm     m/m_D        -a/Ln        -a/Lt        '
              do is = tar_s, tar_e


                write(io_NEO,'(t2,i2,2x,f7.3,2x,6(1pe11.4,2x))') &
                        is,        &
                        sml_charge_num(is), &
                        eq_den(ip,is)     /eq_den(ip,0),      &
                        eq_t_ev(ip,is)    /eq_t_ev(ip,0),     &
                        sml_mass(is)      /sml_prot_mass/2.0, &
                        -eq_dlnn(ip,is)   *a0,     &
                        -eq_dlnT(ip,is)   *a0
              enddo

              write(io_NEO, '(A)') "-------------------------------"
          end do

        close(io_NEO)


  end subroutine eq_record_NEO



  subroutine cal_tau(Za,Zb,ma,mb,na,nb,Ta,Tb, &
                     ln_lambda,tau,nu)
    real(kind=8), intent(in)  :: Za, Zb, ma, mb, na, nb, Ta, Tb
    real(kind=8)              :: Zi, Ze, mi, me, ni, ne, Ti, Te
    real(kind=8)              :: eps0, eV

    real(kind=8), intent(out) :: ln_lambda, tau, nu
    ! a :: colliding
    ! b :: target

    eps0 = 8.8542D-12
    eV = 1.6022D-19
    
    ! NRL Plasma Formulary 2011 version, p.34
    ! It has only difference for Te<Ti*me/mi with NRL Plasma Formulary 2018 version, p.34
    ! me/mi ~ 1/3600 for m_D ==> we can ignore the difference

    if ( (Za.lt.0) .or. (Zb.lt.0) ) then
        if ( (Za.lt.0) .and. (Zb.lt.0) ) then
            ! thermal ee collision
            ln_lambda = 23.5 - log(sqrt(na*1D-6)*(Ta**(-1.25))) &
                      - sqrt(abs(1D-5+((log(Ta)-2.0)**2.0)/16.0))
        else
            ! ie | ei
            if (Za.lt.0) then
                ! a: electron | b: ion
                Ze=Za;me=ma;ne=na;Te=Ta
                Zi=Zb;mi=mb;ni=nb;Ti=Tb
            else
                ! a: ion | b: electron
                Zi=Za;mi=ma;ni=na;Ti=Ta                               
                Ze=Zb;me=mb;ne=nb;Te=Tb
            endif

            if (Ti * me/mi .lt. Te) then
                if(Te .lt. (10.0*(Zi**2.0))) then
                    ln_lambda = 23.0 - log(Zi*sqrt(ne*1D-6/(Te**3.0)))
                else
                    ln_lambda = 24.0 - log(sqrt(ne*1D-6)/Te)
                endif
            else
                ln_lambda = 30.0 - log((Zi**2.0)*sqrt(ni*1D-6/(Ti**3.0))/mi)
            endif
        endif

    else
        ! ii collision
        ln_lambda = 23.0 - log( Za * Zb * (ma + mb)  &
                    / (ma*Tb + mb*Ta)                &
                    * sqrt( na*1D-6*(Za**2.0)/Ta     &
                    + nb*1D-6*(Zb**2.0)/Tb) )
    endif
        
    ! tau_aa_Helander
    tau = 12.0*sqrt(sml_pi)**3.0 / sqrt(2.0)        &
            * sqrt(ma) * sqrt(Ta*eV)**3.0          &
            * eps0**2.0                             &
            / nb / Za**2.0 / Zb**2.0 / eV**4.0   &
            / ln_lambda
    ! nu_aa_Helander
    nu  = 3.0*sqrt(sml_pi)/4.0                          &
          / tau

    ! write(*,*) 'inner info.1', Za, Zb, ma, mb, na, nb, Ta, Tb
    ! write(*,*) 'inner info.2', Zi, Ze, mi, me, ni, ne, Ti, Te
    ! write(*,*) 'inner tau', ln_lambda
    ! write(*,*) 'inner tau', tau
    ! write(*,*) 'inner nu' , nu


  end subroutine cal_tau

  !> Deallocate f0 memory
  subroutine eq_finalize
    implicit none

    deallocate(eq_psin,eq_theta,eq_omega,eq_invJ,eq_I,eq_J,&
               eq_B,eq_Bt,eq_Bp, eq_t_ev,eq_den,eq_dlnn,eq_dlnT,eq_dBdth,&
               eq_dBdpsi,&
               geo_R,geo_Z,geo_psi,&
               geo_FSA_R2,&
               eq_lambda,eq_tau,eq_nu,&
               eq_intJ,eq_phi1,eq_dphidth,&
               eq_phi1_k1,eq_phi1_k2,eq_phi1_k3,&
               eq_dphidth_k1,eq_dphidth_k2,eq_dphidth_k3,eq_dpsidrho,eq_Volume)
    if (opt_global.eq.1) then
      deallocate(eq_profile,prf_den,prf_t_ev,prf_dlnT_dr,prf_dlnn_dr)
    endif

    return

  end subroutine eq_finalize
  
  subroutine spline_psi (x, y, n, xx, yy, nn)
    integer, intent(in) :: n, nn
    real (kind=8), intent(in)   :: x(n),   y(n)
    real (kind=8), intent(in)   :: xx(nn)
    real (kind=8), intent(out)  :: yy(nn)

    real (kind=8) :: b(n), c(n), d(n)
    integer :: i

    call spline (x, y, b, c, d, n)
    
    do i = 1, nn
        call ispline(xx(i), x, y, b, c, d, n, yy(i))
        ! write(*,*) 'i',i,'xx(i)',xx(i),'yy(i)',yy(i)
    end do

  endsubroutine spline_psi

    subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
      implicit none
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer i, j, gap
      double precision h

      gap = n-1
!     check input
      if ( n < 2 ) return
      if ( n < 3 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1)) ! linear interpolation
         c(1) = 0.
         d(1) = 0.
         b(2) = b(1)
         c(2) = 0.
         d(2) = 0.
         return
      end if
!
! step 1: preparation
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, gap
         d(i) = x(i+1) - x(i)
         b(i) = 2.0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      end do
!
! step 2: end conditions 
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if(n /= 3) then
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if
!
! step 3: forward elimination 
!
      do i = 2, n
         h = d(i-1)/b(i-1)
         b(i) = b(i) - h*d(i-1)
         c(i) = c(i) - h*c(i-1)
      end do
!
! step 4: back substitution
!
      c(n) = c(n)/b(n)
      do j = 1, gap
         i = n-j
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
!
! step 5: compute spline coefficients
!
      b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
      do i = 1, gap
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      end subroutine spline

      subroutine ispline(u, x, y, b, c, d, n, uy)
!======================================================================
! subroutine ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! uy = interpolated value at point u
!=======================================================================
      implicit none
      integer n
      double precision  u, uy, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      double precision dx

      save

! if u is ouside the x() interval take a boundary value (left or right)
      if(u <= x(1)) then
         dx = u - x(1)
         uy = y(1) + dx*(b(1) + dx*(c(1) + dx*d(1)))
!         uy = y(1)
         return
      end if
      if(u >= x(n)) then
         dx = u - x(n)
         uy = y(n) + dx*(b(n) + dx*(c(n) + dx*d(n)))
!         uy = y(n)
         return
      end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
      i = 1
      j = n+1
      do while (j > i+1)
         k = (i+j)/2
         if(u < x(k)) then
            j=k
         else
            i=k
         end if
      end do
!*
!  evaluate spline interpolation
!*
      dx = u - x(i)
      uy = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      end subroutine ispline



end module FP4D_equilibrium
