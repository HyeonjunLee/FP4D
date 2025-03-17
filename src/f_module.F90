module f_module

  use FP4D_globals
  use FP4D_math

  contains

  !> Initialize global f0 variables, allocate memory, and initialize distribution function; simple domain decomposition
  !! the analytic form of the distribution function is velocity normalized to the thermal speed v_th = sqrt(T/m)
  !! energy = 0.5 * [ (v_para-u)^2/alpha + v_perp^2 ] // f(v_para,v_perp) = v_perp*n/[T sqrt(alpha)] exp(-energy)
  !! with a velocity volume element of T/sqrt(2*pi) dv_para dv_perp

  subroutine f_init()

    use mpi
    
    implicit none

    ! real(kind=8), dimension(is_s:is_e,ip_s:ip_e) :: lorentz
    real(kind=8), dimension(tar_s:tar_e,ip_s:ip_e) :: lorentz

    integer :: is, ip, ith, ir, iz, mype1, ierr
    real(kind=8) :: energy, renergy, theta_p,v2
    real(kind=8) :: drift_factor, df0_dpsi

    do mype1=0,sml_totalpe
      call mpi_barrier(sml_comm,ierr)
    enddo

     write(*,*) 'start f_init'

    ! Real space

    ! Velocity grid spacing
    norm_dr = f_smu_max/real(nvrm1,8)
    norm_dz = f_vp_max/real(nvzm1/2,8)

    allocate(norm_vr (ir_s:ir_e),norm_vz (iz_s:iz_e),&
             norm_vr0(ir_s:ir_e),&
            !  vth(ip_s:ip_e, is_s:is_e),&
            ! vel_volume(ir_s:ir_e, ip_s:ip_e, is_s:is_e))

             vth(ip_s:ip_e, tar_s:tar_e),&
             vel_volume(ir_s:ir_e, ip_s:ip_e, tar_s:tar_e))

    
    norm_vr  = 0D0; norm_vz  = 0D0
    norm_vr0 = 0D0;
    vth = 0D0;
    vel_volume = 0D0;   

    do iz = iz_s, iz_e
      norm_vz (iz) = real(iz,8) * norm_dz
    enddo

    do ir = ir_s, ir_e
      norm_vr (ir) = real(ir,8) * norm_dr
      norm_vr0(ir) = real(ir,8) * norm_dr
    enddo
    norm_vr(ir_s) = norm_vr(ir_s) + norm_dr/3D0 
    norm_vr(ir_e) = norm_vr(ir_e) - norm_dr/3D0 

    do is = tar_s, tar_e
      do ip = ip_s, ip_e
        vth(ip,is)=sqrt(eq_t_ev(ip,is)*sml_ev2j/sml_mass(is))
        vel_volume(:,ip,is) = sml_2pi * norm_vr * norm_dr * norm_dz * vth(ip,is)**3D0
      enddo
    enddo
    vel_volume(ir_s,:,:) = vel_volume(ir_s,:,:) / 2D0
    vel_volume(ir_e,:,:) = vel_volume(ir_e,:,:) / 2D0


    !-----------------------------------------------------------------------
    ! START Rotation
    !
    ! define the 'small_n0' regardless of  opt_rotation to calculate phi1
    allocate(small_n0 (ith_s:ith_e, ip_s:ip_e, tar_s:tar_e))

    if (opt_rotation .ne. 0) then

          allocate(omega_rot        (ip_s:ip_e))
          allocate(phi_rot          (ith_s:ith_e, ip_s:ip_e))
          allocate(phi_rot_star     (ith_s:ith_e, ip_s:ip_e))
          allocate(phi_rot_avg      (ip_s:ip_e))  
          allocate(phi_rot_star_avg (ip_s:ip_e))  
          allocate(den_rot          (ip_s:ip_e, tar_s:tar_e))
          allocate(lambda_rot       (ith_s:ith_e, ip_s:ip_e, tar_s:tar_e))
          allocate(dlambda_rot_dth  (ith_s:ith_e, ip_s:ip_e, tar_s:tar_e))


          omega_rot=0D0; ! get from input ??
          phi_rot = 0D0; phi_rot_avg = 0D0;
          phi_rot_star = 0D0; phi_rot_star_avg = 0D0;
          den_rot = 0D0;

          small_n0 = 0D0;

          !tmp
          do ip = ip_s, ip_e
              omega_rot(ip) = omega_rot_in
          enddo
          write(*,*) 'omega_rot', omega_rot
    endif

    

    select case (opt_rotation)
      case(1)
          ! Belli_2009_PPCF eq(13)
          ! only for (e,i) case or (e,i,W) w/ trace limit
          ! isp=0 -> e , isp=1 -> D
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              ! eq(13) case.
              ! However, we can calcualte newton's method
              phi_rot(ith,ip) = omega_rot(ip)**2.0 / 2.0                &
                              * (1.0/vth(ip,1)**2.0-1.0/vth(ip,0)**2.0) &
                              * (geo_R(ip,ith)**2.0 - geo_FSA_R2(ip))   &
                              / (1.0/eq_t_ev(ip,0) + 1.0/eq_t_ev(ip,1)) ! Assume Z=1
          enddo
          enddo

          write(*,*) 'phi_rot', phi_rot

      case(2)
          ! Belli_2009_PPCF eq(43,44)
          ! Calculate the phi_rot_* using 
          ! Refer 'neo_rotation.f90'

          ! Step 1. Calculate the phi_rot_star using Newton method
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              ! refer : neo_rotation.f90
              call cal_phi_rot_star(ip,ith,&
                                    phi_rot_star(ith,ip))
          enddo
          enddo
          write(*,*) 'phi_rot_star', phi_rot_star

          ! Step 2. Calculate the phi_rot_star_avg
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              phi_rot_star_avg(ip) = phi_rot_star_avg(ip) &
                              + phi_rot_star(ith,ip) &
                              * eq_J(ith,ip) * eq_dtheta &
                              / eq_intJ(ip)
          enddo
          enddo
          write(*,*) 'phi_rot_star_avg', phi_rot_star_avg

          ! Step 3. phi_rot = phi_rot_star - phi_rot_star_avg
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              phi_rot(ith,ip) = phi_rot_star(ith,ip) - phi_rot_star_avg(ip)
          enddo
          enddo
          write(*,*) 'phi_rot', phi_rot

          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              phi_rot_avg(ip) = phi_rot_avg(ip) &
                              + phi_rot(ith,ip) &
                              * eq_J(ith,ip) * eq_dtheta &
                              / eq_intJ(ip)
          enddo
          enddo
          write(*,*) 'phi_rot_avg', phi_rot_avg

          

      case(0)
        ! Nothing :: weak rotation
      case default
          write(*,*) 'error :: opt_rotation'
    end select

    write(*,*) 'before smll_n0'
    ! Define the 'smll_n0' to calcualte the phi_1.
    if (opt_rotation .ne. 0) then
          do is = tar_s, tar_e
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              small_n0(ith,ip,is) = eq_den(ip,is) &
                                  * exp(-sml_charge_num(is) / eq_t_ev(ip,is) &
                                        * phi_rot_star(ith,ip)               &
                                        + omega_rot(ip)**2.0 / 2.0           &
                                        * 1.0/vth(ip,is)**2.0                &
                                        * (geo_R(ip,ith)**2.0 - geo_R(ip,ith_s)**2.0) ) ! R(th=0)

          enddo
          enddo
          enddo
    elseif (opt_rotation .eq. 0) then
          write(*,*) 'opt_rotation == 0'

          do is = tar_s, tar_e
          do ip = ip_s, ip_e
          do ith = ith_s, ith_e
              small_n0(ith,ip,is) = eq_den(ip,is)
          enddo
          enddo
          enddo
    endif

    write(*,*) 'after smll_n0'

    if (opt_rotation .ne. 0) then
        do is = tar_s, tar_e
        do ip = ip_s, ip_e
            do ith = ith_s, ith_e
                ! Z*e/T[J] = Z/T[eV]
                lambda_rot(ith,ip,is) = sml_charge_num(is) / eq_t_ev(ip,is) * phi_rot(ith,ip)           &
                                      - omega_rot(ip)**2.0 * geo_R(ip,ith)**2.0 / 2.0 / vth(ip,is)**2.0
            enddo
            
            ! @ theta=0 --> ith_s
            ! eq_den = n0_a(th=0) in rotation.
            den_rot(ip,is) = eq_den(ip,is) * exp( lambda_rot(ith_s,ip,is) )

            call cal_dfdth6(lambda_rot(:,ip,is), dlambda_rot_dth(:,ip,is)) 

        enddo
        enddo

        write(*,*) 'den_rot', den_rot

    endif
    !
    ! END Rotation
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    ! drift component
    !
    allocate(drift_psi   (ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             drift_theta (ith_s:ith_e, ip_s:ip_e, is_s:is_e))   
    allocate(drift_rot_1  (ith_s:ith_e, ip_s:ip_e, is_s:is_e))
    allocate(drift_rot_2  (ith_s:ith_e, ip_s:ip_e, is_s:is_e))


    write(*,*) 'before drift factor'

    do is = is_s, is_e
    do ip = ip_s, ip_e
    do ith = ith_s, ith_e
        ! Jacobian for (psi,theta,phi) not (rho,theta,phi)
        drift_factor = 1D0/eq_omega(ith,ip,is)/eq_B(ith,ip) &
                                *eq_I(ith,ip)*eq_invJ(ith,ip)

        drift_psi  (ith,ip,is) = - drift_factor * eq_dBdth (ith,ip)/eq_B(ith,ip)

        ! Need to know the definition of dBdpsi. psi ?? rho ??
        drift_theta(ith,ip,is) = + drift_factor * eq_dBdpsi(ith,ip)/eq_B(ith,ip)

        if (opt_rotation.ne.0) then
            drift_rot_1  (ith,ip,is) = - drift_factor * dlambda_rot_dth (ith,ip,is)
            drift_rot_2  (ith,ip,is) = 1D0/eq_omega(ith,ip,is) *eq_invJ(ith,ip)      &
                                    * geo_dRdth(ip,ith) * omega_rot(ip) * geo_R(ip,ith)
        endif

    enddo
    enddo
    enddo

    write(*,*) 'after drift factor'
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    ! Allocate memory
    !
    allocate(Ff      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             tmp_Ff  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),& !youknow
             g1      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             k1      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             k2      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             k3      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             k4      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             g1_k1   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),& ! youknow
             g1_k2   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             g1_k3   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             g1_k4   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             f1_tmp  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             f1      (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             f1_an   (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             f1_neo  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             f0_M    (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             f0_M_tar(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, tar_s:tar_e),& ! youknow
             g1_tar  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, tar_s:tar_e),& ! youknow : opt_restart=3
             f1_tar  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, tar_s:tar_e),& ! youknow : opt_restart=3
             
            !  col_power(ip_s:ip_e, ith_s:ith_e, is_s:is_e, tar_s:tar_e),& ! youknow
             col_power(ith_s:ith_e, ip_s:ip_e, is_s:is_e, tar_s:tar_e),& ! JPL reorder 
             col_dens (ith_s:ith_e, ip_s:ip_e, is_s:is_e, tar_s:tar_e),& ! JPL reorder 
             col_inner_info (ith_s:ith_e, ip_s:ip_e, is_s:is_e, tar_s:tar_e, 0:2, 1:col_iteration),& ! JPL reorder 

             ! TORIC youknow for half grid
             qld_B (ip_s:ip_e, ith_s:ith_e, 1:nvrm1, 1:nvzm1),  &
             qld_C (ip_s:ip_e, ith_s:ith_e, 1:nvrm1, 1:nvzm1),  &
             qld_F (ip_s:ip_e, ith_s:ith_e, 1:nvrm1, 1:nvzm1),  &
             num_n   (ith_s:ith_e,ip_s:ip_e,is_s:is_e), num_T(ith_s:ith_e,ip_s:ip_e,is_s:is_e),&
             num_p   (ith_s:ith_e,ip_s:ip_e,is_s:is_e), &


             FSAn(ip_s:ip_e,is_s:is_e),FSAT(ip_s:ip_e,is_s:is_e),FSAp(ip_s:ip_e,is_s:is_e),&
             dfdt_c  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             dfdt_q  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             dfdt_p  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             dfdt_imp(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             dfdt_d  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&

             dfdt_e  (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&

             frv     (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e))

    g1=0D0; f1=0D0; f1_an=0D0; f1_neo=0D0; f0_M=0D0; f0_M_tar=0D0; g1_tar=0D0; f1_tar=0D0
    k1=0D0; k2=0d0; k3=0d0; k4=0d0
    f1_tmp=0D0;
    g1_k1=0D0; g1_k2=0D0; g1_k3=0D0; g1_k4=0D0
    num_n=0D0; num_T=0D0; num_p=0D0;  dfdt_e=0D0
    qld_B=0D0; qld_C=0D0; qld_F=0D0;
    dfdt_imp=0d0
    size_f0 = nvz * nvr * nth * nspecies
    dfdt_c=0D0; dfdt_q=0D0; dfdt_p=0d0; dfdt_d=0d0
    FSAn=0.0;FSAT=0.0;FSAp=0.0;
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    ! Define the init f0_M & f0_drift considering rotation
    !
    ! define the target Maxwellian for collision module
    do is = tar_s, tar_e ! 0,f_isp
    do ip = ip_s, ip_e
    do ith = ith_s, ith_e
    do ir = ir_s, ir_e
    do iz = iz_s, iz_e

        energy    = 0.5D0*((norm_vz(iz)-sml_flow(is))**2D0/f_velo_aniso + norm_vr0(ir)**2D0)

        if (opt_rotation .eq. 0) then
            f0_M_tar(iz,ir,ith,ip,is) = eq_den(ip,is)*exp(-energy) ! It doesn't matter if unnecessary repetition occurs.
        else
            ! rotation on
            energy = energy + lambda_rot(ith,ip,is)
            f0_M_tar(iz,ir,ith,ip,is) = den_rot(ip,is)*exp(-energy) ! It doesn't matter if unnecessary repetition occurs.
        endif

    enddo ! do iz
    enddo ! do ir
    enddo ! do ith
    enddo ! do ip
    enddo ! do is
    
    do is = is_s, is_e
    do ip = ip_s, ip_e
      theta_p=(vth(ip,is)/sml_cvel)**2D0
      ! write (*,*) 'is', 'ip', is, ip, eq_den(ip,is)
      do ith = ith_s, ith_e
      do ir = ir_s, ir_e
      do iz = iz_s, iz_e
        if (eq_t_ev(ip,is).gt. 51097) then  ! corresponds to 1e5 eV
            !Relativistic version
            v2 = norm_vr0(ir)*norm_vr0(ir)+norm_vz(iz)*norm_vz(iz)
            renergy = v2*vth(ip,is)**2D0
            renergy = sqrt(1+renergy/sml_cvel**2D0) !relativistic energy
            lorentz(is,ip) = renergy ! lorentz factor
            !===Maxwellian=====!
            Ff  (iz,ir,ith,ip,is) = eq_den(ip,is) * exp(-renergy/theta_p)/(4.0*3.14159*theta_p*1.6248) ! Since Bessel function is not accurate, coefficient is multiplied.!Value 1.6248 corresponds to the value of Bessel function when T=0.511Mev 
            f0_M(iz,ir,ith,ip,is) = eq_den(ip,is) * exp(-renergy/theta_p)/(4.0*3.14159*theta_p*1.6248) 
        else
            v2 = norm_vr0(ir)*norm_vr0(ir)+norm_vz(iz)*norm_vz(iz)
            energy    = 0.5D0*((norm_vz(iz)-sml_flow(is))**2D0/f_velo_aniso + norm_vr0(ir)**2D0)

            if (opt_rotation .eq. 0) then
                Ff      (iz,ir,ith,ip,is) = eq_den(ip,is)*exp(-energy)
                f0_M    (iz,ir,ith,ip,is) = eq_den(ip,is)*exp(-energy)
                df0_dpsi =  (eq_dlnn(ip,is) + (0.5*v2-1.5)*eq_dlnT(ip,is) ) &
                          / (eq_Bp(ith,ip)*geo_R(ip,ith)) & ! d/dpsi = 1/(dpsi/dr) * d/dr, dpsi/dr = grad psi = Bp*R
                          * Ff(iz,ir,ith,ip,is)
            elseif (opt_rotation .ne. 0) then

                ! NEED TO CONSIDER 'dR/dr * dPhi_0/dr' in d(ln N)/dr
                ! temporary set these values are zero
                energy = energy + lambda_rot(ith,ip,is)
                Ff      (iz,ir,ith,ip,is) = den_rot(ip,is)*exp(-energy)
                f0_M    (iz,ir,ith,ip,is) = den_rot(ip,is)*exp(-energy)
                df0_dpsi =  (eq_dlnn(ip,is) + (0.5*v2-1.5+lambda_rot(ith,ip,is))*eq_dlnT(ip,is) ) &
                          / (eq_Bp(ith,ip)*geo_R(ip,ith)) & ! d/dpsi = 1/(dpsi/dr) * d/dr, dpsi/dr = grad psi = Bp*R
                          * Ff(iz,ir,ith,ip,is)
            endif
            
            ! only for w/o rotation ... ?
            f1_an(iz,ir,ith,ip,is) = - eq_I(ith,ip) * vth(ip,is)*norm_vz(iz) &
                                    /eq_omega(ith,ip,is) &
                                    * df0_dpsi

            if (opt_0M.ne.0) then
                if (opt_rotation .eq. 0) then
                    ! dfdt_d is divided by conv_factor
                    dfdt_d(iz,ir,ith,ip,is) = - drift_psi(ith,ip,is) & ! (-) sign due to drift from LHS
                                              * vth(ip,is)*vth(ip,is) &
                                              * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir)/2D0) &
                                              * df0_dpsi
                elseif (opt_rotation .ne. 0) then
                    dfdt_d(iz,ir,ith,ip,is) = + drift_psi(ith,ip,is)   & 
                                              * vth(ip,is)*vth(ip,is)  &
                                              * (norm_vz(iz)*norm_vz(iz)+norm_vr0(ir)*norm_vr0(ir)/2D0) &
                                              + drift_rot_1(ith,ip,is) &
                                              * vth(ip,is)*vth(ip,is)  &
                                              + drift_rot_2(ith,ip,is) &
                                              * vth(ip,is) * norm_vz(iz) * 2D0
                    dfdt_d(iz,ir,ith,ip,is) = - dfdt_d(iz,ir,ith,ip,is) & ! (-) sign due to drift from LHS
                                              * df0_dpsi
                endif

                if (opt_time_solver.eq.-1) then
                  dfdt_d(iz,ir,ith,ip,is) =  dfdt_d(iz,ir,ith,ip,is) * sml_dt * eq_tau(ip,opt_time_species,is) !FP4D_v1.1_steady
                else
                  dfdt_d(iz,ir,ith,ip,is) =  dfdt_d(iz,ir,ith,ip,is) * sml_dt
                endif
                                          
            endif



        endif
      enddo ! do iz
      enddo ! do ir
      enddo ! do ith
    enddo ! do ip
    enddo ! do is
    !-----------------------------------------------------------------------

    return

  end subroutine f_init

  subroutine f1_Maxwellian(f1_in, f1_out)
    implicit none
    real(kind=8), intent(in)  :: f1_in (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)
    real(kind=8), intent(out) :: f1_out(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e)
    real(kind=8) :: int_den (ith_s:ith_e, ip_s:ip_e, is_s:is_e)
    real(kind=8) :: FSA_den (ip_s:ip_e, is_s:is_e)
    real(kind=8) :: energy, conv_factor
    integer :: is, ip, ith, ir, iz


    int_den = 0D0
    FSA_den = 0D0
    do is = is_s, is_e
    do ip = ip_s, ip_e
      conv_factor = vth(ip,is)*sqrt(6.2831853071795862D0)
      conv_factor = 1D0/conv_factor**3
      do ith = ith_s, ith_e
        do ir=ir_s,ir_e
        do iz=iz_s,iz_e
              int_den(ith,ip,is) = int_den(ith,ip,is) + f1_in(iz,ir,ith,ip,is) * vel_volume(ir,ip,is) * conv_factor
        enddo ! do iz
        enddo ! do ir
      enddo ! ith
      FSA_den(ip,is) = compute_FSA(int_den(:,ip,is),ip)
    enddo ! do ip
    enddo ! do is

    write(*,*) 'f1_FSA_den', FSA_den

    do is = is_s, is_e
    do ip = ip_s, ip_e
      do ith = ith_s, ith_e
      do ir = ir_s, ir_e
      do iz = iz_s, iz_e
          energy    = 0.5D0*((norm_vz(iz)-sml_flow(is))**2D0/f_velo_aniso + norm_vr0(ir)**2D0)

          if (opt_rotation .eq. 0) then
              f1_out (iz,ir,ith,ip,is) = f1_in(iz,ir,ith,ip,is) - FSA_den(ip,is)*exp(-energy)
          elseif (opt_rotation .ne. 0) then
              ! NEED TO CONSIDER 'dR/dr * dPhi_0/dr' in d(ln N)/dr
              ! temporary set these values are zero
              energy = energy + lambda_rot(ith,ip,is)
              f1_out (iz,ir,ith,ip,is) = f1_in(iz,ir,ith,ip,is) - FSA_den(ip,is)*exp(-energy)*exp( lambda_rot(ith_s,ip,is) )
          endif
          
      enddo ! do iz
      enddo ! do ir
      enddo ! do ith
    enddo ! do ip
    enddo ! do is
  end subroutine f1_Maxwellian

  subroutine stsl_init(total_pe,nthreads,mype,source_den,source_t_ev,source_width)

    use mpi
    implicit none
    integer, intent(in) :: total_pe, nthreads, mype
    real (kind=8) :: conv_factor ! youknow

    integer :: is, ip, ith, ir, iz, mype1, ierr
    real (kind=8) :: energy
    real (kind=8), dimension(is_s:is_e), intent(in) :: source_den, source_t_ev
    real (kind=8) :: source_energy, source_width,  loss_energy!HJL

    
    ! Allocate memory
    allocate(df_sts (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             df_stl (iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             sts_int(ith_s:ith_e,ip_s:ip_e,is_s:is_e))
    df_sts=0d0; df_stl=0d0;sts_int=0d0;

    if ((opt_NB.eq.1).or.(opt_APS.eq.1)) then
      if (NB_sp.eq.APS_sp) then
        if (mype.eq.0) then
          write(*,*) 'NB_sp and APS_sp are same. Check the input!'
        end if
      end if
      ! Steady-State Source
      
      do is = is_s, is_e
      do ip = ip_s, ip_e

        conv_factor = 1.D0/sqrt((sml_2pi*eq_t_ev(ip,is)*sml_ev2j / sml_mass(is))**3.0) !youknow  

        do ith=ith_s, ith_e
        do ir = ir_s, ir_e
        do iz=iz_s,iz_e
          source_energy = sqrt(norm_vz(iz)**2D0 + norm_vr0(ir)**2D0) &
                         -sqrt(2D0*source_t_ev(is)/eq_t_ev(ip,is))
          source_energy=0.5D0 * source_energy**2D0
          if (is.eq.APS_sp) then
            df_sts(iz,ir,ith,ip,is) = df_sts(iz,ir,ith,ip,is)+ source_den(is)*exp (-source_energy)*sml_dt
          endif
          if (is.eq.NB_sp) then
            ! 'sml_dt' will be considered in fp4d.F90
            ! df_sts(iz,ir,ith,ip,is) = df_sts(iz,ir,ith,ip,is)+frv(iz,ir,ith,ip)
            ! df_sts(iz,ir,ith,ip,is) = df_sts(iz,ir,ith,ip,is)+frv(iz,ir,ith,ip)*1D-4
            ! df_sts(iz,ir,ith,ip,is) = df_sts(iz,ir,ith,ip,is)+frv(iz,ir,ith,ip)*1D-1 
            if (opt_NB_aux .eq. 1) then ! default
              df_sts(iz,ir,ith,ip,is) = df_sts(iz,ir,ith,ip,is)+frv(iz,ir,ith,ip)*NB_pwrscale
            elseif (opt_NB_aux .eq. 2) then
              ! NuBDeC source(weight) unit : [#/s], FP4D source unit : [#/(m)^3/(m/s)^3/s]
              ! Also, df_sts is calcualted @ fp4d.F90
              ! In fp4d.F90, the distribution functions are normalized by conv_factor

              df_sts(iz,ir,ith,ip,is) = df_sts(iz,ir,ith,ip,is)+frv(iz,ir,ith,ip)*NB_pwrscale &
                                                                  /conv_factor &
                                                                  /eq_Volume(ith,ip)    &
                                                                  /vel_volume(ir,ip,is)
              ! if ( (ip .eq. 9) .and. (ith .eq. 2) .and. (ir .eq. 48) .and. (iz .eq. 50)) then
              !   write(*,*) 'conv_factor', conv_factor
              !   write(*,*) 'eq_Volume', eq_Volume(ith,ip)
              !   write(*,*) 'vel_volume' ,vel_volume(ir,ip,is)
              ! endif

            endif
          endif

          sts_int(ith,ip,is) = sts_int(ith,ip,is) + df_sts(iz,ir,ith,ip,is) &
                                                      * conv_factor &
                                                      * vel_volume(ir,ip,is)
          ! 근데 이게 문제가 conv_factor 나눠져서 들어가는게 정상이라. 고려해줘야 할 듯.?? NBI 일 때랑 alpha 일 때랑 따로 고려해야할듯.
                                                      
        enddo
        enddo
        enddo
      enddo
      enddo
      write(*,*) 'df_sts', df_sts(50, 48, 2, 9, NB_sp)
      write(*,*) 'frv',       frv(50, 48, 2, 9)

    !Steady-State Loss
      do is = is_s, is_e
      do ip = ip_s, ip_e
      do ith = ith_s, ith_e
      do ir = ir_s, ir_e
      do iz = iz_s, iz_e
          loss_energy = 0.5D0*(norm_vz(iz)**2D0+ norm_vr0(ir)**2D0)
          if (is.eq.APS_sp) then
            df_stl(iz,ir,ith,ip,is) = df_stl(iz,ir,ith,ip,is)+source_den(is)*exp(-loss_energy)
          endif

          if (is.eq.NB_sp) then
            ! youknow : NB_sink_opt
            ! NB_sink_opt=0 --> original
            ! NB_sink_opt=1,2,3 --> set the loss term in first species (Background Deuterium)
            select case(NB_sink_opt)
              case(0) ! original
                  df_stl(iz,ir,ith,ip,is)  = sts_int(ith,ip,is)*exp(-loss_energy)
              case(1) !!! NB_sp = 1 but loss_sp = 0
                  df_stl(iz,ith,ir,ip, NB_sink_sp) = 1.2D-3 * sts_int(ith,ip,is)*exp(-1D-2*loss_energy)
              case(2)
                  df_stl(iz,ith,ir,ip, NB_sink_sp) = 0.6D-3 * sts_int(ith,ip,is)*exp(-0.5D-2*loss_energy)
              case(3)
                  df_stl(iz,ith,ir,ip, NB_sink_sp) = 2.4D-3 * sts_int(ith,ip,is)*exp(-2D-2*loss_energy)
              case(-1)
                  exit
              case default
                  write(*,*) 'error-NB_sink_opt'
              end select
          endif

      enddo
      enddo
      enddo
      enddo
      enddo
    endif
    return

  end subroutine stsl_init

  subroutine loss_check(check_ok,loss_dt)
  implicit none
  integer :: check_ok, sum_ok
  integer :: is, ip, ith, ir, iz
  real (kind=8) :: loss_dt
  sum_ok=0;
  do is=is_s, is_e
  do ip=ip_s, ip_e
  do ir=ir_s, ir_e
  do ith=ith_s, ith_e
  do iz=iz_s, iz_e
    if (Ff(iz,ir,ith,ip,is) .ge. df_stl(iz,ir,ith,ip,is)*loss_dt) then
      sum_ok=sum_ok+1;
    endif
  enddo
  enddo
  enddo
  enddo
  enddo
  ! write(*,*) sum_ok, size_f0
  if (sum_ok.eq.size_f0) then
    check_ok=1
  else
    check_ok=0
  endif
  write(*,*) sum_ok, size_f0, check_ok
  return

  end subroutine loss_check


  subroutine cal_phi_rot_star(ip,ith,x)
    integer, intent(in) :: ip, ith
    integer :: is, n, nmax
    real(kind=8) :: numer, denom    
    real(kind=8) :: x0, fac
    real(kind=8), intent(out) :: x ! phi_rot_star

    ! refer : neo_rotation.f90
    ! aniso 는 또 따로 고려해줘야 함... lam_rot_aniso 참조
    ! ae_flag

    ! initial guess for Phi_rot_star(psi,theta)
    x = 0.05

    nmax = 200

    if (ith.eq.ith_s) then
      write(*,*) 'sml_charge_num', sml_charge_num
      write(*,*) 'eq_den', eq_den
      write(*,*) 'eq_t_ev', eq_t_ev
      write(*,*) 'omega_rot', omega_rot
      write(*,*) 'vth', vth
      write(*,*) 'geo_R', geo_R

   endif

    do n=1,nmax
        ! use Newton's method to solve the quasi-neutrality relation        
        numer = 0D0
        denom = 0D0
        do is = tar_s, tar_e ! 이거 헷갈리네. is_s로 할지, is_e로 할지.
                             ! 우선 general species 없이모두 포함이니...
          fac = -sml_charge_num(is) / eq_t_ev(ip,is) * x  &
              + omega_rot(ip)**2.0 / 2.0                  &
              * 1.0/vth(ip,is)**2.0                       &
              * (geo_R(ip,ith)**2.0 - geo_R(ip,ith_s)**2.0)    ! R(th=0)

          fac = sml_charge_num(is) * eq_den(ip,is) * exp(fac)
               
          numer = numer + fac
          ! d(n0a)/d(phi_star)=(-Z/T[eV])*n0a
          denom = denom - sml_charge_num(is) / eq_t_ev(ip,is) * fac 
        enddo

        x0 = x
        x  = x0 - numer/denom
        
        if (abs(x-x0) < 1.0e-12) exit
        
        if(n > nmax) exit
        
    enddo

    if(n > nmax) then
        write(*,*) 'error :: phi_rot_star failed to converge @ ith',ith
        return
    endif

  end subroutine cal_phi_rot_star



  !> Deallocate f0 memory
  subroutine f_finalize
    implicit none

    deallocate(vth,norm_vr,norm_vz)
    deallocate(norm_vr0)

    deallocate(drift_psi,drift_theta)
    deallocate(drift_rot_1,drift_rot_2)

    deallocate(Ff,dfdt_c,dfdt_q,dfdt_p,dfdt_imp,f1,k1,k2,k3,k4,&
               f1_an,f1_neo,dfdt_d,df_sts,sts_int,df_stl,dfdt_e,frv,g1,f0_M,f0_M_tar,g1_tar,f1_tar,&
               tmp_Ff,qld_B,qld_C,qld_F,& !youknow
               g1_k1,g1_k2,g1_k3,g1_k4,&
               f1_tmp,&
               num_n,num_T,num_p,FSAn,FSAT,FSAp)

    ! define the 'small_n0' regardless of  opt_rotation to calculate phi1
    deallocate(small_n0)

    if (opt_rotation .ne. 0) then
          deallocate(omega_rot)
          deallocate(phi_rot)
          deallocate(phi_rot_star)
          deallocate(phi_rot_avg)  
          deallocate(phi_rot_star_avg)  
          deallocate(den_rot)
          deallocate(lambda_rot)
          deallocate(dlambda_rot_dth)

    endif

    return

  end subroutine f_finalize

end module f_module
