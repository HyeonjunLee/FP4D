module FP4D_init

  use FP4D_globals
  use FP4D_timer_lib

  contains

  subroutine FP4D_check()
    implicit none
    real (kind=8) :: QN_sum, ion_sum
    integer :: isp


    !-----------------------------------------------------------------------
    ! Check miscellaneous conditions
    if (opt_logaritm_NEO) then
      if (sml_charge_num(0) .gt. 0D0 ) then
        STOP 'Stopping due to the first element is not an electron'
      endif
    endif
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Check the Quasi-Neutrality given input
    if (sml_mype==0) then
        QN_sum = 0D0
        ion_sum = 0D0

        do isp=0, sml_nsp
          QN_sum = QN_sum + sml_den(isp)*sml_charge_num(isp)
          if ( sml_charge_num(isp) .gt. 0D0) then
            ion_sum = ion_sum + sml_den(isp)*sml_charge_num(isp)
          endif
        enddo

        if (QN_sum .eq. 0D0) then
          write(*,*) 'Quasi-Neutrality is satisfied'
        else
          write(*,*) 'ERROR :: Quasi-Neutrality is not satisfied'
          write(*,*) 'Current  electron dens :', ion_sum - QN_sum
          write(*,*) 'Expected electron dens :', ion_sum
          STOP 'Stopping due to quasi-neutrality'
        endif
    endif
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Check the time-stepping variables
    if (step_Col .eq. 1) then
        start_step_Col = sml_nstep + 1
    else
      if (mod(start_step_Col, step_Col) .ne. 0) then
        STOP 'Stopping due to wrong relation [start_step_Col] with [step_Col]'
      endif
      if (mod(sml_nstep, step_Col) .ne. 0) then
        STOP 'Stopping due to wrong relation [sml_nstep] with [step_Col]'
      endif
    endif

    if (on_RHS .eq. cycle_RHS) then
      write(*,*) 'RHS calculate every time-step'
    else
      write(*,*) 'RHS do not calculate every time-step'
      if (start_step_RHS .lt. start_step_Col) then
        STOP 'Stopping :: [start_step_RHS] have to be eqaul of bigger than [start_step_Col]'
      endif
    endif

    if (mod(start_step_RHS, step_Col).ne.0) then
      STOP 'Stopping due to wrong [start_step_RHS]'
    endif

    if (mod(sml_nstep, out_step) .ne. 0) then
      STOP 'Stopping due to wrong [out_step]'
    endif
    !-----------------------------------------------------------------------

  end subroutine FP4D_check


  subroutine FP4D_timer_init()
    implicit none

    ! w/o mype == 0
    call timer_lib_init('TOTAL')
    call timer_lib_init('RHS')
    call timer_lib_init('Col')
    call timer_lib_init('Remove_homogeneous')
    call timer_lib_init('Potent')
    call timer_lib_init('hdf5')
    call timer_lib_init('Post')

  end subroutine FP4D_timer_init


  subroutine FP4D_runfile_create()
    implicit none
    integer :: unit
    character(len=20) :: timestamp
    character(len=100) :: hostname
    !-----------------------------------------------------------------------
    ! out.fp4d.run
    if (sml_mype==0) then
        open(unit=io_FP4Dout, file=runfile_FP4Dout, status='replace')
        close(io_FP4Dout)

        call getenv('FP4D_SYSTEM', hostname)

        if (trim(hostname) .eq. '') then
          call getenv('HOSTNAME', hostname)
        endif

        call cal_timestamp(timestamp)

        open(unit=io_FP4Dout, file=runfile_FP4Dout, status='old', position='append')     
          write(io_FP4Dout, '(A)') "-------------------------------------"
          write(io_FP4Dout, '(A)') "RUNNING ON : " // trim(hostname)
          write(io_FP4Dout, '(A)') "-------------------------------------"
          write(io_FP4Dout, '(A)') "START Time: " // trim(timestamp)
          write(io_FP4Dout, '(A)') "-------------------------------------"
        close(io_FP4Dout)
    endif

    !-----------------------------------------------------------------------    
  end subroutine FP4D_runfile_create

  subroutine FP4D_output_open()
    implicit none
    integer :: unit
    character (len=200) :: filename
    integer :: isp
    !-----------------------------------------------------------------------
    ! OPEN output_FSA*_**.dat
    if (sml_mype==0) then

        ! Prepare the output file
        ! Only rank 0 writes because all grid points are identical anyway.
        ! At each time step, the density, parallel and perpendicular and total temperature as well as the entropy are printed out for physics verification.
        do isp=sml_isp,sml_nsp
            write(filename,'(a,i0.2)') 'out.fp4d.FSAn_',isp
            unit=100+isp
            open(unit=100+isp,file=filename,status="replace",action="write")
            write(unit,'(a,i4)') '# FP4D collision output, species ',isp
            write(unit,'(a,e12.2,a,e12.2)') '# Mass: ',sml_mass(isp),' Charge:',sml_charge(isp)
        enddo

        do isp=sml_isp,sml_nsp
            write(filename,'(a,i0.2)') 'out.fp4d.FSAT_',isp
            unit=200+isp
            open(unit=200+isp,file=filename,status="replace",action="write")
            write(unit,'(a,i4)') '# FP4D collision output, species ',isp
            write(unit,'(a,e12.2,a,e12.2)') '# Mass: ',sml_mass(isp),' Charge:',sml_charge(isp)
        enddo
        do isp=sml_isp,sml_nsp
            write(filename,'(a,i0.2)') 'out.fp4d.FSAp_',isp
            unit=300+isp
            open(unit=300+isp,file=filename,status="replace",action="write")
            write(unit,'(a,i4)') '# FP4D collision output, species ',isp
            write(unit,'(a,e12.2,a,e12.2)') '# Mass: ',sml_mass(isp),'Charge:',sml_charge(isp)
        enddo
    endif
    !-----------------------------------------------------------------------
  end subroutine FP4D_output_open

  subroutine FP4D_output_close()
    implicit none
    integer :: unit
    integer :: isp
    !-----------------------------------------------------------------------
    ! CLOSE output_FSA*_**.dat
    if (sml_mype==0) then
        do isp=sml_isp,sml_nsp
        close(100+isp)
        close(200+isp)
        close(300+isp)
        enddo
    endif
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Record the elapsed time to out.fp4d.run
    !
    if (sml_mype==0) then

        open(unit=io_FP4Dout, file=runfile_FP4Dout, status='old', position='append')     
          write(io_FP4Dout, '(A)') "-------------------------------------"
          write(io_FP4Dout, '(A)') "ELAPSED TIME (s) "
          write(io_FP4Dout, '(A)') "-------------------------------------"
          write(io_FP4Dout, '(A,F12.6)') "TOTAL      : ", timer_lib_time('TOTAL')
          write(io_FP4Dout, '(A,F12.6)') "RHS        : ", timer_lib_time('RHS')
          write(io_FP4Dout, '(A,F12.6)') "Col        : ", timer_lib_time('Col')
          write(io_FP4Dout, '(A,F12.6)') "Potent     : ", timer_lib_time('Potent')
          write(io_FP4Dout, '(A,F12.6)') "hdf5       : ", timer_lib_time('hdf5')
          write(io_FP4Dout, '(A,F12.6)') "Post       : ", timer_lib_time('Post')
        close(io_FP4Dout)
    endif
    !-----------------------------------------------------------------------   

  end subroutine FP4D_output_close


  subroutine cal_timestamp(timestamp)
    implicit none
    integer :: time_vals(8)
    character(len=20) :: timestamp

    call DATE_AND_TIME(values=time_vals)

    ! define the timestamp as follows
    write(timestamp, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
        time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7)

  end subroutine cal_timestamp


  subroutine print_msg(msg)
    implicit none
    character(len=*) :: msg

    print *, "!-------------------------------------"
    print *, "! " // trim(msg)
    print *, "!-------------------------------------"

  end subroutine print_msg

end module FP4D_init