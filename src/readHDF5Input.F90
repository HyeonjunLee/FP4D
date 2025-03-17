module readHDF5Input
use FP4D_globals
use HDF5
use mpi

implicit none

  integer :: n_particle

  interface readVariable
    module procedure readVariable_integer
    module procedure readVariable_scalar
    module procedure readVariable_1d
    module procedure readVariable_2d
    module procedure readVariable_3d
    module procedure readVariable_5d
    module procedure readVariable_6d
  end interface readVariable

  integer(HID_T), private :: fileID, groupID, maingroupID


contains
  subroutine setup_hdf(NB_num_in)
    implicit none
    integer :: NB_num_in

    n_particle = NB_num_in

    write(*,*) 'NB_num_in', NB_num_in

  end subroutine setup_hdf

  subroutine checkPrevData(nname)
    implicit none
    character(len=16) :: pname
    character(len=16),intent(out)::nname
    print *, 'Restart Option is activated !'
    print *, 'You must not change any input variable except a few time-relevant variables!'
    print *, 'We are finding the last data ...'
    !Check the dataXX.h5
    call checkDataName(pname,nname)
    !Open the last dataNN.h5
    call openInputFile(trim(pname),'/')
    print *, 'after  openInputFile'
    !Read and Import
    call readPrevData()
    print *, 'after  readPrevData'
    !Close the data
    call closeInputFile()
    print *, 'after  closeInputFile'
    return
  end subroutine checkPrevData

  subroutine checkDataName(pastname,newname)
    implicit none
    integer :: dnumber
    character(len=16), intent(out) :: newname,pastname
    character(len=3) :: x1
    logical :: there
    newname='data00.h5'
    do dnumber=0,99
      pastname=newname
      write(x1,'(I0.1)') dnumber
      if (dnumber.lT.10) then
        newname='data0'//trim(x1)//'.h5'
      else
        newname='data'//trim(x1)//'.h5'
      endif
      inquire(file=trim(newname),exist=there)
      if (there==.true.) then
      elseif ((dnumber.ne.0).and.(there==.false.)) then
        print *, trim(pastname), ' is the last data file!'
        print *, trim(newname), ' will be started!' 
        exit
      elseif ((dnumber.eq.0).and.(there==.false.)) then
        print *, 'There is no last data! Please set opt_restart=0 and try again!'
        stop
      endif
    enddo

    return
  end subroutine checkDataName

  subroutine openInputFile(fileName, groupName)
    implicit none
    character(len=*), intent(in) :: fileName, groupName
    character(len=100) :: group
    integer :: HDF5Error
    
    call h5open_f(HDF5Error)
       if (HDF5Error < 0) then
          print *,"Error initializing HDF5."
          stop
       end if

    ! Open input file
    call h5fopen_f(trim(fileName), H5F_ACC_RDONLY_F, fileID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening file:", fileName
      stop
    end if
       
    ! Open group
    call h5gopen_f(fileID, trim(groupName), groupID, HDF5Error)  !ym
    if (HDF5Error < 0) then                                          !ym
      print *,"Error opening group:",groupName                       !ym
      stop
    end if
  
  end subroutine openInputFile

  subroutine readPrevData()                                        !!!!HJ!!
    use FP4D_globals, only : opt_restart
    use FP4D_globals, only : f1_neo

    implicit none
    real(kind=8), dimension(:,:,:,:,:,:), allocatable :: g1_load, f1_load,Ff_load
    real(kind=8), dimension(:,:,:,:,:), allocatable :: df_sts_load                                 
    real(kind=8), dimension(:,:,:), allocatable :: phi1_load
    integer :: ndstep_load, is, ip, ith, ir, iz

    call readVariable(ndstep_load,'ndstep')
    allocate(Ff_load(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep_load), &
             g1_load(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep_load),&
             f1_load(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e, 0:ndstep_load),&
             df_sts_load(iz_s:iz_e, ir_s:ir_e, ith_s:ith_e, ip_s:ip_e, is_s:is_e),&
             phi1_load(ith_s:ith_e, ip_s:ip_e, 0:ndstep_load))
    call readVariable(Ff_load,'Ff')
    call readVariable(g1_load,'g1')
    call readVariable(f1_load,'f1')
    call readVariable(phi1_load,'phi1')
    call readVariable(df_sts_load,'df_sts')
    call readVariable(sml_time,'end_time')

    select case (opt_restart)
      case (1) ! origin
        do is=is_s,is_e
        do ip=ip_s,ip_e
        do ith=ith_s,ith_e
        do ir=ir_s,ir_e
        do iz=iz_s,iz_e
            Ff(iz,ir,ith,ip,is)=Ff_load(iz,ir,ith,ip,is,ndstep_load)
            g1(iz,ir,ith,ip,is)=g1_load(iz,ir,ith,ip,is,ndstep_load)
            f1(iz,ir,ith,ip,is)=f1_load(iz,ir,ith,ip,is,ndstep_load)
            f1_neo(iz,ir,ith,ip,is)=f1_load(iz,ir,ith,ip,is,ndstep_load)
            df_sts(iz,ir,ith,ip,is)=df_sts_load(iz,ir,ith,ip,is) 
        enddo
        enddo
        enddo
        enddo
        enddo
        do ip=ip_s,ip_e
        do ith=ith_s,ith_e
            eq_phi1(ith,ip)=phi1_load(ith,ip,ndstep_load)
        enddo
        enddo
      case (2) ! only call the f1 for opt_QLD_dist .eq. 4
        do is=is_s,is_e
        do ip=ip_s,ip_e
        do ith=ith_s,ith_e
        do ir=ir_s,ir_e
        do iz=iz_s,iz_e
            f1_neo(iz,ir,ith,ip,is)=f1_load(iz,ir,ith,ip,is,ndstep_load)
        enddo
        enddo
        enddo
        enddo
        enddo
      case default
        write(*,*) 'error :: opt_restart'
    end select
    deallocate(Ff_load,g1_load,f1_load,df_sts_load,phi1_load)

  end subroutine
 
  subroutine readHDF5()                                        !!!!ymstart!!!!!
    implicit none                                                   
    real(kind=8), dimension(:,:), allocatable :: variable             
    real(kind=8), dimension(1:12,1:n_particle) :: Beam_deposition          
    real(kind=8), dimension(1:n_particle) :: NB_R, Phi,NB_Z, Vpar,Vper, r_process,NB_theta,psi_n,idx_psi
    real(kind=8), dimension(1:n_particle) :: Vpar_n, Vper_n, weight
    integer, dimension(1:n_particle) :: idx_th
    real(kind=8) :: NB_vth
    real(kind=8) :: maxR, minR, maxZ, minZ, Vpar_max, Vper_max, percent
    real(kind=8) :: r_zero, min_d, d
    integer :: var_R, var_Z, idx_Vpar,idx_Vper,m,ith,ip,min_ith,min_ip
    integer :: NVpar,NVper,NVparm1,NVperm1
    character(len=30) :: dataname 
    integer :: n_theta, n_psi, theta_num, psi_num	!!!!!
    real(kind=8), dimension(:,:), allocatable :: raw_data, res	!!!!!!

   
    allocate(variable(1:n_particle,1:12))

    dataname = 'Beam_Deposition'
    call readVariable(variable,dataname)
    Beam_deposition = reshape(variable,[12,n_particle])

    write(*,*) 'nparticles' , n_particle

 
    NB_R = Beam_deposition(1,:)
    Phi = Beam_deposition(2,:)
    NB_Z = Beam_deposition(3,:)
    Vpar = Beam_deposition(4,:)   !Parallel velocity
    Vper = Beam_deposition(5,:)   !Perpendicular velocity
    weight = Beam_deposition(8,:)
   !NB_theta calculation
  write(*,*)'Beammmmmmmmm',weight(1),weight(2)
    r_zero = geo_R(1,(nth)/4)  ! youknow : ??? 
    r_process = NB_R - r_zero

   do m =1,n_particle
     
    if ((NB_R(m).GT.r_zero).AND.(NB_Z(m).GT.0.0))then
     NB_theta(m) = atan(NB_Z(m)/r_process(m))
    elseif ((NB_R(m).LT.r_zero).AND.(NB_Z(m).GT.0.0))then
     NB_theta(m) = atan(NB_Z(m)/r_process(m)) + sml_pi
    elseif ((NB_R(m).LT.r_zero).AND.(NB_Z(m).LT.0.0))then
     NB_theta(m) = atan(NB_Z(m)/r_process(m)) + sml_pi
    elseif ((NB_R(m).GT.r_zero).AND.(NB_Z(m).LT.0.0))then
     NB_theta(m) = atan(NB_Z(m)/r_process(m)) + sml_pi*2.0
    endif
 
   enddo

  !idx_psi calculation
      
   do m = 1,n_particle
    idx_th(m) = int(floor(NB_theta(m)/( eq_dtheta )))
   enddo

  
   do m = 1,n_particle
    min_d = 10000.0
    if (idx_th(m).eq.(nth-1)) then
      do ith = ith_s,ith_e
        do ip = ip_s, ip_e
          d = (geo_R(ip,ith) - NB_R(m))**2 + (geo_Z(ip,ith) - NB_Z(m))**2
          if(d.LT.min_d) then
           min_d = d
           min_ith = ith
           min_ip  = ip
          endif
        enddo
      enddo

    else
      do ith = idx_th(m), idx_th(m)+1, 1
        do ip = ip_s, ip_e
          d = (geo_R(ip,ith) - NB_R(m))**2 + (geo_Z(ip,ith) - NB_Z(m))**2
          if(d.LT.min_d) then
          min_d = d
          min_ith = ith
          min_ip  = ip
         endif
        enddo
      enddo
    endif
   ! idx_psi(m) = geo_psi(min_ip,min_ith)
   idx_psi(m)= min_ip
  enddo
  
  
  !!!!!!!!!
  do m = 1, n_particle

    NB_vth = vth(idx_psi(m),NB_sp)
    Vpar_n(m) = Vpar(m)/NB_vth   !Parallel velocity
    Vper_n(m) = Vper(m)/NB_vth   !Parallel velocity
  enddo


! write(*,*) 'vth', vth
! write(*,*) 'eq_t_ev', eq_t_ev

! write(*,*) 'NB_vth', NB_vth

! write(*,*) 'Vpar_n', Vpar_n
  !grid
  
  Vpar_max = ceiling(maxval(Vpar_n))
  Vper_max = ceiling(maxval(Vper_n))
  
  write(*,*) 'vpar_min', minval(Vpar_n)
  write(*,*) 'vpar_max', maxval(Vpar_n)

  write(*,*) 'vperp_min', minval(Vper_n)
  write(*,*) 'vperp_max', maxval(Vper_n)

  ! if (f_vp_max .le. )

  ! write(*,*) 'vpar_min',minval(Vpar_n)         
  ! write(*,*) 'vper_max',Vper_max         

  NVpar = nvzm1 !2*f_nvp 
  NVper = nvrm1 !f_nmu 

  NVparm1 = NVpar + 1 
  NVperm1 = NVper + 1 
  
 !make f from nubdec
 do m = 1,n_particle
 idx_Vpar = int(floor((Vpar_n(m)+Vpar_max)*NVpar/(Vpar_max)/2.0)-nvzm1/2)
 idx_Vper = int(floor((Vper_n(m))*NVper/(Vper_max)) )
 frv( idx_Vpar, idx_Vper, idx_th(m),idx_psi(m)) = frv( idx_Vpar, idx_Vper, idx_th(m),idx_psi(m))  + weight(m)

 end do
 
  end subroutine readHDF5                                 !!!!!!!ymend!!!!


  subroutine closeInputFile()

    integer :: HDF5Error

    ! Close the group
    call h5gclose_f(groupID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing input group"
      stop
    end if

    ! Close the file
    call h5fclose_f(fileID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing input file"
      stop
    end if

  end subroutine closeInputFile

  subroutine readVariable_integer(variable,varname)

    integer, intent(inout) :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: HDF5Error

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_INTEGER, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_integer

  subroutine readVariable_scalar(variable,varname)

    real(kind=8), intent(inout) :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(0) :: dimensions
    integer :: HDF5Error

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_scalar

  subroutine readVariable_1d(variable,varname)

    real(kind=8), intent(inout), dimension(:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(1) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_1d

  subroutine readVariable_2d(variable,varname)

    real(kind=8), intent(inout), dimension(:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(2) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    write(*,*) 'shape(variable)', shape(variable)
    write(*,*) 'dimensions', dimensions


    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,  dimensions
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_2d

  subroutine readVariable_3d(variable,varname)

    real(kind=8), intent(inout), dimension(:,:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(3) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,  dimensions
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine readVariable_3d

  subroutine readVariable_5d(variable,varname)

    real(kind=8), intent(inout), dimension(:,:,:,:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(5) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,  dimensions
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine 

  subroutine readVariable_6d(variable,varname)

    real(kind=8), intent(inout), dimension(:,:,:,:,:,:), allocatable :: variable
    character(len=*), intent(in) :: varname
    integer(HID_T) :: dataSetID
    integer(HSIZE_T), dimension(6) :: dimensions
    integer :: HDF5Error

    if (.not. allocated(variable)) then
      print *,"Tried to read into unallocated array:",varname
      stop
    end if

    dimensions = shape(variable)

    ! Open Dataset
    call h5dopen_f(groupID, varname, dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error opening dataset for variable:",varname
      stop
    end if
    
    ! Read variable
    call h5dread_f(dataSetID, H5T_NATIVE_DOUBLE, variable, dimensions, HDF5Error)
    if (HDF5Error < 0) then
      print *,  dimensions
      print *,"Error reading variable:",varname
      stop
    end if

    ! Close Dataset
    call h5dclose_f(dataSetID, HDF5Error)
    if (HDF5Error < 0) then
      print *,"Error closing dataset for variable:",varname
      stop
    end if

  end subroutine 
end module readHDF5Input

