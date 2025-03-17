module readNetCDFInput

use netcdf
use FP4D_globals

implicit none

integer, parameter :: onedim=1, twodim=2, threedim=3, fourdim=4, fivedim=5

integer :: tstart(onedim),tcount(onedim)
integer :: xtstart(twodim),xtcount(twodim), spptstart(threedim), spptcount(threedim)
integer :: xdim, tdim, spdim, pdim
integer :: npid, densid, massid

integer :: f_fdim(fivedim)
integer :: f_smu_maxdim, f_vp_maxdim
integer :: f_nvpsdim, f_nodesdim, f_nmudim, f_nspsdim
integer :: f_ftstart(fivedim), f_ftcount(fivedim)
integer :: f_fid, f_df0gid, df_stsid, df_stspid, df_stsrfid, df_stlid

character (len = *), parameter :: UNITS = "units"
character (len = *), parameter :: Ex_UNITS = ""
character (len = *), parameter :: Vx_UNITS = ""
character (len = *), parameter :: Pos_UNITS = ""

interface rank
   module procedure rank_real
   module procedure rank_integer
   module procedure rank_character
   module procedure rank_scalar
   module procedure rank_1d
   module procedure rank_1d_integer
   !module procedure rank_1d_nonalloc
   module procedure rank_2d
   module procedure rank_3d
   module procedure rank_4d
   module procedure rank_5d
end interface rank

contains

! -----------------------------------------------------------------------------------
! youknow_TORIC
  subroutine openNCDF()

    implicit none
    integer :: status
    integer :: ncid_B, ncid_C, ncid_F
    integer :: Qldce_LD4DID2_B, Qldce_LD4DID2_C, Qldce_LD4DID2_F
    integer, parameter :: file_name_len = 80
    character (file_name_len) :: file_name_B, file_name_C, file_name_F
! qld_B,C,F are defined in f_module.F90

! QLE_coeff-B
    file_name_B = 'toric_qlde_B.cdf'

    status = nf90_open(file_name_B, NF90_NOWRITE, ncid_B)

    status = nf90_inq_varid(ncid_B, "Qldce_LD4D_B", Qldce_LD4DID2_B)

    status = nf90_get_var(ncid_B, Qldce_LD4DID2_B, qld_B)

    write(*,*) 'READ_TORIC-NETCDF_TEST'
    write(*,*) 'lbound', lbound(qld_B,1), lbound(qld_B,2),lbound(qld_B,3),lbound(qld_B,4)
    write(*,*) 'ubound', ubound(qld_B,1), ubound(qld_B,2),ubound(qld_B,3),ubound(qld_B,4)
    status = nf90_close(ncid_B)

! QLE_coeff-C
    file_name_C = 'toric_qlde_C.cdf'

    status = nf90_open(file_name_C, NF90_NOWRITE, ncid_C)

    status = nf90_inq_varid(ncid_C, "Qldce_LD4D_C", Qldce_LD4DID2_C)

    status = nf90_get_var(ncid_C, Qldce_LD4DID2_C, qld_C)

    status = nf90_close(ncid_C)

! QLE_coeff-F
    file_name_F = 'toric_qlde_F.cdf'

    status = nf90_open(file_name_F, NF90_NOWRITE, ncid_F)

    status = nf90_inq_varid(ncid_F, "Qldce_LD4D_F", Qldce_LD4DID2_F)

    status = nf90_get_var(ncid_F, Qldce_LD4DID2_F, qld_F)

    status = nf90_close(ncid_F)
  end subroutine openNCDF

  subroutine closeInputNCDF()

    implicit none


  end subroutine closeInputNCDF


  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  integer function rank_real(A)
    real (kind=8), intent(in) :: A
    rank_real=size(shape(A))
    return
  end function rank_real

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
    real, intent(in) :: A
    rank_scalar=size(shape(A))
    return
  end function rank_scalar

  integer function rank_1d(A)
    real, dimension(:), allocatable, intent(in) :: A
    rank_1d=size(shape(A))
    return
  end function rank_1d

  integer function rank_1d_integer(A)
    integer, dimension(:), allocatable, intent(in) :: A
    rank_1d_integer=size(shape(A))
    return
  end function rank_1d_integer

  
  ! conflicts with rank_1d
  !integer function rank_1d_nonalloc(A)
  !  PetscScalar, dimension(:), intent(in) :: A
  !  rank_1d_nonalloc=size(shape(A))
  !  return
  !end function rank_1d_nonalloc

  integer function rank_2d(A)
    real, dimension(:,:), allocatable, intent(in) :: A
    rank_2d=size(shape(A))
    return
  end function rank_2d

  integer function rank_3d(A)
    real, dimension(:,:,:), allocatable, intent(in) :: A
    rank_3d=size(shape(A))
    return
  end function rank_3d

  integer function rank_4d(A)
    real, dimension(:,:,:,:), allocatable, intent(in) :: A
    rank_4d=size(shape(A))
    return
  end function rank_4d


  integer function rank_5d(A)
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: A
    rank_5d=size(shape(A))
    return
  end function rank_5d


end module readNetCDFInput
