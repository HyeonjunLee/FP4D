module eqmodule
  
  use iso_c_binding
  use FP4D_globals
  implicit none
  ! ********************************************************

  ! !!!! inputs !!!!!! --> defined @ FP4D_globals.F90
  ! public :: xx0,xx1,q0,q1,q2,q3,q4,a0,R0,B0,ase_q0,ase_a0,ase_R0,ase_B0
  ! public :: ase_kappa,ase_del,ase_shift
  ! public :: npsi,nth,ifd0,ifd1,op_geo_info
  ! public :: ipe ! youknow

  ! integer :: npsi,nth,ifd0,ifd1,op_geo_info
  ! integer :: ipe ! youknow
  ! real *8 :: xx0,xx1,q0,q1,q2,q3,q4,a0,R0,B0,ase_q0,ase_a0,ase_R0,ase_B0
  ! real *8 :: ase_kappa,ase_del,ase_shift

  !!! Outputs !!!
  public :: R_arr,Z_arr,th_arr,psi_arr,dpdr_arr,dpdz_arr,dpdrr_arr,dpdzz_arr,dpdrz_arr,B_arr,Bp_arr,I_arr,J_arr,gradpsi2,gradth2,gradpsith,dBdp_arr,dBdt_arr,temp_arr,temp2_arr
   
  real *8,dimension(:,:),allocatable :: R_arr,Z_arr,th_arr,psi_arr,dpdr_arr,dpdz_arr,dpdrr_arr,dpdzz_arr,dpdrz_arr,B_arr,Bp_arr,I_arr,J_arr,gradpsi2,gradth2,gradpsith,dBdp_arr,dBdt_arr,temp_arr,temp2_arr


      interface 
        subroutine init_diag(ptr,length) bind(C)
        use iso_c_binding
        real(c_double),intent(out),dimension(*) :: ptr
        integer(c_int), value :: length
!         c_ptr :: R_arr
!         R_arr = c_null_ptr 
        end subroutine

        subroutine farr_to_carr(intpar,dblpar) bind(C)
        use iso_c_binding
        integer(c_int),intent(in),dimension(*) :: intpar
        real(c_double),intent(in),dimension(*) :: dblpar
        end subroutine

      end interface 

contains

  ! ------------------------------------------------------------------------
 
  subroutine set_eq

  
   implicit none

   logical :: file_exists
   character (20) :: filein
   real(c_double),dimension(:),allocatable,target :: dblin
   integer(c_int),dimension(:),allocatable,target :: intin
   integer(c_int) :: istat



  write(*,*) 'youknow-R0', R0


  call alloc1

! fortran array to carray 
  allocate(intin(5) )
  allocate(dblin(17))

  intin(1)=npsi
  intin(2)=nth
  intin(3)=ifd0
  intin(4)=ifd1
  intin(5)=op_geo_info
  intin(6)=ipe

  dblin(1)=xx0
  dblin(2)=xx1
  dblin(3)=q0
  dblin(4)=q1
  dblin(5)=q2
  dblin(6)=q3
  dblin(7)=q4
  dblin(8)=a0
  dblin(9)=R0
  dblin(10)=B0
  dblin(11)=ase_q0
  dblin(12)=ase_a0
  dblin(13)=ase_R0
  dblin(14)=ase_B0
  dblin(15)=ase_kappa
  dblin(16)=ase_del
  dblin(17)=ase_shift

  ! fortran array to c array
  call farr_to_carr(intin,dblin) 
  
  ! main c functions
  call set_my_constants()  
  call reset_flag0() 
  call set_equilibrium_geometry()
  call set_rz_simulation_domain()
  call setup_xx_theta_grid()
  
  ! c array to fortran array
  call carr_to_farr
 
  end subroutine set_eq


  subroutine carr_to_farr
      real(c_double),dimension(:),allocatable,target :: Output
      integer(c_int) :: istat,n1,n2,n3
      integer(c_int) :: rows,cols
      integer :: i,j,idx,nwr,set_eq

      rows=npsi+1
      cols=nth

      ! idx of carray starts from 0
      allocate( Output(0:rows*cols*20-1),stat=istat )
      if ( istat /= 0 ) then
         ! error handling elided; error stop for this simple example
         error stop
      end if

      write(*,*) 'Output Index Check', LBOUND(Output, 1), UBOUND(Output, 1)

      !Output = 0.0_c_double


      call init_diag(Output,size(Output))

      nwr=(npsi+1)*nth

      do i = 1, npsi + 1
          do j = 0, nth-1
              idx = (i - 1) * nth + j

              ! write(*,*) 'idx check',idx,i,j

              R_arr(i, j)       = Output(idx)
              Z_arr(i, j)       = Output(idx + 1 * nwr)
              th_arr(i, j)      = Output(idx + 2 * nwr)
              psi_arr(i, j)     = Output(idx + 3 * nwr)
              dpdr_arr(i, j)    = Output(idx + 4 * nwr)
              dpdz_arr(i, j)    = Output(idx + 5 * nwr)
              dpdrr_arr(i, j)   = Output(idx + 6 * nwr)
              dpdzz_arr(i, j)   = Output(idx + 7 * nwr)
              dpdrz_arr(i, j)   = Output(idx + 8 * nwr)
              B_arr(i, j)       = Output(idx + 9 * nwr)
              Bp_arr(i, j)      = Output(idx + 10 * nwr)
              I_arr(i, j)       = Output(idx + 11 * nwr)
              J_arr(i, j)       = Output(idx + 12 * nwr)
              gradpsi2(i, j)    = Output(idx + 13 * nwr)
              gradth2(i, j)     = Output(idx + 14 * nwr)
              gradpsith(i, j)   = Output(idx + 15 * nwr)
              dBdp_arr(i, j)    = Output(idx + 16 * nwr)
              dBdt_arr(i, j)    = Output(idx + 17 * nwr)
              temp_arr(i, j)    = Output(idx + 18 * nwr)
              temp2_arr(i, j)   = Output(idx + 19 * nwr)
          end do
      end do

      ! write(*,*) 'whole check', Output
      ! write(*,*) 'temp2_arr', temp2_arr


!      print *, 'R_arr',R_arr
!      print *, 'temp2_arr',temp2_arr

      return

   end subroutine


  subroutine alloc1
   implicit none
   allocate(R_arr(1:npsi+1,0:nth-1),Z_arr(1:npsi+1,0:nth-1),th_arr(1:npsi+1,0:nth-1),psi_arr(1:npsi+1,0:nth-1),dpdr_arr(1:npsi+1,0:nth-1),dpdz_arr(1:npsi+1,0:nth-1),dpdrr_arr(1:npsi+1,0:nth-1),dpdzz_arr(1:npsi+1,0:nth-1),dpdrz_arr(1:npsi+1,0:nth-1),B_arr(1:npsi+1,0:nth-1),Bp_arr(1:npsi+1,0:nth-1),I_arr(1:npsi+1,0:nth-1),J_arr(1:npsi+1,0:nth-1),gradpsi2(1:npsi+1,0:nth-1),gradth2(1:npsi+1,0:nth-1),gradpsith(1:npsi+1,0:nth-1),dBdp_arr(1:npsi+1,0:nth-1),dBdt_arr(1:npsi+1,0:nth-1),temp_arr(1:npsi+1,0:nth-1),temp2_arr(1:npsi+1,0:nth-1))
 
  end subroutine alloc1

  subroutine deallocate_eqmodule
  implicit none
  !deallocate (global_arrays) 

      write(*,*) ''
      write(*,*) 'eqmodule :: DEALLOCATE'
      write(*,*) ''

      deallocate(R_arr, Z_arr, th_arr, psi_arr, dpdr_arr, dpdz_arr, dpdrr_arr, dpdzz_arr, dpdrz_arr, &
              B_arr, Bp_arr, I_arr, J_arr, gradpsi2, gradth2, gradpsith, dBdp_arr, dBdt_arr, &
              temp_arr, temp2_arr)

  end subroutine deallocate_eqmodule

end module eqmodule



