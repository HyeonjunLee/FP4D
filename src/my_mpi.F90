module my_mpi_module

  contains
  !> XGC common MPI init routine.
  !! This function will call MPI_init and set the total number of
  !! processes and the rank IDs.
  subroutine MY_MPI_INIT
    use FP4D_globals
    use mpi
    implicit none
    integer :: ierror

    call mpi_init(ierror)
    call mpi_comm_rank(MPI_COMM_WORLD,sml_mype,ierror)
    call mpi_comm_size(mpi_comm_world,sml_totalpe,ierror)
    call mpi_comm_dup(mpi_comm_world,sml_comm,ierror)
    sml_comm_null = mpi_comm_null

    write(*,*) 'total mpi,myid',sml_totalpe,sml_mype
    allocate(istart1(0:sml_totalpe-1),iend1(0:sml_totalpe-1))
    allocate(istart2(0:sml_totalpe-1),iend2(0:sml_totalpe-1))
  end subroutine MY_MPI_INIT


  subroutine get_MPI_grid
    use FP4D_globals
    use mpi
    implicit none
    integer :: ierror, i, j, temp1, temp2
    integer :: npsi_part,npsi_rem,nth_part,nth_rem

    np1=sml_totalpe !/col_f_nthreads       
    np2=np1/npsi
    if (np1>0) then
      if (np2>0) then !2D mpi paralel on psi and theta
            np1=npsi
            nptot_mpi=np1*np2
            allocate(ip1_mpi(nptot_mpi))
            nth_part=nth/np2
            nth_rem=nth-np2*nth_part
            temp1=1
            do i=1, npsi
              istart1((i-1)*np2:i*np2-1)=temp1
                    temp1=temp1+1
              iend1((i-1)*np2:i*np2-1)=temp1-1
              temp2=1
              do j=1, np2
                  ip1_mpi((i-1)*np2+(j-1)+1)=((i-1)*np2+(j-1))   !mype starting from 0
                  istart2(((i-1)*np2+(j-1)):((i-1)*np2+j)-1)=temp2
                  if (j<nth_rem+1) then
                    temp2=temp2+nth_part+1
                  else
                    temp2=temp2+nth_part
                  endif
                  iend2(((i-1)*np2+(j-1)):((i-1)*np2+j)-1)=temp2-1
              enddo
            enddo
      else  ! 1D mpis parallel on psi 
            nptot_mpi=np1
            allocate(ip1_mpi(nptot_mpi)) 
            istart2(:)=1 
            iend2(:)=nth
            temp1=1  
            npsi_part=npsi/np1
            npsi_rem=npsi-npsi_part*np1 ! remaininng from np1 
            do i=1, np1
              ip1_mpi(i)=(i-1)+1
              istart1(i-1)=temp1
              if (i<npsi_rem+1) then
                    temp1=temp1+npsi_part+1
              else
                    temp1=temp1+npsi_part
              endif
              iend1(i-1)=temp1-1
            enddo
      endif
    else 
      nptot_mpi=0
      istart1(:)=1
      iend1(:)=npsi
      istart2(:)=1
      iend2(:)=nth
    end if  

    if (sml_mype.eq.0) write(*,*) 'MPI partition on psi and theta',np1, np2
    if (sml_mype.eq.0) write(*,*) 'istart1', istart1
    if (sml_mype.eq.0) write(*,*) 'iend1', iend1
    if (sml_mype.eq.0) write(*,*) 'istart2', istart2
    if (sml_mype.eq.0) write(*,*) 'iend2', iend2
    if (sml_mype.eq.0) write(*,*) 'ip1_mpi', ip1_mpi

  end subroutine get_MPI_grid

  !> Finalize MPI
  subroutine MY_MPI_FINALIZE
    use FP4D_globals
    use mpi
    implicit none
    integer :: ierr

    call mpi_barrier(sml_comm,ierr)
    ! print *,'ierr barrier',ierr
    call mpi_finalize(ierr)
    write(*,*) 'Process', sml_mype, 'finalized with ierr=', ierr

  end subroutine MY_MPI_FINALIZE

end module my_mpi_module
