	module blacs_mod
!	----------------------------
!	interface to BLACS  routines
!	----------------------------


        implicit none
      
        public :: init_blacs, exit_blacs
        public :: icontext, iam, ip_start, ip_end
        public :: isroot, rsrc, csrc
        public :: nprow, npcol, myrow, mycol, nprocs
        public :: parallel_partition
        public :: global_collect, global_sum

        public ::global_sum_scalar_r8, global_sum_scalar_z16, &
        global_sum_vector_r8,  global_sum_vector_z16, &
        global_sum_matrix_r8,  global_sum_matrix_z16

        public ::  global_collect_r8,    global_collect_z16, &
        global_collect_2d_r8, global_collect_2d_z16,&
        global_collect_3d_r8, global_collect_3d_z16


        integer :: icontext=-1,myrow, mycol, nprow,npcol
        integer :: rsrc=0, csrc=0
        logical :: isroot
        integer :: nprocs, iam
        integer, dimension(:), allocatable :: ip_start, ip_end !


        private

        interface

	subroutine zgesd2d(icontext, m,n,A,lda,rdest,cdest)
        implicit none
	integer icontext,m,n,lda,rdest,cdest
	complex  A(lda,*)
	end subroutine zgesd2d

	subroutine dgesd2d(icontext, m,n,A,lda,rdest,cdest)
        implicit none
	integer icontext,m,n,lda,rdest,cdest
	real  A(lda,*)
	end subroutine dgesd2d

	subroutine zgerv2d(icontext, m,n,A,lda,rsrc,csrc)
        implicit none
	integer icontext,m,n,lda,rsrc,csrc
	complex  A(lda,*)
	end subroutine zgerv2d

	subroutine dgerv2d(icontext, m,n,A,lda,rsrc,csrc)
        implicit none
	integer icontext,m,n,lda,rsrc,csrc
	real  A(lda,*)
	end subroutine dgerv2d

	subroutine zgsum2d(icontext,scope,top,m,n,A,lda,rdest,cdest)
        implicit none
	integer icontext,m,n,lda,rdest,cdest
	complex  A(*)
	character scope
	character top
	end subroutine zgsum2d

	subroutine dgsum2d(icontext,scope,top,m,n,A,lda,rdest,cdest)
        implicit none
	integer icontext,m,n,lda,rdest,cdest
	real  A(*)
	character scope
	character top
	end subroutine dgsum2d



	subroutine blacs_pcoord(icontext,pnum,prow,pcol)
        implicit none
	integer icontext,pnum,prow,pcol
	end subroutine blacs_pcoord

	subroutine blacs_pinfo(mypnum,nprocs)
        implicit none
	integer mypnum,nprocs
	end subroutine blacs_pinfo

	subroutine blacs_get(icontext,what,val)
        implicit none
	integer icontext,what,val
	end subroutine blacs_get


	subroutine blacs_gridinit(icontext,order,nprow,npcol)
        implicit none
	integer icontext,nprow,npcol
	character order
	end subroutine blacs_gridinit


	subroutine blacs_gridinfo(icontext,nprow,npcol,myprow,mypcol)
        implicit none
	integer icontext,nprow,npcol,myprow,mypcol
	end subroutine blacs_gridinfo

	end interface

        interface global_sum
        module procedure global_sum_scalar_r8,  &
        global_sum_scalar_z16,  global_sum_vector_r8, &  
        global_sum_vector_z16,  global_sum_matrix_r8, & 
        global_sum_matrix_z16
        end interface


        interface global_collect
        module procedure  global_collect_r8, global_collect_z16, &
        global_collect_2d_r8, global_collect_2d_z16, &
        global_collect_3d_r8, global_collect_3d_z16
        end interface

        !------------------------------
        contains



        subroutine init_blacs(icontext_inout)

        integer , intent(inout) :: icontext_inout

        logical :: isok
        integer :: istart

        icontext = icontext_inout

        call blacs_pinfo(iam,nprocs)
        isok = (nprocs > 0)
        istart = nint(sqrt(real(nprocs)))
        if (nprocs.gt.1) then
           istart = min( nprocs-1,istart)
        endif
        do npcol = istart,1,-1
          nprow = nprocs/npcol
          if (nprow*npcol .eq. nprocs) exit
        enddo
        isok = (nprow*npcol .eq. nprocs)


! CALL BLACS_GET (icontxt, what, val) 
! what=0 means  icontext is ignored, system default context is returned 
! in val.
        call blacs_get(-1,0,icontext)
        call blacs_gridinit( icontext, 'C',nprow,npcol)

        call blacs_gridinfo( icontext, nprow,npcol, myrow,mycol)


        allocate(ip_start(0:nprocs-1),ip_end(0:nprocs-1))
        ip_start=0
        ip_end=0

        isroot = (myrow.eq.rsrc).and.(mycol.eq.csrc)

        return

        end subroutine
 

        subroutine exit_blacs(icontinue_in)

        integer, intent(in),optional :: icontinue_in
        integer :: icontinue

        call  blacs_barrier( icontext, 'A')

        icontinue = 0
        if (present(icontinue_in)) then
          icontinue = icontinue_in
        endif

        deallocate(ip_start,ip_end)
        call blacs_exit( icontinue )
        return
        end subroutine exit_blacs


!       =========
!       parallel_partition
!       =========
        subroutine parallel_partition( nsize, istart, iend )

        integer, intent(in) :: nsize
        integer, dimension(0:(nprocs-1)), intent(out) :: istart,iend

        integer, dimension(0:(nprocs-1)) :: isize
        integer :: iremain,iproc
        logical :: isok

        iremain = mod(nsize,nprocs)
        isize(0:(nprocs-1)) = int(nsize/nprocs)
        if (iremain .gt. 0) then
          isize(0:(iremain-1)) = isize(0:(iremain-1)) + 1
        endif

!       ------------
!       double check
!       ------------
        if (sum(isize(0:(nprocs-1))).ne.nsize) then
           write(*,*) 'nsize ',nsize
           write(*,*) 'isize(:) ', isize(:)
           stop '** error in parallel_partition ** '
        endif

        istart(0) = 1
        iend(0) = isize(0)
        do iproc=1,(nprocs-1)
           istart(iproc) = iend(iproc-1) + 1
           iend(iproc) = istart(iproc) + isize(iproc)-1
        enddo

!       ------------
!       double check
!       ------------
        isok = (iend(nprocs-1).eq.nsize)
        if (.not.isok) then
           write(*,*) 'nsize ', nsize
           write(*,*) 'istart(:) ', istart(:)
           write(*,*) 'end(:) ', iend(:)
           write(*,*) 'isize(:) ',isize(:)
           stop '** error in partitioin ** '
        endif
        
        do iproc=0,(nprocs-1)
          isok = (iend(iproc)-istart(iproc)+1) .eq. isize(iproc)
          if (.not.isok) then
             write(*,*) 'iproc,istart,iend,isize ', &
              iproc,istart(iproc),iend(iproc),isize(iproc)
             stop '** error in parallel_partition ** '
          endif
        enddo


        return
        end subroutine parallel_partition



!       ==========
!       global sum
!       ==========


        subroutine global_sum_scalar_r8(total)

        real , intent(inout) :: total
        real , dimension(1) :: tarray

        tarray(1) = total
        call dgsum2d( icontext, 'A', ' ', 1,1,tarray,1,-1,-1)
        total = tarray(1)

        return
        end subroutine global_sum_scalar_r8


        subroutine global_sum_scalar_z16(total)

        complex , intent(inout) :: total
        complex , dimension(1) :: tarray

        tarray(1) = total
        call zgsum2d( icontext, 'A', ' ', 1,1,tarray,1,-1,-1)
        total = tarray(1)

        return
        end subroutine global_sum_scalar_z16

        subroutine global_sum_vector_r8(total)

        real , dimension(:), intent(inout) :: total

        integer :: m,n,lda

        m = size(total,1)
        n = 1
        lda = m

        call dgsum2d( icontext,'A', ' ', &              
              m,n,total,lda,            -1,-1)

        return
        end subroutine global_sum_vector_r8

        subroutine global_sum_vector_z16(total)

        complex , dimension(:), intent(inout) :: total

        integer :: m,n,lda

        m = size(total,1)
        n = 1
        lda = m

        call zgsum2d( icontext,'A', ' ', &   
              m,n,total,lda,            -1,-1)

        return
        end subroutine global_sum_vector_z16

        subroutine global_sum_matrix_r8(total)

        real , dimension(:,:), intent(inout) :: total

        integer :: m,n,lda

        m = size(total,1)
        n = size(total,2)
        lda = m

        call dgsum2d( icontext,'A',' ',   & 
             m,n,total,lda,          -1,-1 )

        return
        end subroutine global_sum_matrix_r8

        subroutine global_sum_matrix_z16(total)

        complex , dimension(:,:) :: total

        integer :: m,n,lda

        m = size(total,1)
        n = size(total,2)
        lda = m

        call zgsum2d( icontext,'A',' ', &
              m,n,total,lda,   -1,-1 )

        return
        end subroutine global_sum_matrix_z16

!       ---------------
!       global_collect
!       ---------------

        subroutine global_collect_r8(vec)

        real , dimension(:),intent(inout) :: vec

        integer, dimension(0:(nprocs-1)) :: istart,iend
        real , dimension(:), allocatable :: tmp
        integer :: isize,lb,ub,  iistart,iiend

        call parallel_partition(size(vec),istart,iend)

        lb = lbound(vec,1)
        ub = ubound(vec,1)

        allocate( tmp(lb:ub) )
        tmp(:) = 0.0

        iistart = (lb-1) + istart(iam)
        iiend = (lb-1) + iend(iam)
        isize = iiend-iistart + 1

        if (isize.ge.1) then
          tmp( iistart:iiend ) = vec( iistart:iiend )
        endif

        call global_sum( tmp )
        vec(:) = tmp(:)

        deallocate(tmp)
        return
        end subroutine global_collect_r8

        

        subroutine global_collect_z16(vec)

        complex , dimension(:),intent(inout) :: vec

        integer, dimension(0:(nprocs-1)) :: istart,iend
        complex , dimension(:), allocatable :: tmp
        integer :: isize,lb,ub,   iistart,iiend

        isize = size(vec)
        call parallel_partition(isize,istart,iend)

        lb = lbound(vec,1)
        ub = ubound(vec,1)

        allocate( tmp(lb:ub) )
        tmp(:) = 0.0

        iistart = (lb-1) + istart(iam)
        iiend = (lb-1) + iend(iam)

        isize = iiend-iistart+1
        if (isize.ge.1) then
          tmp( iistart:iiend ) = vec( iistart:iiend )
        endif
        call global_sum( tmp )
        vec(:) = tmp(:)

        deallocate(tmp)
        return
        end subroutine global_collect_z16




        subroutine global_collect_2d_1_r8( mat, ldim )

!       ----------------------------
!       collect one vector at a time
!       ----------------------------
        real , dimension(:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        real , dimension(:), allocatable :: tmp
        integer :: i,  lb,ub
        
        lb = lbound(mat,ldim)
        ub = ubound(mat,ldim)
        allocate( tmp(lb:ub) )

        select case (ldim)
        case (1)
          do i=lbound(mat,2),ubound(mat,2)
              tmp(lb:ub) = mat(lb:ub,i)
              call global_collect( tmp )
              mat(lb:ub,i) = tmp(lb:ub)
          enddo
        case (2)
          do i=lbound(mat,1),ubound(mat,1)
              tmp(lb:ub) = mat(i,lb:ub)
              call global_collect( tmp )
              mat(i,lb:ub) = tmp(lb:ub)
          enddo
        case default
          write(*,*) 'global_collect: invalid ldim ',ldim
          stop '** error in global_collect ** '
        end select

        deallocate( tmp )
       
        return
        end subroutine global_collect_2d_1_r8




        subroutine global_collect_2d_1_z16( mat, ldim )

!       ----------------------------
!       collect one vector at a time
!       ----------------------------

        complex , dimension(:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        complex , dimension(:), allocatable :: tmp
        integer :: i,  lb,ub
        
        lb = lbound(mat,ldim)
        ub = ubound(mat,ldim)
        allocate( tmp(lb:ub) )

        select case (ldim)
        case (1)
          do i=lbound(mat,2),ubound(mat,2)
              tmp(lb:ub) = mat(lb:ub,i)
              call global_collect( tmp )
              mat(lb:ub,i) = tmp(lb:ub)
          enddo
        case (2)
          do i=lbound(mat,1),ubound(mat,1)
              tmp(lb:ub) = mat(i,lb:ub)
              call global_collect( tmp )
              mat(i,lb:ub) = tmp(lb:ub)
          enddo
        case default
          write(*,*) 'global_collect: invalid ldim ',ldim
          stop '** error in global_collect ** '
        end select

        deallocate( tmp )
       
        return
        end subroutine global_collect_2d_1_z16




        subroutine global_collect_3d_1_z16( mat, ldim )

!       -----------------------------------------
!       perform collection one vector at the time
!       -----------------------------------------

        complex , dimension(:,:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        complex , dimension(:), allocatable :: tmp
        integer :: i,j,    lb,ub
        
        lb = lbound(mat,ldim)
        ub = ubound(mat,ldim)
        allocate( tmp(lb:ub) )

        select case (ldim)
        case (1)

          do j=lbound(mat,3),ubound(mat,3)
          do i=lbound(mat,2),ubound(mat,2)
              tmp(lb:ub) = mat(lb:ub,i,j)
              call global_collect( tmp )
              mat(lb:ub,i,j) = tmp(lb:ub)
          enddo
          enddo

        case (2)

          do j=lbound(mat,3),ubound(mat,3)
          do i=lbound(mat,1),ubound(mat,1)
              tmp(lb:ub) = mat(i,lb:ub,j)
              call global_collect( tmp )
              mat(i,lb:ub,j) = tmp(lb:ub)
          enddo
          enddo

        case (3)

          do j=lbound(mat,2),ubound(mat,2)
          do i=lbound(mat,1),ubound(mat,1)
             tmp(lb:ub) = mat(i,j,lb:ub)
             call global_collect( tmp )
             mat(i,j,lb:ub) = tmp(lb:ub)
          enddo
          enddo

        case default
          write(*,*) 'global_collect: invalid ldim ',ldim
          stop '** error in global_collect ** '
        end select
       
        deallocate( tmp )

        return
        end subroutine global_collect_3d_1_z16




        subroutine global_collect_3d_1_r8( mat, ldim )

!       -----------------------------------------
!       perform collection one vector at the time
!       -----------------------------------------

        real , dimension(:,:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        real , dimension(:), allocatable :: tmp
        integer :: i,j,    lb,ub
        
        lb = lbound(mat,ldim)
        ub = ubound(mat,ldim)
        allocate( tmp(lb:ub) )

        select case (ldim)
        case (1)

          do j=lbound(mat,3),ubound(mat,3)
          do i=lbound(mat,2),ubound(mat,2)
              tmp(lb:ub) = mat(lb:ub,i,j)
              call global_collect( tmp )
              mat(lb:ub,i,j) = tmp(lb:ub)
          enddo
          enddo

        case (2)

          do j=lbound(mat,3),ubound(mat,3)
          do i=lbound(mat,1),ubound(mat,1)
              tmp(lb:ub) = mat(i,lb:ub,j)
              call global_collect( tmp )
              mat(i,lb:ub,j) = tmp(lb:ub)
          enddo
          enddo

        case (3)

          do j=lbound(mat,2),ubound(mat,2)
          do i=lbound(mat,1),ubound(mat,1)
             tmp(lb:ub) = mat(i,j,lb:ub)
             call global_collect( tmp )
             mat(i,j,lb:ub) = tmp(lb:ub)
          enddo
          enddo

        case default
          write(*,*) 'global_collect: invalid ldim ',ldim
          stop '** error in global_collect ** '
        end select

        deallocate( tmp )
       
        return
        end subroutine global_collect_3d_1_r8




        subroutine global_collect_3d_r8( mat, ldim )

!       -----------------------------------------
!       perform collection one matrix  at the time
!       -----------------------------------------

        real , dimension(:,:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        real , dimension(:,:), allocatable :: tmp
        integer :: j

        select case (ldim)
        case (1,2)
           allocate( tmp(size(mat,1),size(mat,2)) )
           do j=lbound(mat,3),ubound(mat,3)
              tmp(:,:) = mat(:,:,j)
              call global_collect( tmp, ldim )
              mat(:,:,j) = tmp(:,:)
           enddo
        case (3)
           allocate( tmp(size(mat,1),size(mat,3)) )
           do j=lbound(mat,2),ubound(mat,2)
               tmp(:,:)  = mat(:,j,:)
               call global_collect( tmp, 2 )
               mat(:,j,:) = tmp(:,:)
           enddo
        case default
           write(*,*) 'global_collect:invalid ldim ',ldim
           stop '** error in global_collect '
        end select

        deallocate( tmp )
       
        return
        end subroutine global_collect_3d_r8




        subroutine global_collect_3d_z16( mat, ldim )

!       -----------------------------------------
!       perform collection one matrix  at the time
!       -----------------------------------------

        complex , dimension(:,:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        complex , dimension(:,:), allocatable :: tmp
        integer :: j

        select case (ldim)
        case (1,2)
           allocate( tmp(size(mat,1),size(mat,2)) )
           do j=lbound(mat,3),ubound(mat,3)
              tmp(:,:) = mat(:,:,j)
              call global_collect( tmp, ldim )
              mat(:,:,j) = tmp(:,:)
           enddo
        case (3)
           allocate( tmp(size(mat,1),size(mat,3)) )
           do j=lbound(mat,2),ubound(mat,2)
               tmp(:,:)  = mat(:,j,:)
               call global_collect( tmp, 2 )
               mat(:,j,:) = tmp(:,:)
           enddo
        case default
           write(*,*) 'global_collect:invalid ldim ',ldim
           stop '** error in global_collect '
        end select

        deallocate( tmp )
       
        return
        end subroutine global_collect_3d_z16






        subroutine global_collect_2d_r8( mat, ldim )

!       ----------------------------
!       collect one matrix at a time
!       ----------------------------
        real , dimension(:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        real , dimension(:,:), allocatable :: tmp
        integer :: lb1,ub1,lb2,ub2
        integer :: jstart,jend,jsize
        integer, dimension(0:(nprocs-1)) :: istart,iend

        lb1 = lbound(mat,1)
        ub1 = ubound(mat,1)
        lb2 = lbound(mat,2)
        ub2 = ubound(mat,2)

        allocate( tmp(lb1:ub1,lb2:ub2) )
        tmp(:,:) = 0.0

        call parallel_partition( size(mat,ldim), istart,iend )

        select case (ldim)
        case (1)
            jstart = (lb1-1) + istart(iam)
            jend = (lb1-1) + iend(iam)
            jsize = jend-jstart+1
            if (jsize.ge.1) then
              tmp( jstart:jend,: ) = mat( jstart:jend,:)
            endif
        case (2)
            jstart = (lb2-1) + istart(iam)
            jend =   (lb2-1) + iend(iam)
            jsize = jend-jstart+1
            if (jsize.ge.1) then
               tmp( :, jstart:jend ) = mat(:, jstart:jend )
            endif
        case default
            write(*,*) 'global_collect:invalid ldim ',ldim
            stop '** error in global_collect ** '
        end select

        call global_sum( tmp )
        mat(:,:) = tmp(:,:)
        deallocate( tmp )
       
        return
        end subroutine global_collect_2d_r8




        subroutine global_collect_2d_z16( mat, ldim )

!       ----------------------------
!       collect one matrix at a time
!       ----------------------------
        complex , dimension(:,:), intent(inout) :: mat
        integer, intent(in) :: ldim

        complex , dimension(:,:), allocatable :: tmp
        integer :: lb1,ub1,lb2,ub2
        integer :: jstart,jend,jsize
        integer, dimension(0:(nprocs-1)) :: istart,iend

        lb1 = lbound(mat,1)
        ub1 = ubound(mat,1)
        lb2 = lbound(mat,2)
        ub2 = ubound(mat,2)

        allocate( tmp(lb1:ub1,lb2:ub2) )
        tmp(:,:) = 0.0

        call parallel_partition( size(mat,ldim), istart,iend )

        select case (ldim)
        case (1)
            jstart = (lb1-1) + istart(iam)
            jend = (lb1-1) + iend(iam)
            jsize = jend-jstart+1
            if (jsize.ge.1) then
              tmp( jstart:jend,: ) = mat( jstart:jend,:)
            endif
        case (2)
            jstart = (lb2-1) + istart(iam)
            jend =   (lb2-1) + iend(iam)
            jsize = jend-jstart+1
            if (jsize.ge.1) then
               tmp( :, jstart:jend ) = mat(:, jstart:jend )
            endif
        case default
            write(*,*) 'global_collect:invalid ldim ',ldim
            stop '** error in global_collect ** '
        end select

        call global_sum( tmp )
        mat(:,:) = tmp(:,:)
        deallocate( tmp )
       
        return
        end subroutine global_collect_2d_z16
 
    end module blacs_mod
