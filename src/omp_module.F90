!> OpenMP module
!! Provides a helper function for OpenMP use
!! that divides a loop into chunks for OMP parallelization
module omp_module
  contains
    !> Split loop range into uniformly sized patches
    subroutine split_indices(total,num_pieces,ibeg,iend)
      implicit none

      integer :: total
      integer :: num_pieces
      integer :: ibeg(num_pieces), iend(num_pieces)
      integer :: itmp1, itmp2, ioffset, i

      if (num_pieces > 0) then
         itmp1 = total/num_pieces
         itmp2 = mod(total,num_pieces)
         ioffset = 0
         do i=1,itmp2
            ibeg(i) = ioffset + 1
            iend(i) = ioffset + (itmp1+1)
            ioffset = iend(i)
         enddo
         do i=itmp2+1,num_pieces
            ibeg(i) = ioffset + 1
            if (ibeg(i) > total) then
               iend(i) = ibeg(i) - 1
            else
               iend(i) = ioffset + itmp1
               ioffset = iend(i)
            endif
         enddo
      endif

    end subroutine split_indices
end module omp_module

