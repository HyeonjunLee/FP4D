module FP4D_math
  use FP4D_globals

  contains

    subroutine cal_dfdr(f_in,dfdr)
      implicit none
      real(kind=8), dimension(ir_s  :ir_e  ) :: f_in, dfdr
      real(kind=8), dimension(ir_s-3:ir_e+3) :: f
      real(kind=8) :: tmpr1, tmpr2
      integer :: index_I
      ! For finite difference parameter
      real (kind=8) :: delr

      dfdr = 0D0
      delr = 60D0*norm_dr
      f(ir_s:ir_e)=f_in
      f(ir_s-1)=f_in(ir_s); f(ir_e+1)=f_in(ir_e-1)
      f(ir_s-2)=f_in(ir_s); f(ir_e+2)=f_in(ir_e-2)
      f(ir_s-3)=f_in(ir_s); f(ir_e+3)=f_in(ir_e-3)

      do index_I=ir_s,ir_e
            tmpr1 =      f(index_I-3) &
                  + 45D0*f(index_I-1) &
                  +  9D0*f(index_I+2)
            tmpr2 =  9D0*f(index_I-2) &
                  + 45D0*f(index_I+1) &
                  +      f(index_I+3)
            dfdr(index_I) = (tmpr2-tmpr1)/delr
      enddo
      return

    end subroutine cal_dfdr

    subroutine cal_dfdz(f_in,dfdz)
      implicit none
      real(kind=8), dimension(iz_s  :iz_e  ) :: f_in, dfdz
      real(kind=8), dimension(iz_s-3:iz_s+3) :: f
      real(kind=8) :: tmpz1,tmpz2
      integer :: index_J
      ! For finite difference parameter
      real (kind=8) :: delz

      dfdz=0D0
      delz  = 60D0*norm_dz
      f(iz_s:iz_e)=f_in
      f(iz_s-1)=f_in(iz_s+1); f(iz_s+1)=f_in(iz_s-1)
      f(iz_s-2)=f_in(iz_s+2); f(iz_s+2)=f_in(iz_s-2)
      f(iz_s-3)=f_in(iz_s+3); f(iz_s+3)=f_in(iz_s-3)

        do index_J=iz_s,iz_e
            tmpz1 =      f(index_J-3) &
                  + 45D0*f(index_J-1) &
                  +  9D0*f(index_J+2)
            tmpz2 =  9D0*f(index_J-2) &
                  + 45D0*f(index_J+1) &
                  +      f(index_J+3)
            dfdz(index_J) = (tmpz2-tmpz1)/delz
        enddo

      return

    end subroutine cal_dfdz

    ! subroutine cal_dfdth2(f_in,dfdth)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdth
    !   real(kind=8) :: tmpt1,tmpt2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delth
    !   dfdth=0D0
    !   delth =eq_dtheta
    !   do index_H=ith_s,ith_e
    !     do index_I=ir_s,ir_e
    !       do index_J=iz_s,iz_e
    !         if (index_H.eq.0) then
    !           tmpt1 = 1.0/2.0*f_in(ith_e,index_I,index_J)
    !           tmpt2 = 1.0/2.0*f_in(index_H+1,index_I,index_J)
    !         else if (index_H.eq.ith_e)then
    !           tmpt1 = 1.0/2.0*f_in(index_H-1,index_I,index_J)
    !           tmpt2 = 1.0/2.0*f_in(        0,index_I,index_J)
    !         else
    !           tmpt1 = 1.0/2.0*f_in(index_H-1,index_I,index_J)
    !           tmpt2 = 1.0/2.0*f_in(index_H+1,index_I,index_J)
    !         endif
    !         dfdth(index_H,index_I,index_J) = (tmpt2-tmpt1)/delth
    !       enddo
    !     enddo
    !   enddo 
      
    !   return

    ! end subroutine cal_dfdth2

    ! subroutine cal_dfdth4(f_in,dfdth)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdth
    !   real(kind=8) :: tmpt1,tmpt2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delth
    !   dfdth=0D0
    !   delth =eq_dtheta
    !   do index_H=ith_s,ith_e
    !     do index_I=ir_s,ir_e
    !       do index_J=iz_s,iz_e
    !         if (index_H.eq.0) then
    !           tmpt1 =  2.0/3.0*f_in(ith_e  ,index_I,index_J) &
    !                 + 1.0/12.0*f_in(  index_H+2,index_I,index_J)
    !           tmpt2 = 1.0/12.0*f_in(ith_e-1,index_I,index_J) &
    !                 +  2.0/3.0*f_in(  index_H+1,index_I,index_J)
    !         else if (index_H.eq.1) then
    !           tmpt1 =  2.0/3.0*f_in(index_H-1,index_I,index_J) &
    !                 + 1.0/12.0*f_in(index_H+2,index_I,index_J)
    !           tmpt2 = 1.0/12.0*f_in(ith_e,index_I,index_J) &
    !                 +  2.0/3.0*f_in(index_H+1,index_I,index_J)
    !         else if (index_H.eq.ith_e)then
    !           tmpt1 =  2.0/3.0*f_in(index_H-1,index_I,index_J) &
    !                 + 1.0/12.0*f_in(        1,index_I,index_J)
    !           tmpt2 = 1.0/12.0*f_in(index_H-2,index_I,index_J) &
    !                 +  2.0/3.0*f_in(        0,index_I,index_J)
    !         else if (index_H.eq.(ith_e-1))then
    !           tmpt1 =  2.0/3.0*f_in(index_H-1,index_I,index_J) &
    !                 + 1.0/12.0*f_in(        0,index_I,index_J)
    !           tmpt2 = 1.0/12.0*f_in(index_H-2,index_I,index_J) &
    !                 +  2.0/3.0*f_in(index_H+1,index_I,index_J)
    !         else
    !           tmpt1 =  2.0/3.0*f_in(index_H-1,index_I,index_J) &
    !                 + 1.0/12.0*f_in(index_H+2,index_I,index_J)
    !           tmpt2 = 1.0/12.0*f_in(index_H-2,index_I,index_J) &
    !                 +  2.0/3.0*f_in(index_H+1,index_I,index_J)
    !         endif
    !         dfdth(index_H,index_I,index_J) = (tmpt2-tmpt1)/delth
    !       enddo
    !     enddo
    !   enddo 
      
    !   return

    ! end subroutine cal_dfdth4


   subroutine cal_dfdth6(f_in,dfdth)
      implicit none
      real(kind=8), dimension(ith_s:ith_e) :: f_in, dfdth
      real(kind=8) :: tmpt1,tmpt2
      integer :: index_H
      ! For finite difference parameter
      real (kind=8) :: delth




      dfdth=0D0
      delth =eq_dtheta*60.0

      do index_H=ith_s,ith_e

           if (index_H.eq.0) then
              tmpt1 =      f_in(ith_e-2) &
                    + 45D0*f_in(ith_e  ) &
                    +  9D0*f_in(ith_s+2)
              tmpt2 =  9D0*f_in(ith_e-1) &
                    + 45D0*f_in(ith_s+1) &
                    +      f_in(ith_s+3)
            else if (index_H.eq.1) then
              tmpt1 =      f_in(ith_e-1) &
                    + 45D0*f_in(ith_s+0) &
                    +  9D0*f_in(ith_s+3)
              tmpt2 =  9D0*f_in(ith_e  ) &
                    + 45D0*f_in(ith_s+2) &
                    +      f_in(ith_s+4)
            else if (index_H.eq.2) then
              tmpt1 =      f_in(ith_e  ) &
                    + 45D0*f_in(ith_s+1) &
                    +  9D0*f_in(ith_s+4)
              tmpt2 =  9D0*f_in(ith_s+0) &
                    + 45D0*f_in(ith_s+3) &
                    +      f_in(ith_s+5)
            else if (index_H.eq.ith_e)then
              tmpt1 =      f_in(ith_e-3) &
                    + 45D0*f_in(ith_e-1) &
                    +  9D0*f_in(ith_s+1)
              tmpt2 =  9D0*f_in(ith_e-2) &
                    + 45D0*f_in(ith_s+0) &
                    +      f_in(ith_s+2)
            else if (index_H.eq.(ith_e-1))then
              tmpt1 =      f_in(ith_e-4) &
                    + 45D0*f_in(ith_e-2) &
                    +  9D0*f_in(ith_s+0)
              tmpt2 =  9D0*f_in(ith_e-3) &
                    + 45D0*f_in(ith_e  ) &
                    +      f_in(ith_s+1)
            else if (index_H.eq.(ith_e-2))then
              tmpt1 =      f_in(ith_e-5) &
                    + 45D0*f_in(ith_e-3) &
                    +  9D0*f_in(ith_e  )
              tmpt2 =  9D0*f_in(ith_e-4) &
                    + 45D0*f_in(ith_e-1) &
                    +      f_in(ith_s+0)
            else
              tmpt1 =      f_in(index_H-3) &
                    + 45D0*f_in(index_H-1) &
                    +  9D0*f_in(index_H+2)
              tmpt2 =  9D0*f_in(index_H-2) &
                    + 45D0*f_in(index_H+1) &
                    +      f_in(index_H+3)
            endif
            dfdth(index_H) = (tmpt2-tmpt1)/delth

      enddo 
      
      return

    end subroutine cal_dfdth6

    subroutine cal_dfdz8(f_in,dfdz)
      implicit none
      real(kind=8), dimension (iz_s  :iz_e ) :: f_in, dfdz
      real(kind=8), dimension (iz_s-4:iz_e+4) :: f
      real(kind=8) :: tmpz1,tmpz2
      integer :: index_J
      ! For finite difference parameter
      real (kind=8) :: delz
      dfdz = 0D0; f=0d0
      delz = norm_dz
      f(iz_s:iz_e)=f_in
!      f(iz_s-1)=f_in(iz_s); f(iz_s+1)=f_in(iz_s)
!      f(iz_s-2)=f_in(iz_s); f(iz_s+2)=f_in(iz_s)
!      f(iz_s-3)=f_in(iz_s); f(iz_s+3)=f_in(iz_s)
!      f(iz_s-4)=f_in(iz_s); f(iz_s+4)=f_in(iz_s)
          do index_J=iz_s,iz_e
            tmpz1 = 4.0/105.0*f(index_J-3) &
                  +   4.0/5.0*f(index_J-1) &
                  +   1.0/5.0*f(index_J+2) &
                  + 1.0/280.0*f(index_J+4)
            tmpz2 = 1.0/280.0*f(index_J-4) &
                  +   1.0/5.0*f(index_J-2) &
                  +   4.0/5.0*f(index_J+1) &
                  + 4.0/105.0*f(index_J+3)
            dfdz(index_J) = (tmpz2-tmpz1)/delz
          enddo
      return

    end subroutine cal_dfdz8

    ! subroutine cal_dfdz8f(f_in,dfdz)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdz
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e+6) :: f
    !   real(kind=8) :: tmpz1,tmpz2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delz
    !   dfdz=0d0; f=0d0;
    !   delz =norm_dz;
    !   f(:,ir_s:ir_e,iz_s:iz_e)=f_in
    !   do index_H=ith_s,ith_e
    !     do index_J=iz_s,iz_e
    !       do index_I=ir_s,ir_e
    !         tmpz1 = 49.0/20.0*f(index_H,index_I,index_J  ) &
    !               +  15.0/2.0*f(index_H,index_I,index_J+2) &
    !               +  15.0/4.0*f(index_H,index_I,index_J+4) &
    !               +   1.0/6.0*f(index_H,index_I,index_J+6)
    !         tmpz2 =       6.0*f(index_H,index_I,index_J+1) &
    !               +  20.0/3.0*f(index_H,index_I,index_J+3) &
    !               +   6.0/5.0*f(index_H,index_I,index_J+5)
    !         dfdz(index_H,index_I,index_J) = (tmpz2-tmpz1)/delz
    !       enddo
    !     enddo
    !   enddo 
    !   return

    ! end subroutine cal_dfdz8f

    ! subroutine cal_dfdz8b(f_in,dfdz)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdz
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s-6:iz_s) :: f
    !   real(kind=8) :: tmpz1,tmpz2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delz
    !   dfdz=0d0; f=0d0;
    !   delz =norm_dz;
    !   f(:,ir_s:ir_e,iz_s:iz_e)=f_in
    !   do index_H=ith_s,ith_e
    !     do index_J=iz_s,iz_e
    !       do index_I=ir_s,ir_e
    !         tmpz1 = 49.0/20.0*f(index_H,index_I,index_J  ) &
    !               +  15.0/2.0*f(index_H,index_I,index_J-2) &
    !               +  15.0/4.0*f(index_H,index_I,index_J-4) &
    !               +   1.0/6.0*f(index_H,index_I,index_J-6)
    !         tmpz2 =       6.0*f(index_H,index_I,index_J-1) &
    !               +  20.0/3.0*f(index_H,index_I,index_J-3) &
    !               +   6.0/5.0*f(index_H,index_I,index_J-5)
    !         dfdz(index_H,index_I,index_J) = (tmpz1-tmpz2)/delz
    !       enddo
    !     enddo
    !   enddo 
    !   return

    ! end subroutine cal_dfdz8b


    ! subroutine cal_dfdr8f(f_in,dfdr)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdr
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e+6,iz_s:iz_e) :: f
    !   real(kind=8) :: tmpr1,tmpr2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delr
    !   dfdr=0d0; f=0d0;
    !   delr =norm_dr;
    !   f(:,ir_s:ir_e,iz_s:iz_e)=f_in
    !   do index_H=ith_s,ith_e
    !     do index_J=iz_s,iz_e
    !       do index_I=ir_s,ir_e
    !         tmpr1 = 49.0/20.0*f(index_H,index_I  ,index_J) &
    !               +  15.0/2.0*f(index_H,index_I+2,index_J) &
    !               +  15.0/4.0*f(index_H,index_I+4,index_J) &
    !               +   1.0/6.0*f(index_H,index_I+6,index_J)
    !         tmpr2 =       6.0*f(index_H,index_I+1,index_J) &
    !               +  20.0/3.0*f(index_H,index_I+3,index_J) &
    !               +   6.0/5.0*f(index_H,index_I+5,index_J)
    !         dfdr(index_H,index_I,index_J) = (tmpr2-tmpr1)/delr
    !       enddo
    !     enddo
    !   enddo 
    !   return

    ! end subroutine cal_dfdr8f

    ! subroutine cal_dfdr8b(f_in,dfdr)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdr
    !   real(kind=8), dimension(ith_s:ith_e,-6:ir_e,iz_s:iz_e) :: f
    !   real(kind=8) :: tmpr1,tmpr2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delr
    !   dfdr=0d0; f=0d0;
    !   delr =norm_dr;
    !   f(:,ir_s:ir_e,iz_s:iz_e)=f_in
    !   f(:,-1,:)=f_in(:,0,:); f(:,-4,:)=f_in(:,0,:)
    !   f(:,-2,:)=f_in(:,0,:); f(:,-5,:)=f_in(:,0,:)
    !   f(:,-3,:)=f_in(:,0,:); f(:,-6,:)=f_in(:,0,:)
    !   do index_H=ith_s,ith_e
    !     do index_J=iz_s,iz_e
    !       do index_I=ir_s,ir_e
    !         tmpr1 = 49.0/20.0*f(index_H,index_I  ,index_J) &
    !               +  15.0/2.0*f(index_H,index_I-2,index_J) &
    !               +  15.0/4.0*f(index_H,index_I-4,index_J) &
    !               +   1.0/6.0*f(index_H,index_I-6,index_J)
    !         tmpr2 =       6.0*f(index_H,index_I-1,index_J) &
    !               +  20.0/3.0*f(index_H,index_I-3,index_J) &
    !               +   6.0/5.0*f(index_H,index_I-5,index_J)
    !         dfdr(index_H,index_I,index_J) = (tmpr1-tmpr2)/delr
    !       enddo
    !     enddo
    !   enddo 
    !   return

    ! end subroutine cal_dfdr8b



    subroutine cal_dfdr8(f_in,dfdr)
      implicit none
      real(kind=8), dimension(ir_s:ir_e) :: f_in, dfdr
      real(kind=8), dimension(-4:ir_e+4) :: f
      real(kind=8) :: tmpr1,tmpr2
      integer :: index_I
      ! For finite difference parameter
      real (kind=8) :: delr
      dfdr=0d0; f=0d0;
      delr =norm_dr;
      f(ir_s:ir_e)=f_in
      f(-1)=f_in(0); f(-3)=f_in(0); 
      f(-2)=f_in(0); f(-4)=f_in(0);
!     f(:,ir_e+1,:)=f_in(:,ir_e,:); f(:,ir_e+3,:)=f_in(:,ir_e,:)
!     f(:,ir_e+2,:)=f_in(:,ir_e,:); f(:,ir_e+4,:)=f_in(:,ir_e,:)

          do index_I=ir_s,ir_e
            tmpr1 = 4.0/105.0*f(index_I-3) &
                  +   4.0/5.0*f(index_I-1) &
                  +   1.0/5.0*f(index_I+2) &
                  + 1.0/280.0*f(index_I+4)
            tmpr2 = 1.0/280.0*f(index_I-4) &
                  +   1.0/5.0*f(index_I-2) &
                  +   4.0/5.0*f(index_I+1) &
                  + 4.0/105.0*f(index_I+3)
            dfdr(index_I) = (tmpr2-tmpr1)/delr
          enddo

      return

    end subroutine cal_dfdr8

!    subroutine cal_dfdth8(f_in,dfdth)
!       implicit none
!       real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, dfdth
!       real(kind=8) :: tmpt1,tmpt2
!       integer :: index_H,index_I,index_J
!       ! For finite difference parameter
!       real (kind=8) :: delth
!       dfdth=0D0
!       delth =eq_dtheta
!       do index_H=ith_s,ith_e
!         do index_I=ir_s,ir_e
!           do index_J=iz_s,iz_e
!             if (index_H.eq.0) then
!               tmpt1 = 4.0/105.0*f_in(ith_e-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(ith_e,index_I,index_J) &
!                     +   1.0/5.0*f_in(  index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(  index_H+4,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(ith_e-3,index_I,index_J) &
!                     +   1.0/5.0*f_in(ith_e-1,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(  index_H+3,index_I,index_J)
!             else if (index_H.eq.1) then
!               tmpt1 = 4.0/105.0*f_in(ith_e-1,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(  index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(  index_H+4,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(ith_e-2,index_I,index_J) &
!                     +   1.0/5.0*f_in(  ith_e,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(  index_H+3,index_I,index_J)
!             else if (index_H.eq.2) then
!               tmpt1 = 4.0/105.0*f_in(  ith_e,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(  index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(  index_H+4,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(ith_e-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(  index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(  index_H+3,index_I,index_J)
!             else if (index_H.eq.3) then
!               tmpt1 = 4.0/105.0*f_in(  index_H-3,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(  index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(  index_H+4,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(  ith_e,index_I,index_J) &
!                     +   1.0/5.0*f_in(  index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(  index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(  index_H+3,index_I,index_J)
!             else if (index_H.eq.ith_e) then
!               tmpt1 = 4.0/105.0*f_in(index_H-3,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(        1,index_I,index_J) &
!                     + 1.0/280.0*f_in(        3,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(index_H-4,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(        0,index_I,index_J) &
!                     + 4.0/105.0*f_in(        2,index_I,index_J)
!             else if (index_H.eq.(ith_e-1)) then
!               tmpt1 = 4.0/105.0*f_in(index_H-3,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(        0,index_I,index_J) &
!                     + 1.0/280.0*f_in(        2,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(index_H-4,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(        1,index_I,index_J)
!             else if (index_H.eq.(ith_e-2)) then
!               tmpt1 = 4.0/105.0*f_in(index_H-3,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(        1,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(index_H-4,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(        0,index_I,index_J)
!             else if (index_H.eq.(ith_e-3)) then
!               tmpt1 = 4.0/105.0*f_in(index_H-3,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(        0,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(index_H-4,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(index_H+3,index_I,index_J)
!             else
!               tmpt1 = 4.0/105.0*f_in(index_H-3,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H-1,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H+2,index_I,index_J) &
!                     + 1.0/280.0*f_in(index_H+4,index_I,index_J)
!               tmpt2 = 1.0/280.0*f_in(index_H-4,index_I,index_J) &
!                     +   1.0/5.0*f_in(index_H-2,index_I,index_J) &
!                     +   4.0/5.0*f_in(index_H+1,index_I,index_J) &
!                     + 4.0/105.0*f_in(index_H+3,index_I,index_J)
!             endif
!             dfdth(index_H,index_I,index_J) = (tmpt2-tmpt1)/delth
!           enddo
!         enddo
!       enddo 
      
!       return

!     end subroutine cal_dfdth8

    subroutine cal_d2fdr28(f_in,d2fdr2)
      implicit none
      real(kind=8), dimension(ir_s  :ir_e  ) :: f_in, d2fdr2
      real(kind=8), dimension(ir_s-4:ir_e+4) :: f
      real(kind=8) :: tmpr3, tmpr4
      integer :: index_I
      ! For finite difference parameter
      real (kind=8) :: delr
      d2fdr2=0D0;f=0d0;
      delr  = norm_dr
      f(ir_s:ir_e)=f_in
      f(-1)=f_in(0); f(-3)=f_in(0)
      f(-2)=f_in(0); f(-4)=f_in(0)
!     f(:,ir_e+1,:)=f_in(:,ir_e,:); f(:,ir_e+3,:)=f_in(:,ir_e,:)
!     f(:,ir_e+2,:)=f_in(:,ir_e,:); f(:,ir_e+4,:)=f_in(:,ir_e,:)
        do index_I=ir_s,ir_e
            tmpr3 = 1.0/560.0*f(index_I-4) &
                  +   1.0/5.0*f(index_I-2) &
                  +205.0/72.0*f(index_I  ) &
                  +   1.0/5.0*f(index_I+2) &
                  + 1.0/560.0*f(index_I+4)
            tmpr4 = 8.0/315.0*f(index_I-3) &
                  +   8.0/5.0*f(index_I-1) &
                  +   8.0/5.0*f(index_I+1) &
                  + 8.0/315.0*f(index_I+3)
            d2fdr2(index_I) =(tmpr4-tmpr3)/delr/delr
        enddo
      
      return

    end subroutine cal_d2fdr28

    subroutine cal_d2fdz28(f_in,d2fdz2)
      implicit none
      real(kind=8), dimension(iz_s  :iz_e  ) :: f_in,d2fdz2
      real(kind=8), dimension(iz_s-4:iz_e+4) :: f
      real(kind=8) :: tmpz3,tmpz4
      integer :: index_J
      ! For finite difference parameter
      real (kind=8) :: delz
      d2fdz2=0D0; f=0d0
      delz  = norm_dz
      f(iz_s:iz_e)=f_in
!      f(iz_s-1)=f_in(iz_s); f(iz_s+1)=f_in(iz_s)
!      f(iz_s-2)=f_in(iz_s); f(iz_s+2)=f_in(iz_s)
!      f(iz_s-3)=f_in(iz_s); f(iz_s+3)=f_in(iz_s)
!      f(iz_s-4)=f_in(iz_s); f(iz_s+4)=f_in(iz_s)
  
          do index_J=iz_s,iz_e
            tmpz3 = 1.0/560.0*f(index_J-4) &
                  +   1.0/5.0*f(index_J-2) &
                  +205.0/72.0*f(index_J) &
                  +   1.0/5.0*f(index_J+2) &
                  + 1.0/560.0*f(index_J+4)
            tmpz4 = 8.0/315.0*f(index_J-3) &
                  +   8.0/5.0*f(index_J-1) &
                  +   8.0/5.0*f(index_J+1) &
                  + 8.0/315.0*f(index_J+3)
            d2fdz2(index_J) = (tmpz4-tmpz3)/delz/delz
          enddo

      
      return

    end subroutine cal_d2fdz28


    ! subroutine cal_d2fdr2(f_in,d2fdr2)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, d2fdr2
    !   real(kind=8), dimension(ith_s:ith_e,-3:ir_e+3,iz_s:iz_e) :: f
    !   real(kind=8) :: tmpr3, tmpr4
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delr
    !   d2fdr2=0D0
    !   delr  = 60D0*norm_dr
    !   f(:,ir_s:ir_e,iz_s:iz_e)=f_in
    !   f(:,-1,:)=f_in(:,1,:); f(:,ir_e+1,:)=f_in(:,ir_e-1,:)
    !   f(:,-2,:)=f_in(:,2,:); f(:,ir_e+2,:)=f_in(:,ir_e-2,:)
    !   f(:,-3,:)=f_in(:,3,:); f(:,ir_e+3,:)=f_in(:,ir_e-3,:)  
    !   do index_H=ith_s,ith_e
    !     do index_I=ir_s,ir_e
    !       do index_J=iz_s,iz_e
    !         tmpr3 = 54D0*f(index_H,index_I-2,index_J) &
    !               +980D0*f(index_H,index_I  ,index_J) &
    !               + 54D0*f(index_H,index_I+2,index_J)
    !         tmpr4 =  4D0*f(index_H,index_I-3,index_J) &
    !               +540D0*f(index_H,index_I-1,index_J) &
    !               +540D0*f(index_H,index_I+1,index_J) &
    !               +  4D0*f(index_H,index_I+3,index_J)
    !         d2fdr2(index_H,index_I,index_J) =(tmpr4-tmpr3)/delr/delr*10.0
    !       enddo
    !     enddo
    !   enddo 
      
    !   return

    ! end subroutine cal_d2fdr2

    ! subroutine cal_d2fdz2(f_in,d2fdz2)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in,d2fdz2
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s-3:iz_s+3) :: f
    !   real(kind=8) :: tmpz3,tmpz4
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delz
    !   d2fdz2=0D0
    !   delz  = 60D0*norm_dz
    !   f(:,ir_s:ir_e,iz_s:iz_e)=f_in
    !   f(:,:,iz_s-1)=f_in(:,:,iz_s+1); f(:,:,iz_s+1)=f_in(:,:,iz_s-1)
    !   f(:,:,iz_s-2)=f_in(:,:,iz_s+2); f(:,:,iz_s+2)=f_in(:,:,iz_s-2)
    !   f(:,:,iz_s-3)=f_in(:,:,iz_s+3); f(:,:,iz_s+3)=f_in(:,:,iz_s-3)
    !  do index_H=ith_s,ith_e
    !     do index_I=ir_s,ir_e
    !       do index_J=iz_s,iz_e
    !         tmpz3 = 54D0*f(index_H,index_I,index_J-2) &
    !               +980D0*f(index_H,index_I,index_J  ) &
    !               + 54D0*f(index_H,index_I,index_J+2)
    !         tmpz4 =  4D0*f(index_H,index_I,index_J-3) &
    !               +540D0*f(index_H,index_I,index_J-1) &
    !               +540D0*f(index_H,index_I,index_J+1) &
    !               +  4D0*f(index_H,index_I,index_J+3)
    !         d2fdz2(index_H,index_I,index_J) = (tmpz4-tmpz3)/delz/delz*10.0
    !       enddo
    !     enddo
    !   enddo 
      
    !   return

    ! end subroutine cal_d2fdz2

    subroutine cal_dfdths(f_in,dfdth2)
      implicit none
      real(kind=8), dimension(ith_s:ith_e) :: f_in, dfdth2
      real(kind=8) :: tmpt1,tmpt2
      integer :: index_H
      ! For finite difference parameter
      real (kind=8) :: delth
      dfdth2=0D0
      delth=eq_dtheta*60.0
      do index_H=ith_s,ith_e
            if (index_H.eq.0) then
              tmpt1 =      f_in(ith_e-2) &
                    + 15D0*f_in(ith_e  ) &
                    + 15D0*f_in(index_H+1  ) &
                    +      f_in(index_H+3  )
              tmpt2 =  6D0*f_in(ith_e-1) &
                    + 20D0*f_in(index_H    ) &
                    +  6D0*f_in(index_H+2  )      
            else if (index_H.eq.1) then
              tmpt1 =      f_in(ith_e-1) &
                    + 15D0*f_in(          0) &
                    + 15D0*f_in(index_H+1  ) &
                    +      f_in(index_H+3  )
              tmpt2 =  6D0*f_in(ith_e  ) &
                    + 20D0*f_in(index_H    ) &
                    +  6D0*f_in(index_H+2  )      
            else if (index_H.eq.2) then
              tmpt1 =      f_in(ith_e) &
                    + 15D0*f_in(index_H-1) &
                    + 15D0*f_in(index_H+1) &
                    +      f_in(index_H+3)
              tmpt2 =  6D0*f_in(index_H-2) &
                    + 20D0*f_in(index_H  ) &
                    +  6D0*f_in(index_H+2)
            else if (index_H.eq.ith_e)then
              tmpt1 =      f_in(index_H-3) &
                    + 15D0*f_in(index_H-1) &
                    + 15D0*f_in(        0) &
                    +      f_in(        2)
              tmpt2 =  6D0*f_in(index_H-2) &
                    + 20D0*f_in(index_H  ) &
                    +  6D0*f_in(        1)
            else if (index_H.eq.(ith_e-1))then
              tmpt1 =      f_in(index_H-3) &
                    + 15D0*f_in(index_H-1) &
                    + 15D0*f_in(index_H+1) &
                    +      f_in(        1)
              tmpt2 =  6D0*f_in(index_H-2) &
                    + 20D0*f_in(index_H  ) &
                    +  6D0*f_in(        0)
            else if (index_H.eq.(ith_e-2))then
              tmpt1 =      f_in(index_H-3) &
                    + 15D0*f_in(index_H-1) &
                    + 15D0*f_in(index_H+1) &
                    +      f_in(        0)
              tmpt2 =  6D0*f_in(index_H-2) &
                    + 20D0*f_in(index_H  ) &
                    +  6D0*f_in(index_H+2)
            else
              tmpt1 =      f_in(index_H-3) &
                    + 15D0*f_in(index_H-1) &
                    + 15D0*f_in(index_H+1) &
                    +      f_in(index_H+3)
              tmpt2 =  6D0*f_in(index_H-2) &
                    + 20D0*f_in(index_H  ) &
                    +  6D0*f_in(index_H+2)
            endif
            dfdth2(index_H) = (tmpt2-tmpt1)/delth
      enddo 
      
      return

    end subroutine cal_dfdths

    ! subroutine cal_d2fdth2(f_in,d2fdth2)
    !   implicit none
    !   real(kind=8), dimension(ith_s:ith_e,ir_s:ir_e,iz_s:iz_e) :: f_in, d2fdth2
    !   real(kind=8) :: tmpt1,tmpt2
    !   integer :: index_H,index_I,index_J
    !   ! For finite difference parameter
    !   real (kind=8) :: delth
    !   d2fdth2=0D0
    !   delth = 280D0*eq_dtheta
    !   do index_H=ith_s,ith_e
    !     do index_I=ir_s,ir_e
    !       do index_J=iz_s,iz_e
    !         if (index_H.eq.0) then
    !           tmpt1 =  8.0*f_in(ith_e-2,index_I,index_J) &
    !                 + 56.0*f_in(ith_e  ,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1  ,index_I,index_J) &
    !                 +  8.0*f_in(index_H+3  ,index_I,index_J)
    !           tmpt2 =      f_in(ith_e-3,index_I,index_J) &
    !                 + 28.0*f_in(ith_e-1,index_I,index_J) &
    !                 + 70.0*f_in(index_H    ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2  ,index_I,index_J) &
    !                 +      f_in(index_H+4  ,index_I,index_J)
    !         else if (index_H.eq.1) then
    !           tmpt1 =  8.0*f_in(ith_e-1,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1  ,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1  ,index_I,index_J) &
    !                 +  8.0*f_in(index_H+3  ,index_I,index_J)
    !           tmpt2 =      f_in(ith_e-2,index_I,index_J) &
    !                 + 28.0*f_in(ith_e  ,index_I,index_J) &
    !                 + 70.0*f_in(index_H    ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2  ,index_I,index_J) &
    !                 +      f_in(index_H+4  ,index_I,index_J)
    !         else if (index_H.eq.2) then
    !           tmpt1 =  8.0*f_in(ith_e  ,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1  ,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1  ,index_I,index_J) &
    !                 +  8.0*f_in(index_H+3  ,index_I,index_J)
    !           tmpt2 =      f_in(ith_e-1,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2  ,index_I,index_J) &
    !                 + 70.0*f_in(index_H    ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2  ,index_I,index_J) &
    !                 +      f_in(index_H+4  ,index_I,index_J)
    !         else if (index_H.eq.3) then
    !           tmpt1 =  8.0*f_in(index_H-3  ,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1  ,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1  ,index_I,index_J) &
    !                 +  8.0*f_in(index_H+3  ,index_I,index_J)
    !           tmpt2 =      f_in(ith_e  ,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2  ,index_I,index_J) &
    !                 + 70.0*f_in(index_H    ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2  ,index_I,index_J) &
    !                 +      f_in(index_H+4  ,index_I,index_J)
    !         else if (index_H.eq.ith_e) then
    !           tmpt1 =  8.0*f_in(index_H-3,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1,index_I,index_J) &
    !                 + 56.0*f_in(0        ,index_I,index_J) &
    !                 +  8.0*f_in(2        ,index_I,index_J)
    !           tmpt2 =      f_in(index_H-4,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2,index_I,index_J) &
    !                 + 70.0*f_in(index_H  ,index_I,index_J) &
    !                 + 28.0*f_in(1        ,index_I,index_J) &
    !                 +      f_in(3        ,index_I,index_J) 
    !         else if (index_H.eq.(ith_e-1)) then
    !           tmpt1 =  8.0*f_in(index_H-3,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1,index_I,index_J) &
    !                 +  8.0*f_in(1        ,index_I,index_J)
    !           tmpt2 =      f_in(index_H-4,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2,index_I,index_J) &
    !                 + 70.0*f_in(index_H  ,index_I,index_J) &
    !                 + 28.0*f_in(0        ,index_I,index_J) &
    !                 +      f_in(2        ,index_I,index_J) 
    !         else if (index_H.eq.(ith_e-2)) then
    !           tmpt1 =  8.0*f_in(index_H-3,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1,index_I,index_J) &
    !                 +  8.0*f_in(0        ,index_I,index_J)
    !           tmpt2 =      f_in(index_H-4,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2,index_I,index_J) &
    !                 + 70.0*f_in(index_H  ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2,index_I,index_J) &
    !                 +      f_in(1        ,index_I,index_J) 
    !         else if (index_H.eq.(ith_e-3)) then
    !           tmpt1 =  8.0*f_in(index_H-3,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1,index_I,index_J) &
    !                 +  8.0*f_in(index_H+3,index_I,index_J)
    !           tmpt2 =      f_in(index_H-4,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2,index_I,index_J) &
    !                 + 70.0*f_in(index_H  ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2,index_I,index_J) &
    !                 +      f_in(0        ,index_I,index_J) 
    !         else
    !           tmpt1 =  8.0*f_in(index_H-3,index_I,index_J) &
    !                 + 56.0*f_in(index_H-1,index_I,index_J) &
    !                 + 56.0*f_in(index_H+1,index_I,index_J) &
    !                 +  8.0*f_in(index_H+3,index_I,index_J)
    !           tmpt2 =      f_in(index_H-4,index_I,index_J) &
    !                 + 28.0*f_in(index_H-2,index_I,index_J) &
    !                 + 70.0*f_in(index_H  ,index_I,index_J) &
    !                 + 28.0*f_in(index_H+2,index_I,index_J) &
    !                 +      f_in(index_H+4,index_I,index_J)
    !         endif
    !        d2fdth2(index_H,index_I,index_J) = (tmpt2-tmpt1)/delth
    !       enddo
    !     enddo
    !   enddo 
      
    !   return

    ! end subroutine cal_d2fdth2


    subroutine cal_dphidth(phival,dphidth_out)
      implicit none
      real(kind=8), intent(in) , dimension(ith_s:ith_e) :: phival
      real(kind=8), intent(out), dimension(ith_s:ith_e) :: dphidth_out
      real(kind=8) :: tmpt1,tmpt2,delth
      integer :: index_H
      ! For finite difference parameter
      dphidth_out=0D0
      delth = 5D0*eq_dtheta
      do index_H=ith_s,ith_e
            if (index_H.eq.0) then
              tmpt1 = 4.0/21.0*phival(ith_e-2) &
                    +      4.0*phival(ith_e) &
                    +      1.0*phival(  index_H+2) &
                    + 1.0/56.0*phival(  index_H+4)
              tmpt2 = 1.0/56.0*phival(ith_e-3) &
                    +      1.0*phival(ith_e-1) &
                    +      4.0*phival(  index_H+1) &
                    + 4.0/21.0*phival(  index_H+3)
            else if (index_H.eq.1) then
              tmpt1 = 4.0/21.0*phival(ith_e-1) &
                    +      4.0*phival(  index_H-1) &
                    +      1.0*phival(  index_H+2) &
                    + 1.0/56.0*phival(  index_H+4)
              tmpt2 = 1.0/56.0*phival(ith_e-2) &
                    +      1.0*phival(  ith_e) &
                    +      4.0*phival(  index_H+1) &
                    + 4.0/21.0*phival(  index_H+3)
            else if (index_H.eq.2) then
              tmpt1 = 4.0/21.0*phival(  ith_e) &
                    +      4.0*phival(  index_H-1) &
                    +      1.0*phival(  index_H+2) &
                    + 1.0/56.0*phival(  index_H+4)
              tmpt2 = 1.0/56.0*phival(ith_e-1) &
                    +      1.0*phival(  index_H-2) &
                    +      4.0*phival(  index_H+1) &
                    + 4.0/21.0*phival(  index_H+3)
            else if (index_H.eq.3) then
              tmpt1 = 4.0/21.0*phival(  index_H-3) &
                    +      4.0*phival(  index_H-1) &
                    +      1.0*phival(  index_H+2) &
                    + 1.0/56.0*phival(  index_H+4)
              tmpt2 = 1.0/56.0*phival(  ith_e) &
                    +      1.0*phival(  index_H-2) &
                    +      4.0*phival(  index_H+1) &
                    + 4.0/21.0*phival(  index_H+3)
            else if (index_H.eq.ith_e) then
              tmpt1 = 4.0/21.0*phival(index_H-3) &
                    +      4.0*phival(index_H-1) &
                    +      1.0*phival(        1) &
                    + 1.0/56.0*phival(        3)
              tmpt2 = 1.0/56.0*phival(index_H-4) &
                    +      1.0*phival(index_H-2) &
                    +      4.0*phival(        0) &
                    + 4.0/21.0*phival(        2)
            else if (index_H.eq.(ith_e-1)) then
              tmpt1 = 4.0/21.0*phival(index_H-3) &
                    +      4.0*phival(index_H-1) &
                    +      1.0*phival(        0) &
                    + 1.0/56.0*phival(        2)
              tmpt2 = 1.0/56.0*phival(index_H-4) &
                    +      1.0*phival(index_H-2) &
                    +      4.0*phival(index_H+1) &
                    + 4.0/21.0*phival(        1)
            else if (index_H.eq.(ith_e-2)) then
              tmpt1 = 4.0/21.0*phival(index_H-3) &
                    +      4.0*phival(index_H-1) &
                    +      1.0*phival(index_H+2) &
                    + 1.0/56.0*phival(        1)
              tmpt2 = 1.0/56.0*phival(index_H-4) &
                    +      1.0*phival(index_H-2) &
                    +      4.0*phival(index_H+1) &
                    + 4.0/21.0*phival(        0)
            else if (index_H.eq.(ith_e-3)) then
              tmpt1 = 4.0/21.0*phival(index_H-3) &
                    +      4.0*phival(index_H-1) &
                    +      1.0*phival(index_H+2) &
                    + 1.0/56.0*phival(        0)
              tmpt2 = 1.0/56.0*phival(index_H-4) &
                    +      1.0*phival(index_H-2) &
                    +      4.0*phival(index_H+1) &
                    + 4.0/21.0*phival(index_H+3)
            else
              tmpt1 = 4.0/21.0*phival(index_H-3) &
                    +      4.0*phival(index_H-1) &
                    +      1.0*phival(index_H+2) &
                    + 1.0/56.0*phival(index_H+4)
              tmpt2 = 1.0/56.0*phival(index_H-4) &
                    +      1.0*phival(index_H-2) &
                    +      4.0*phival(index_H+1) &
                    + 4.0/21.0*phival(index_H+3)
            endif
            dphidth_out(index_H) = (tmpt2-tmpt1)/delth
      enddo 
      
      return

    end subroutine cal_dphidth


   function compute_FSA(x,ip) result(FSA_val)
      ! NEED TO KNOW THE 'ip'
      implicit none
      real(kind=8), intent(in) :: x(ith_s:ith_e)
      real(kind=8) :: FSA_val
      integer, intent(in) :: ip
      integer :: ith

      FSA_val = 0D0
      do ith = ith_s, ith_e
        FSA_val = FSA_val + x(ith) * eq_J(ith,ip) * eq_dtheta
      enddo
      FSA_val = FSA_val / eq_intJ(ip)
    end function compute_FSA


   function delta_gaussian(x,x0,width) result(result)
      use FP4D_globals, only:sml_pi
      implicit none
      ! delta(x-x0) = 1/width/sqrt(pi) * exp(-(x-x0)^2/width^2) for x in [-infty,infty]
      real(kind=8), intent(in) :: x, x0, width
      real(kind=8) :: result
      result = (x-x0)/width
      result = 1D0/width/sqrt(sml_pi) &
            * exp(- result*result)
    end function delta_gaussian
    
    function delta_rectangular(x,x0,width) result(result)
      implicit none
      ! delta(x-x0) = 1/dx for x in [x0-dx/2,x0+dx/2]
      ! delta(x-x0) = 0    for the others
      real(kind=8), intent(in) :: x, x0, width
      real(kind=8) :: result

      if (x .eq. x0) then
            result = 1D0/width
      else
            result = 0D0
      endif
    end function delta_rectangular

end module FP4D_math
