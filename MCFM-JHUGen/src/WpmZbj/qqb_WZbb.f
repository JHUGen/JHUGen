      subroutine qqb_WZbb(p,msq)
      implicit none
!---  Author: R.K. Ellis February 2013
      include 'constants.f'
      include 'zprods_com.f'
      include 'nwz.f'
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),WZbbmsq

      call spinoru(8,p,za,zb)

c--- initialize matrix elements to zero
      msq(:,:)=0d0

      if (nwz .eq. 1) then
      msq(+2,-1)=WZbbmsq(1,2,3,4,5,6,7,8)
      msq(-1,+2)=WZbbmsq(2,1,3,4,5,6,7,8)

      msq(+4,-3)=msq(+2,-1)
      msq(-3,+4)=msq(-1,+2)

      elseif (nwz .eq. -1) then

      msq(+1,-2)=WZbbmsq(1,2,3,4,5,6,7,8)
      msq(-2,+1)=WZbbmsq(2,1,3,4,5,6,7,8)

      msq(+3,-4)=msq(+1,-2)
      msq(-4,+3)=msq(-2,+1)

      endif

      return
      end


