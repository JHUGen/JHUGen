      subroutine qqb_WZbj(p,msq)
      implicit none
      include 'types.f'

!---  Author: R.K. Ellis February 2013
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'nwz.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),WZbbmsq

      call spinoru(8,p,za,zb)
      msq(:,:)=0._dp

      if (nwz == 1) then
      msq(+2,+5)=WZbbmsq(1,8,3,4,5,6,7,2)
      msq(+4,+5)=msq(+2,+5)
      msq(+5,+2)=WZbbmsq(2,8,3,4,5,6,7,1)
      msq(+5,+4)=msq(+5,+2)

      msq(-1,+5)=WZbbmsq(8,1,3,4,5,6,7,2)
      msq(-3,+5)=msq(-1,+5)
      msq(+5,-1)=WZbbmsq(8,2,3,4,5,6,7,1)
      msq(+5,-3)=msq(+5,-1)

      msq(+2,-5)=WZbbmsq(1,7,3,4,5,6,8,2)
      msq(+4,-5)=msq(+2,+5)
      msq(-5,+2)=WZbbmsq(2,7,3,4,5,6,8,1)
      msq(-5,+4)=msq(+5,+2)

      msq(-1,-5)=WZbbmsq(7,1,3,4,5,6,8,2)
      msq(-3,-5)=msq(-1,+5)
      msq(-5,-1)=WZbbmsq(7,2,3,4,5,6,8,1)
      msq(-5,-3)=msq(+5,-1)

      elseif (nwz == -1) then
      msq(+1,+5)=WZbbmsq(1,8,3,4,5,6,7,2)
      msq(+3,+5)=msq(+1,+5)
      msq(+5,+1)=WZbbmsq(2,8,3,4,5,6,7,1)
      msq(+5,+3)=msq(+5,+1)

      msq(-2,+5)=WZbbmsq(8,1,3,4,5,6,7,2)
      msq(-4,+5)=msq(-2,+5)
      msq(+5,-2)=WZbbmsq(8,2,3,4,5,6,7,1)
      msq(+5,-4)=msq(+5,-2)

      msq(+1,-5)=WZbbmsq(1,7,3,4,5,6,8,2)
      msq(+3,-5)=msq(+1,+5)
      msq(-5,+1)=WZbbmsq(2,7,3,4,5,6,8,1)
      msq(-5,+3)=msq(+5,+1)

      msq(-2,-5)=WZbbmsq(7,1,3,4,5,6,8,2)
      msq(-4,-5)=msq(-2,+5)
      msq(-5,-2)=WZbbmsq(7,2,3,4,5,6,8,1)
      msq(-5,-4)=msq(+5,-2)

      endif


      return
      end

