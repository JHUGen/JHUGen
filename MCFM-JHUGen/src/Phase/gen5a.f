      subroutine gen5a(r,p5,wt,*)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'debug.f'
      real(dp):: p5(mxpart,4),p4(mxpart,4),r(mxdim),wt,wt4
c----although wt4 is generated we will not use it
      call gen4(r,p4,wt4,*999)      
c     write(6,*) 's45',2._dp
c    .*(p4(4,4)*p4(5,4)-p4(4,1)*p4(5,1)-p4(4,2)*p4(5,2)-p4(4,3)*p4(5,3))
      if (debug) write(6,*) 'wt4 in gen5a',wt4

c----this generates the full weight from both branchings
      call gen5from4(p4,r(11),r(12),r(13),p5,wt,*999)      
c     write(6,*) 's45',2._dp
c    .*(p5(4,4)*p5(5,4)-p5(4,1)*p5(5,1)-p5(4,2)*p5(5,2)-p5(4,3)*p5(5,3))
c      pause
      return
 999  wt=0._dp
      return 1
      end
