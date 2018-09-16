      subroutine gen4a(r,p4,wt,*)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'debug.f'
      real(dp):: p4(mxpart,4),p3(mxpart,4),r(mxdim),wt,wt3
c----although wt4 is generated we will not use it
      call gen3(r,p3,wt3,*999)      
      if (debug) write(6,*) 'wt3 in gen4a',wt3
c----this generates the full weight from both branchings
      call gen4from3(p3,r(8),r(9),r(10),p4,wt,*999)     
c     write(6,*) 's45',2._dp
c    .*(p5(4,4)*p5(5,4)-p5(4,1)*p5(5,1)-p5(4,2)*p5(5,2)-p5(4,3)*p5(5,3))
c      pause
      return
 999  wt=0._dp
      return 1
      end
