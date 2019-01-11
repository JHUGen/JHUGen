      subroutine pvswap(q1,q2,m1s,m2s)
      implicit none
      double precision q1(4),q2(4),qtemp(4),m1s,m2s,mtemps
      integer nu
      do nu=1,4
      qtemp(nu)=q1(nu)
      q1(nu)=q2(nu)
      q2(nu)=qtemp(nu)
      enddo
      mtemps=m1s
      m1s=m2s
      m2s=mtemps
      return
      end
