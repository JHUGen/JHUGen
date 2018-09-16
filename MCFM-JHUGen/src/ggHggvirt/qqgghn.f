      subroutine qqgghn(p1,p2,p3,p4,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
      implicit none
      include 'types.f'

c---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p5+p6)+g(p3)+g(p4)
c   with momentum 4 contracted with the vector n(mu),
c   separated into colour orderings AB, BA and the sum
c     calculated by the program qqgghn.frm
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4
      real(dp):: p(mxpart,4),n(4),nDn,nDp1,nDp2,nDp3,nDp4
      real(dp):: s123,s124,qqgghn_ab,qqgghn_ba,qqgghn_sym
      nDp1=n(4)*p(p1,4)-n(3)*p(p1,3)-n(2)*p(p1,2)-n(1)*p(p1,1)
      nDp2=n(4)*p(p2,4)-n(3)*p(p2,3)-n(2)*p(p2,2)-n(1)*p(p2,1)
      nDp3=n(4)*p(p3,4)-n(3)*p(p3,3)-n(2)*p(p3,2)-n(1)*p(p3,1)
      nDp4=n(4)*p(p4,4)-n(3)*p(p4,3)-n(2)*p(p4,2)-n(1)*p(p4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      s123=s(p1,p2)+s(p2,p3)+s(p3,p1)
      s124=s(p1,p2)+s(p2,p4)+s(p4,p1)

      call checkndotp(p,n,p4)
c--- appropriate scale is approx 1.e-3_dp*energy(incoming)
c--- so of order(1) for the Tevatron
c      if (abs(nDp4)>1.e-3_dp*sqrt(abs(p(1,4)*n(4)))) then
c        write(*,*) 'Error in qqgghn for :',p1,p2,p3,p4
c        write(*,*) 'cutoff',1.e-3_dp*sqrt(abs(p(1,4)*n(4)))
c        write(6,*) 'nDp4',nDp4
c        call flush(6)
c        stop
c      endif

      qqgghn_ab=  + s123**(-2)*xn*nDn * ( 6._dp*s(p1,p4)**2 + 12._dp*s(p1
     &    ,p4)*s(p2,p4) + 12._dp*s(p1,p4)*s(p3,p4) + 6._dp*s(p2,p4)**2 +
     &    12._dp*s(p2,p4)*s(p3,p4) + 6._dp*s(p3,p4)**2 + 8._dp/(s(p1,p2))*
     &    s(p1,p3)*s(p1,p4)**2 + 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(
     &    p2,p4) + 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4) + 8._dp/(
     &    s(p1,p2))*s(p1,p3)*s(p2,p4)**2 + 16._dp/(s(p1,p2))*s(p1,p3)*s(
     &    p2,p4)*s(p3,p4) + 8._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4)**2 + 4._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4)**2 + 8._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2
     &    ))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4)*s(p3,p4) + 4._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p1,p3)**2*s(p2,p4)**2 + 8._dp/(s(p1,p2))/(s(p1,
     &    p2))*s(p1,p3)**2*s(p2,p4)*s(p3,p4) + 4._dp/(s(p1,p2))/(s(p1,p2
     &    ))*s(p1,p3)**2*s(p3,p4)**2 + 2._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4
     &    )**2 + 4._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4)*s(p2,p4) + 4._dp/(s(
     &    p1,p3))*s(p1,p2)*s(p1,p4)*s(p3,p4) + 2._dp/(s(p1,p3))*s(p1,p2)
     &    *s(p2,p4)**2 )
      qqgghn_ab = qqgghn_ab + s123**(-2)*xn*nDn * ( 4._dp/(s(p1,p3))*s(
     &    p1,p2)*s(p2,p4)*s(p3,p4) + 2._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4)
     &    **2 )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDp1*nDp2 * (
     &     - 16._dp*s(p3,p4) - 48._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 80._dp
     &    /(s(p1,p2))*s(p1,p4)*s(p3,p4) - 32._dp/(s(p1,p2))*s(p3,p4)**2
     &     - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p3,p4) - 64._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4) - 32._dp/(s(p1
     &    ,p2))/(s(p1,p2))*s(p1,p4)**2*s(p3,p4) - 16._dp/(s(p1,p2))/(s(
     &    p1,p3))*s(p1,p4)**2*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p3))*s(
     &    p1,p4)*s(p3,p4)**2 - 16._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)**2*
     &    s(p3,p4) + 16._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)*s(p3,p4)**2
     &     - 16._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p3))*s(p3
     &    ,p4)**2 - 16._dp/(s(p2,p4))*s(p1,p3)*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDp1*nDp3 * (
     &     - 16._dp*s(p1,p2) - 32._dp*s(p1,p3) - 32._dp*s(p3,p4) - 16._dp/(
     &    s(p1,p2))*s(p1,p3)**2 - 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4) -
     &    16._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) + 16._dp/(s(p1,p2))*s(p1,p4
     &    )**2 - 32._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4) + 16._dp/(s(p1,p2))
     &    /(s(p1,p3))*s(p1,p4)**3 - 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4
     &    )**2*s(p3,p4) + 16._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4) - 16._dp/(
     &    s(p1,p3))*s(p1,p2)*s(p3,p4) + 32._dp/(s(p1,p3))*s(p1,p4)**2 -
     &    32._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDp1**2 * (  -
     &    48._dp*s(p3,p4) - 48._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 64._dp/(
     &    s(p1,p2))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p3)**2*s(p3,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(
     &    p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2*s(
     &    p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)**2*s(p3,p4) -
     &    16._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4) - 32._dp/(s(p1,p3))*s(p1,p4
     &    )*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDp2*nDp3 * (
     &     - 16._dp*s(p1,p3) - 16._dp*s(p1,p4) - 16._dp/(s(p1,p2))*s(p1,p3
     &    )**2 - 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4) + 16._dp/(s(p1,p2))*
     &    s(p1,p3)*s(p3,p4) + 16._dp/(s(p1,p2))*s(p1,p4)**2 + 32._dp/(s(
     &    p1,p2))*s(p1,p4)*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,
     &    p4)**3 + 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)**2*s(p3,p4) +
     &    16._dp/(s(p1,p3))*s(p1,p4)**2 + 16._dp/(s(p1,p3))*s(p1,p4)*s(p3
     &    ,p4) - 16._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3) - 16._dp/(s(p2,p4))*
     &    s(p1,p3)**2 + 16._dp/(s(p2,p4))*s(p1,p3)*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDp2**2 * (  -
     &    16._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,
     &    p2))*s(p1,p3)**2*s(p3,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,
     &    p3)*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)
     &    **2*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)**2*s(p3,
     &    p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDp3**2 * ( 32._dp
     &    *s(p1,p4) + 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,
     &    p2))*s(p1,p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)**3 +
     &    16._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4) + 32._dp/(s(p1,p3))*s(p1,p4
     &    )**2 )
      qqgghn_ab = qqgghn_ab + s123**(-1)*s124**(-1)*xn*nDn * (  - 8._dp*
     &    s(p1,p2)*s(p1,p4) - 4._dp*s(p1,p2)*s(p3,p4) - 8._dp*s(p1,p3)*s(
     &    p1,p4) + 16._dp*s(p1,p4)*s(p3,p4) + 8._dp*s(p3,p4)**2 + 16._dp/(
     &    s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4) - 8._dp/(s(p1,p2))*s(p1,
     &    p4)*s(p3,p4)**2 - 4._dp/(s(p1,p2))*s(p3,p4)**3 - 8._dp/(s(p1,p2
     &    ))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4)**2 - 2._dp/(s(p1,p3))
     &    *s(p1,p2)**2*s(p1,p4) - 2._dp/(s(p1,p3))*s(p1,p2)**2*s(p3,p4)
     &     + 4._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4)*s(p3,p4) + 4._dp/(s(p1,p3
     &    ))*s(p1,p2)*s(p3,p4)**2 - 2._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4)**
     &    2 - 2._dp/(s(p1,p3))*s(p3,p4)**3 + 2._dp/(s(p2,p4))*s(p1,p2)**3
     &     + 2._dp/(s(p2,p4))*s(p1,p2)**2*s(p1,p3) - 6._dp/(s(p2,p4))*s(
     &    p1,p2)**2*s(p3,p4) - 4._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3)*s(p3,
     &    p4) + 6._dp/(s(p2,p4))*s(p1,p2)*s(p3,p4)**2 + 2._dp/(s(p2,p4))*
     &    s(p1,p3)*s(p3,p4)**2 - 2._dp/(s(p2,p4))*s(p3,p4)**3 )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDp1*nDp2 * ( 16._dp + 16._dp
     &    /(s(p1,p2))*s(p1,p3) - 48._dp/(s(p1,p2))*s(p1,p4) + 16._dp/(s(
     &    p1,p2))*s(p2,p4) - 8._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,p2)
     &    )/(s(p1,p2))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p1,p3)*s(p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p3
     &    ,p4) + 64._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p2,p4) + 32._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p3,p4) + 32._dp/(s(p1,p2))/(
     &    s(p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) + 32._dp/(s(p1,p2))
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 + 16._dp/(s(p1,p2)
     &    )/(s(p1,p3))*s(p1,p4)*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p1,p3))*
     &    s(p1,p4)*s(p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p3))*s(p2,p4)*s(p3,
     &    p4) + 8._dp/(s(p1,p2))/(s(p1,p3))*s(p3,p4)**2 + 32._dp/(s(p1,p2
     &    ))/(s(p1,p3))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) - 8._dp/(s(p1,p2
     &    ))/(s(p2,p4))*s(p1,p4)*s(p3,p4) + 8._dp/(s(p1,p3))*s(p1,p2) -
     &    16._dp/(s(p1,p3))*s(p1,p4) - 16._dp/(s(p1,p3))*s(p3,p4) - 8._dp
     &    /(
     & s(p1,p3))/(s(p2,p4))*s(p1,p4)*s(p3,p4) - 8._dp/(s(p1,p3))/(s(p2,
     &    p4))*s(p3,p4)**2 - 8._dp/(s(p2,p4))*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDp1*nDp3 * ( 24._dp + 48._dp
     &    /(s(p1,p2))*s(p1,p3) - 24._dp/(s(p1,p2))*s(p1,p4) + 8._dp/(s(p1
     &    ,p2))*s(p2,p4) + 32._dp/(s(p1,p2))*s(p3,p4) + 32._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p1,p3)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3
     &    )*s(p2,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4)*s(p2,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)*s(p2,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)**2
     &     - 8._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p2,p4) - 8._dp/(s(p1,
     &    p2))/(s(p1,p3))*s(p2,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p3))*s(p2
     &    ,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p3))/(s(p3,p4))*s(p1,
     &    p4)*s(p2,p4)**2 - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p2,
     &    p4) - 48._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4) - 16._dp/(
     &    s(p1,p2))/(s(p3,p4))*s(p2,p4)**2 - 24._dp/(s(p1,p3))*s(p1,p4)
     &     + 16._dp/(s(p1,p3))*s(p2,p4) + 24._dp/(s(p1,p3))*s(p3,p4) - 8.D
     &    0/(s(p1,p3))/(s(p3,p4))*s(p1,p2)*s(p2,p4) - 16._dp/(s(p1,p3))
     &    /(
     & s(p3,p4))*s(p1,p4)*s(p2,p4) - 16._dp/(s(p3,p4))*s(p2,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDp1**2 * ( 48._dp/(s(p1,p2)
     &    )*s(p2,p4) + 56._dp/(s(p1,p2))*s(p3,p4) + 32._dp/(s(p1,p2))/(s(
     &    p1,p2))*s(p1,p3)*s(p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,
     &    p3)*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p3,p4)
     &     - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p2,p4)**2 - 16._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p2,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(
     &    s(p3,p4))*s(p1,p4)*s(p2,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))
     &    /(s(p3,p4))*s(p2,p4)**3 + 16._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4
     &    )*s(p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p3))*s(p2,p4)*s(p3,p4) + 8
     &    ._dp/(s(p1,p2))/(s(p1,p3))*s(p3,p4)**2 - 16._dp/(s(p1,p2))/(s(p1
     &    ,p3))/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 + 16._dp/(s(p1,p3))*s(p2
     &    ,p4) + 24._dp/(s(p1,p3))*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDp2*nDp3 * ( 48._dp + 64._dp
     &    /(s(p1,p2))*s(p1,p3) + 56._dp/(s(p1,p2))*s(p1,p4) + 32._dp/(s(
     &    p1,p2))/(s(p1,p2))*s(p1,p3)**2 + 32._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)*s(p1,p4)**2 + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p3)*s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p3))*s(p1,
     &    p4)**2 + 8._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p2,p4) - 8._dp
     &    /(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p3,p4) + 16._dp/(s(p1,p2))/(
     &    s(p1,p3))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) + 16._dp/(s(p1,p2))
     &    /(s(p2,p4))*s(p1,p3)**2 - 16._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3
     &    )*s(p3,p4) + 8._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p4)**2 + 16._dp/(
     &    s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 48._dp/(s(p1,p2))/(s(
     &    p3,p4))*s(p1,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(
     &    p2,p4) + 16._dp/(s(p1,p3))*s(p1,p2) + 16._dp/(s(p1,p3))*s(p1,p4
     &    ) + 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p4)**2 + 8._dp/(s(p1,p3))
     &    /(
     & s(p2,p4))*s(p1,p4)*s(p3,p4) + 8._dp/(s(p1,p3))/(s(p3,p4))*s(p1,p2
     &    )*s(p1,p4) + 16._dp/(s(p1,p3))/(s(p3,p4))*s(p1,p4)**2 + 16._dp
     &    /(s(p2,p4))*s(p1,p3) + 8._dp/(s(p2,p4))*s(p1,p4) + 16._dp/(s(p3
     &    ,p4))*s(p1,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDp2**2 * (  - 16._dp/(s(p1,
     &    p2))*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)
     &     - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2 - 16._dp/(s(p1,p2))
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p4)**3 - 16._dp/(s(p1,p2))/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) - 16._dp/(s(p1,p2))/(s(p1
     &    ,p3))*s(p1,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p3,
     &    p4) - 16._dp/(s(p1,p2))/(s(p1,p3))/(s(p3,p4))*s(p1,p4)**3 - 16
     &    ._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))
     &    /(s(p2,p4))*s(p1,p3)*s(p3,p4) - 8._dp/(s(p1,p2))/(s(p2,p4))*s(
     &    p1,p4)*s(p3,p4) - 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p1,p4
     &    ) - 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4) - 16._dp/(s(
     &    p2,p4))*s(p1,p4) - 24._dp/(s(p2,p4))*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDp3**2 * (  - 32._dp/(s(p1,
     &    p2))*s(p1,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,
     &    p3)**2*s(p1,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1
     &    ,p3)**2*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)**2 + 8
     &    ._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p2,p4) - 16._dp/(s(p1,p2))
     &    /(s(p3,p4))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p3)*s(p2,p4) - 32._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p4) - 24._dp/(s(p1,p3))*s(p1,p4) - 8._dp/(s(p1,p3))/(s(p3,p4))
     &    *s(p1,p2)*s(p2,p4) - 16._dp/(s(p1,p3))/(s(p3,p4))*s(p1,p4)*s(
     &    p2,p4) - 8._dp/(s(p3,p4))*s(p1,p4) - 24._dp/(s(p3,p4))*s(p2,p4)
     &     )
      qqgghn_ab = qqgghn_ab + s123**(-1)*xn*nDn * ( 2._dp*s(p1,p2) + 2._dp
     &    *s(p1,p3) + 36._dp*s(p1,p4) + 16._dp*s(p2,p4) + 16._dp*s(p3,p4)
     &     + 36._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4) + 18._dp/(s(p1,p2))*s(p1
     &    ,p3)*s(p2,p4) + 8._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 14._dp/(s(
     &    p1,p2))*s(p1,p4)**2 - 20._dp/(s(p1,p2))*s(p1,p4)*s(p2,p4) - 26
     &    ._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4) - 6._dp/(s(p1,p2))*s(p2,p4)**2
     &     - 10._dp/(s(p1,p2))*s(p2,p4)*s(p3,p4) - 10._dp/(s(p1,p2))*s(p3
     &    ,p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4) +
     &    16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p2,p4) + 8._dp/(s(p1
     &    ,p2))/(s(p1,p2))*s(p1,p3)**2*s(p3,p4) - 8._dp/(s(p1,p2))/(s(p1
     &    ,p2))*s(p1,p3)*s(p1,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,
     &    p3)*s(p1,p4)*s(p2,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*
     &    s(p1,p4)*s(p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)**2*s(p1,p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4)
     &    )*s(p1,p3)**2*s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))
     &    /(
     & s(p3,p4))*s(p1,p3)**2*s(p2,p4)**2 + 2._dp/(s(p1,p2))/(s(p2,p4))*
     &    s(p1,p3)*s(p1,p4)**2 - 2._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)*s(
     &    p3,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2
     &     + 24._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)*s(p2,p4) + 8
     &    ._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p2,p4)**2 - 2._dp/(s(p1,p2
     &    ))/(s(p3,p4))*s(p1,p4)**3 - 6._dp/(s(p1,p2))/(s(p3,p4))*s(p1,
     &    p4)**2*s(p2,p4) - 6._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4
     &    )**2 - 2._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p4)**3 + 12._dp/(s(p1,
     &    p3))*s(p1,p2)*s(p1,p4) + 4._dp/(s(p1,p3))*s(p1,p2)*s(p2,p4) +
     &    8._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4) - 6._dp/(s(p1,p3))*s(p1,p4)
     &    **2 - 6._dp/(s(p1,p3))*s(p1,p4)*s(p2,p4) - 8._dp/(s(p1,p3))*s(
     &    p1,p4)*s(p3,p4) - 2._dp/(s(p1,p3))*s(p2,p4)*s(p3,p4) - 4._dp/(
     &    s(p1,p3))*s(p3,p4)**2 + 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*
     &    s(p1,p4)**2 + 4._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p1,p4)*s(
     &    p3,p4) + 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4)**2 + 4.D
     &    0/(
     & s(p1,p3))/(s(p3,p4))*s(p1,p2)*s(p1,p4)**2 + 4._dp/(s(p1,p3))/(s(
     &    p3,p4))*s(p1,p2)*s(p1,p4)*s(p2,p4) - 2._dp/(s(p1,p3))/(s(p3,p4
     &    ))*s(p1,p4)**3 - 4._dp/(s(p1,p3))/(s(p3,p4))*s(p1,p4)**2*s(p2,
     &    p4) - 2._dp/(s(p1,p3))/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 - 2._dp
     &    /(s(p2,p4))*s(p1,p2)**2 - 2._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3)
     &     + 2._dp/(s(p2,p4))*s(p1,p2)*s(p1,p4) + 6._dp/(s(p2,p4))*s(p1,
     &    p2)*s(p3,p4) + 2._dp/(s(p2,p4))*s(p1,p3)*s(p1,p4) + 4._dp/(s(p2
     &    ,p4))*s(p1,p3)*s(p3,p4) + 4._dp/(s(p2,p4))*s(p1,p4)**2 + 4._dp
     &    /(s(p2,p4))*s(p1,p4)*s(p3,p4) - 2._dp/(s(p2,p4))*s(p3,p4)**2
     &     + 12._dp/(s(p3,p4))*s(p1,p4)**2 + 16._dp/(s(p3,p4))*s(p1,p4)*
     &    s(p2,p4) + 4._dp/(s(p3,p4))*s(p2,p4)**2 )
      qqgghn_ab = qqgghn_ab + s124**(-2)*xn*nDn * ( 2._dp*s(p1,p3)**2 +
     &    4._dp*s(p1,p3)*s(p2,p3) + 4._dp*s(p1,p3)*s(p3,p4) + 2._dp*s(p2,
     &    p3)**2 + 4._dp*s(p2,p3)*s(p3,p4) + 2._dp*s(p3,p4)**2 + 4._dp/(s(
     &    p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4)**2 + 8._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p1,p3)*s(p1,p4)**2*s(p2,p3) + 8._dp/(s(p1,p2))/(
     &    s(p1,p2))*s(p1,p3)*s(p1,p4)**2*s(p3,p4) + 4._dp/(s(p1,p2))/(s(
     &    p1,p2))*s(p1,p4)**2*s(p2,p3)**2 + 8._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p1,p4)**2*s(p2,p3)*s(p3,p4) + 4._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p4)**2*s(p3,p4)**2 + 2._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3)**2
     &     + 4._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3)*s(p2,p3) + 4._dp/(s(p2,p4
     &    ))*s(p1,p2)*s(p1,p3)*s(p3,p4) + 2._dp/(s(p2,p4))*s(p1,p2)*s(p2
     &    ,p3)**2 + 4._dp/(s(p2,p4))*s(p1,p2)*s(p2,p3)*s(p3,p4) + 2._dp/(
     &    s(p2,p4))*s(p1,p2)*s(p3,p4)**2 )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDp1*nDp2 * (  - 16._dp/(s(
     &    p1,p2))*s(p1,p3) - 48._dp/(s(p1,p2))*s(p1,p4) - 16._dp/(s(p1,p2
     &    ))*s(p2,p3) - 32._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,p2))/(
     &    s(p1,p2))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p3)*s(p3,p4) - 64._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2 +
     &    32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p2,p3) + 32._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p4)*s(p3,p4) - 32._dp/(s(p1,p2))/(s(p1,p2
     &    ))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p2,
     &    p4))*s(p1,p3)**2 + 32._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)*s(p3,
     &    p4) + 8._dp/(s(p1,p2))/(s(p2,p4))*s(p2,p3)**2 - 8._dp/(s(p1,p2)
     &    )/(s(p2,p4))*s(p3,p4)**2 - 32._dp/(s(p1,p2))/(s(p3,p4))*s(p1,
     &    p3)*s(p1,p4) - 16._dp/(s(p2,p4))*s(p2,p3) - 16._dp/(s(p2,p4))*
     &    s(p3,p4) )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDp1*nDp3 * (  - 16._dp - 64
     &    ._dp/(s(p1,p2))*s(p1,p4) + 8._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(
     &    p1,p2))/(s(p1,p2))*s(p1,p4)**2 + 32._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p1,p4)*s(p2,p3) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p4)**2*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2
     &     - 48._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 16._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p4)*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p3)**2 -
     &    16._dp/(s(p3,p4))*s(p1,p3) + 16._dp/(s(p3,p4))*s(p1,p4) )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDp1**2 * (  - 32._dp - 64._dp
     &    /(s(p1,p2))*s(p1,p4) + 16._dp/(s(p1,p2))*s(p2,p3) + 32._dp/(s(
     &    p1,p2))*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p3,
     &    p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2 + 32._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p4)*s(p2,p3) + 32._dp/(s(p1,p2))/(s(p1,p2
     &    ))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p2,p3)*
     &    s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(
     &    p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**
     &    2*s(p2,p3) - 32._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) -
     &    16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p3) - 16._dp/(s(p3,
     &    p4))*s(p1,p3) )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDp2*nDp3 * (  - 16._dp + 32
     &    ._dp/(s(p1,p2))*s(p1,p3) + 32._dp/(s(p1,p2))*s(p1,p4) - 32._dp/(
     &    s(p1,p2))*s(p2,p3) - 24._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2
     &    ))*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,
     &    p3)*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1
     &    ,p4)**2*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2 -
     &    16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 16._dp/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p4)*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p3)**2 - 16
     &    ._dp/(s(p2,p4))*s(p1,p2) + 32._dp/(s(p2,p4))*s(p1,p3) - 32._dp/(
     &    s(p2,p4))*s(p2,p3) - 24._dp/(s(p2,p4))*s(p3,p4) - 16._dp/(s(p2,
     &    p4))/(s(p3,p4))*s(p1,p2)*s(p2,p3) - 8._dp/(s(p2,p4))/(s(p3,p4)
     &    )*s(p1,p3)**2 - 8._dp/(s(p2,p4))/(s(p3,p4))*s(p2,p3)**2 - 16._dp
     &    /(s(p3,p4))*s(p2,p3) )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDp2**2 * (  - 16._dp + 16._dp
     &    /(s(p1,p2))*s(p1,p3) + 16._dp/(s(p1,p2))*s(p1,p4) - 32._dp/(s(
     &    p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1
     &    ,p2))*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p4)**2*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3)**2
     &     + 8._dp/(s(p1,p2))/(s(p2,p4))*s(p2,p3)**2 + 16._dp/(s(p1,p2))
     &    /(s(p2,p4))*s(p2,p3)*s(p3,p4) + 8._dp/(s(p1,p2))/(s(p2,p4))*s(
     &    p3,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p3) -
     &    16._dp/(s(p2,p4))*s(p1,p2) + 16._dp/(s(p2,p4))*s(p1,p3) + 8._dp
     &    /(s(p2,p4))/(s(p2,p4))*s(p1,p3)**2 + 8._dp/(s(p2,p4))/(s(p2,p4
     &    ))*s(p2,p3)**2 + 16._dp/(s(p2,p4))/(s(p2,p4))*s(p2,p3)*s(p3,p4
     &    ) + 8._dp/(s(p2,p4))/(s(p2,p4))*s(p3,p4)**2 - 16._dp/(s(p2,p4))
     &    /(s(p3,p4))*s(p1,p2)*s(p2,p3) - 16._dp/(s(p3,p4))*s(p2,p3) )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDp3**2 * (  - 16._dp/(s(p1,
     &    p2))*s(p1,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,
     &    p3)*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1
     &    ,p4)**2*s(p2,p3) + 32._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2 -
     &    16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p3) + 16._dp/(s(p3,
     &    p4))*s(p1,p4) )
      qqgghn_ab = qqgghn_ab + s124**(-1)*xn*nDn * (  - 6._dp*s(p1,p2) +
     &    16._dp*s(p1,p3) + 8._dp*s(p1,p4) + 16._dp*s(p2,p3) + 20._dp*s(p3,
     &    p4) - 4._dp/(s(p1,p2))*s(p1,p3)**2 + 16._dp/(s(p1,p2))*s(p1,p3)
     &    *s(p1,p4) - 8._dp/(s(p1,p2))*s(p1,p3)*s(p2,p3) - 10._dp/(s(p1,
     &    p2))*s(p1,p3)*s(p3,p4) - 8._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4) -
     &    4._dp/(s(p1,p2))*s(p2,p3)**2 - 10._dp/(s(p1,p2))*s(p2,p3)*s(p3,
     &    p4) - 10._dp/(s(p1,p2))*s(p3,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p2
     &    ))*s(p1,p3)**2*s(p1,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3
     &    )*s(p1,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*
     &    s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3
     &    ,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2*s(p2,p3) + 8._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2*s(p3,p4) + 8._dp/(s(p1,p2))
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p1,p4)**2 + 16._dp/(s(p1,
     &    p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2*s(p2,p3) + 8.D
     &    0/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p3)**2 +
     &    2._dp/(
     & s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p3,p4)**2 + 8._dp/(s(p1,p2))/(s(
     &    p3,p4))*s(p1,p3)**2*s(p1,p4) + 8._dp/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)*s(p1,p4)*s(p2,p3) + 2._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4)
     &     - 4._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4) - 4._dp/(s(p2,p4))*s(p1,
     &    p2)**2 + 4._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3) + 12._dp/(s(p2,p4))
     &    *s(p1,p2)*s(p2,p3) + 12._dp/(s(p2,p4))*s(p1,p2)*s(p3,p4) - 2._dp
     &    /(s(p2,p4))*s(p1,p3)**2 - 8._dp/(s(p2,p4))*s(p1,p3)*s(p2,p3)
     &     - 4._dp/(s(p2,p4))*s(p1,p3)*s(p3,p4) - 6._dp/(s(p2,p4))*s(p2,
     &    p3)**2 - 10._dp/(s(p2,p4))*s(p2,p3)*s(p3,p4) - 6._dp/(s(p2,p4))
     &    *s(p3,p4)**2 - 2._dp/(s(p2,p4))/(s(p3,p4))*s(p1,p2)**2*s(p2,p3
     &    ) + 4._dp/(s(p2,p4))/(s(p3,p4))*s(p1,p2)*s(p1,p3)*s(p2,p3) + 4
     &    ._dp/(s(p2,p4))/(s(p3,p4))*s(p1,p2)*s(p2,p3)**2 - 2._dp/(s(p3,p4
     &    ))*s(p1,p2)*s(p1,p3) - 2._dp/(s(p3,p4))*s(p1,p2)*s(p2,p3) + 4.D
     &    0/(s(p3,p4))*s(p1,p3)**2 + 8._dp/(s(p3,p4))*s(p1,p3)*s(p2,p3)
     &     + 4._dp/(s(p3,p4))*s(p2,p3)**2 )
      qqgghn_ab = qqgghn_ab + xn*nDp1*nDp2 * (  - 8._dp/(s(p1,p2)) + 32.D
     &    0/(s(p1,p2))/(s(p1,p2))*s(p1,p3) + 96._dp/(s(p1,p2))/(s(p1,p2)
     &    )*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(
     &    p2,p3) + 64._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p4) + 16._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p3) + 8._dp/(s(p1,p2))
     &    /(s(p2,p4))*s(p1,p4) + 16._dp/(s(p1,p2))/(s(p2,p4))*s(p2,p3)
     &     + 16._dp/(s(p1,p2))/(s(p2,p4))*s(p3,p4) - 8._dp/(s(p1,p3)) - 8
     &    ._dp/(s(p1,p3))/(s(p2,p4))*s(p3,p4) )
      qqgghn_ab = qqgghn_ab + xn*nDp1*nDp3 * ( 16._dp/(s(p1,p2)) - 32._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p3) + 32._dp/(s(p1,p2))/(s(p1,p2))
     &    *s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(
     &    p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p3) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,
     &    p4) - 24._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4) + 8._dp/(s(p1,p2))
     &    /(s(p3,p4))*s(p2,p4) + 8._dp/(s(p1,p3)) + 8._dp/(s(p1,p3))/(s(
     &    p3,p4))*s(p2,p4) )
      qqgghn_ab = qqgghn_ab + xn*nDp1**2 * ( 32._dp/(s(p1,p2)) + 32._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p2,p3) - 64._dp/(s(p1,p2))/(s(p1,p2))*s(p2,p4) - 32._dp/(s(p1
     &    ,p2))/(s(p1,p2))*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3
     &    ,p4))*s(p1,p3)*s(p1,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,
     &    p4))*s(p1,p3)*s(p2,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4
     &    ))*s(p1,p4)*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))
     &    *s(p2,p3)*s(p2,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*
     &    s(p2,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3) )
      qqgghn_ab = qqgghn_ab + xn*nDp2*nDp3 * (  - 32._dp/(s(p1,p2)) + 32
     &    ._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 32._dp
     &    /(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2 - 32._dp/(s(p1,
     &    p2))/(s(p2,p4))*s(p1,p3) - 16._dp/(s(p1,p2))/(s(p2,p4))*s(p1,
     &    p4) + 8._dp/(s(p1,p2))/(s(p2,p4))*s(p3,p4) - 8._dp/(s(p1,p2))/(
     &    s(p2,p4))/(s(p3,p4))*s(p1,p3)**2 - 16._dp/(s(p1,p2))/(s(p2,p4)
     &    )/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 8._dp/(s(p1,p2))/(s(p2,p4))/(
     &    s(p3,p4))*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p2,p4))/(s(p3,p4)
     &    )*s(p1,p4)*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p2,p4))/(s(p3,p4))*
     &    s(p2,p3)**2 - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4) - 32._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p3,p4))*s(
     &    p2,p4) + 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p4) - 8._dp/(s(p1,p3)
     &    )/(s(p3,p4))*s(p1,p4) + 16._dp/(s(p2,p4)) + 16._dp/(s(p2,p4))/(
     &    s(p3,p4))*s(p2,p3) )
      qqgghn_ab = qqgghn_ab + xn*nDp2**2 * ( 8._dp/(s(p1,p2)) - 32._dp/(
     &    s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2 - 16._dp/(s(p1,p2)
     &    )/(s(p2,p4))*s(p1,p3) - 8._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p4)
     &     + 16._dp/(s(p1,p2))/(s(p2,p4))*s(p2,p3) + 16._dp/(s(p1,p2))/(
     &    s(p2,p4))*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p2,p4))/(s(p3,p4))*
     &    s(p1,p4)*s(p2,p3) + 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p4) + 16.D
     &    0/(s(p2,p4)) + 16._dp/(s(p2,p4))/(s(p3,p4))*s(p2,p3) )
      qqgghn_ab = qqgghn_ab + xn*nDp3**2 * ( 8._dp/(s(p1,p2)) + 32._dp/(
     &    s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 16._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p4) + 16._dp/(s(p1,p2))/(s(p3,p4))*s(
     &    p2,p3) + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p4) + 8._dp/(s(p1,p2
     &    ))/(s(p3,p4))/(s(p3,p4))*s(p1,p3)**2 + 16._dp/(s(p1,p2))/(s(p3
     &    ,p4))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 8._dp/(s(p1,p2))/(s(p3,p4
     &    ))/(s(p3,p4))*s(p1,p4)**2 + 8._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,
     &    p4))*s(p2,p3)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(
     &    p2,p3)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p2,
     &    p4)**2 - 8._dp/(s(p1,p3))/(s(p3,p4))*s(p1,p4) )
      qqgghn_ab = qqgghn_ab + xn*nDn * ( 6._dp - 20._dp/(s(p1,p2))*s(p1,
     &    p3) - 42._dp/(s(p1,p2))*s(p1,p4) - 16._dp/(s(p1,p2))*s(p2,p3)
     &     - 24._dp/(s(p1,p2))*s(p2,p4) - 24._dp/(s(p1,p2))*s(p3,p4) + 4.D
     &    0/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2 - 32._dp/(s(p1,p2))/(s(p1,
     &    p2))*s(p1,p3)*s(p1,p4) + 4._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)
     &    **2 - 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p1,
     &    p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p2,
     &    p4) - 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)
     &    **2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4
     &    )*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p4)**2*s(p2,p3) + 4._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))/(
     &    s(p3,p4))*s(p1,p3)**2*s(p2,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p2)
     &    )/(s(p3,p4))/(s(p3,p4))*s(p1,p3)*s(p1,p4)*s(p2,p3)*s(p2,p4)
     &     + 4._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p1,p4)
     &    **2*s(p2,p3)**2 )
      qqgghn_ab = qqgghn_ab + xn*nDn * (  - 2._dp/(s(p1,p2))/(s(p1,p3))*
     &    s(p1,p4)**2 - 2._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p2,p4) +
     &    2._dp/(s(p1,p2))/(s(p1,p3))*s(p1,p4)*s(p3,p4) - 2._dp/(s(p1,p2)
     &    )/(s(p1,p3))/(s(p3,p4))*s(p1,p4)**3 - 4._dp/(s(p1,p2))/(s(p1,
     &    p3))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) - 2._dp/(s(p1,p2))/(s(p1,
     &    p3))/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 + 4._dp/(s(p1,p2))/(s(p2,
     &    p4))*s(p1,p3)*s(p2,p3) - 2._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p4)
     &    **2 + 4._dp/(s(p1,p2))/(s(p2,p4))*s(p1,p4)*s(p2,p3) + 4._dp/(s(
     &    p1,p2))/(s(p2,p4))*s(p2,p3)**2 + 2._dp/(s(p1,p2))/(s(p2,p4))*
     &    s(p2,p3)*s(p3,p4) + 2._dp/(s(p1,p2))/(s(p2,p4))*s(p3,p4)**2 +
     &    4._dp/(s(p1,p2))/(s(p2,p4))/(s(p3,p4))*s(p1,p3)*s(p1,p4)*s(p2,
     &    p3) + 2._dp/(s(p1,p2))/(s(p2,p4))/(s(p3,p4))*s(p1,p4)**2*s(p2,
     &    p3) + 4._dp/(s(p1,p2))/(s(p2,p4))/(s(p3,p4))*s(p1,p4)*s(p2,p3)
     &    **2 - 4._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2 - 10._dp/(s(p1,p2
     &    ))/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 4._dp/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p3)*s(p2,p3) )
      qqgghn_ab = qqgghn_ab + xn*nDn * (  - 10._dp/(s(p1,p2))/(s(p3,p4))
     &    *s(p1,p3)*s(p2,p4) - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2
     &     - 6._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p3) - 24._dp/(s(p1
     &    ,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p3,p4
     &    ))*s(p2,p3)*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p4)**2
     &     - 8._dp/(s(p1,p3))*s(p1,p4) - 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1
     &    ,p4)**2 - 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p4)*s(p3,p4) - 4._dp
     &    /(s(p1,p3))/(s(p3,p4))*s(p1,p4)**2 - 4._dp/(s(p1,p3))/(s(p3,p4
     &    ))*s(p1,p4)*s(p2,p4) + 4._dp/(s(p2,p4))*s(p1,p2) - 4._dp/(s(p2,
     &    p4))*s(p1,p3) - 4._dp/(s(p2,p4))*s(p1,p4) - 12._dp/(s(p2,p4))*
     &    s(p2,p3) - 8._dp/(s(p2,p4))*s(p3,p4) + 2._dp/(s(p2,p4))/(s(p3,
     &    p4))*s(p1,p2)*s(p2,p3) - 4._dp/(s(p2,p4))/(s(p3,p4))*s(p1,p3)*
     &    s(p2,p3) - 2._dp/(s(p2,p4))/(s(p3,p4))*s(p1,p4)*s(p2,p3) - 4._dp
     &    /(s(p2,p4))/(s(p3,p4))*s(p2,p3)**2 + 2._dp/(s(p3,p4))*s(p1,p3)
     &     )

      qqgghn_ba=  + s123**(-2)*xn*nDn * ( 2._dp*s(p1,p4)**2 + 4._dp*s(p1,
     &    p4)*s(p2,p4) + 4._dp*s(p1,p4)*s(p3,p4) + 2._dp*s(p2,p4)**2 + 4.D
     &    0*s(p2,p4)*s(p3,p4) + 2._dp*s(p3,p4)**2 + 4._dp/(s(p1,p2))/(s(
     &    p1,p2))*s(p1,p3)**2*s(p1,p4)**2 + 8._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p1,p3)**2*s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p3)**2*s(p1,p4)*s(p3,p4) + 4._dp/(s(p1,p2))/(s(p1,p2))*s(p1
     &    ,p3)**2*s(p2,p4)**2 + 8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*
     &    s(p2,p4)*s(p3,p4) + 4._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(
     &    p3,p4)**2 + 2._dp/(s(p2,p3))*s(p1,p2)*s(p1,p4)**2 + 4._dp/(s(p2
     &    ,p3))*s(p1,p2)*s(p1,p4)*s(p2,p4) + 4._dp/(s(p2,p3))*s(p1,p2)*
     &    s(p1,p4)*s(p3,p4) + 2._dp/(s(p2,p3))*s(p1,p2)*s(p2,p4)**2 + 4.D
     &    0/(s(p2,p3))*s(p1,p2)*s(p2,p4)*s(p3,p4) + 2._dp/(s(p2,p3))*s(
     &    p1,p2)*s(p3,p4)**2 )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDp1*nDp2 * (
     &     - 16._dp*s(p3,p4) - 80._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 48._dp
     &    /(s(p1,p2))*s(p1,p4)*s(p3,p4) - 32._dp/(s(p1,p2))*s(p3,p4)**2
     &     - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p3,p4) - 64._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4) - 32._dp/(s(p1
     &    ,p2))/(s(p1,p2))*s(p1,p4)**2*s(p3,p4) - 16._dp/(s(p1,p2))/(s(
     &    p1,p4))*s(p1,p3)**2*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p4))*s(
     &    p1,p3)*s(p3,p4)**2 - 16._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)**2*
     &    s(p3,p4) + 16._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)*s(p3,p4)**2
     &     - 16._dp/(s(p1,p4))*s(p1,p3)*s(p3,p4) - 16._dp/(s(p1,p4))*s(p3
     &    ,p4)**2 - 16._dp/(s(p2,p3))*s(p1,p4)*s(p3,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDp1*nDp3 * ( 16
     &    ._dp*s(p1,p2) - 32._dp*s(p1,p3) + 32._dp*s(p1,p4) - 48._dp*s(p3,p4
     &    ) - 16._dp/(s(p1,p2))*s(p1,p3)**2 - 16._dp/(s(p1,p2))*s(p1,p3)*
     &    s(p1,p4) - 16._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) + 16._dp/(s(p1,
     &    p2))*s(p1,p4)**2 - 32._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4) - 16._dp
     &    /(s(p1,p2))/(s(p2,p3))*s(p1,p4)**3 + 16._dp/(s(p1,p2))/(s(p2,
     &    p3))*s(p1,p4)**2*s(p3,p4) - 16._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3
     &    ) - 16._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4) - 16._dp/(s(p1,p4))*s(
     &    p1,p3)**2 - 16._dp/(s(p1,p4))*s(p1,p3)*s(p3,p4) - 16._dp/(s(p2,
     &    p3))*s(p1,p2)*s(p1,p4) - 32._dp/(s(p2,p3))*s(p1,p4)**2 + 16._dp
     &    /(s(p2,p3))*s(p1,p4)*s(p3,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDp1**2 * (  -
     &    48._dp*s(p3,p4) - 64._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 48._dp/(
     &    s(p1,p2))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p3)**2*s(p3,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(
     &    p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2*s(
     &    p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3)**2*s(p3,p4) -
     &    16._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4) - 32._dp/(s(p1,p4))*s(p1,p3
     &    )*s(p3,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDp2*nDp3 * (
     &     - 16._dp*s(p1,p3) + 16._dp*s(p1,p4) + 16._dp*s(p3,p4) - 16._dp/(
     &    s(p1,p2))*s(p1,p3)**2 - 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4) +
     &    16._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) + 16._dp/(s(p1,p2))*s(p1,p4
     &    )**2 + 32._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))
     &    /(s(p2,p3))*s(p1,p4)**3 - 16._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4
     &    )**2*s(p3,p4) - 16._dp/(s(p2,p3))*s(p1,p4)**2 )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDp2**2 * (  -
     &    16._dp/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,
     &    p2))*s(p1,p3)**2*s(p3,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,
     &    p3)*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)
     &    **2*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)**2*s(p3,
     &    p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDp3**2 * ( 16._dp
     &    *s(p1,p2) + 16._dp*s(p1,p3) + 48._dp*s(p1,p4) + 16._dp/(s(p1,p2)
     &    )*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))*s(p1,p4)**2 - 16._dp/(
     &    s(p1,p2))/(s(p2,p3))*s(p1,p4)**3 - 16._dp/(s(p2,p3))*s(p1,p4)
     &    **2 )
      qqgghn_ba = qqgghn_ba + s123**(-1)*s124**(-1)*xn*nDn * (  - 8._dp*
     &    s(p1,p2)*s(p1,p3) - 4._dp*s(p1,p2)*s(p3,p4) - 8._dp*s(p1,p3)*s(
     &    p1,p4) + 16._dp*s(p1,p3)*s(p3,p4) + 8._dp*s(p3,p4)**2 + 16._dp/(
     &    s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4) - 8._dp/(s(p1,p2))*s(p1,
     &    p3)*s(p3,p4)**2 - 4._dp/(s(p1,p2))*s(p3,p4)**3 - 8._dp/(s(p1,p2
     &    ))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4)**2 - 2._dp/(s(p1,p4))
     &    *s(p1,p2)**2*s(p1,p3) - 2._dp/(s(p1,p4))*s(p1,p2)**2*s(p3,p4)
     &     + 4._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3)*s(p3,p4) + 4._dp/(s(p1,p4
     &    ))*s(p1,p2)*s(p3,p4)**2 - 2._dp/(s(p1,p4))*s(p1,p3)*s(p3,p4)**
     &    2 - 2._dp/(s(p1,p4))*s(p3,p4)**3 + 2._dp/(s(p2,p3))*s(p1,p2)**3
     &     + 2._dp/(s(p2,p3))*s(p1,p2)**2*s(p1,p4) - 6._dp/(s(p2,p3))*s(
     &    p1,p2)**2*s(p3,p4) - 4._dp/(s(p2,p3))*s(p1,p2)*s(p1,p4)*s(p3,
     &    p4) + 6._dp/(s(p2,p3))*s(p1,p2)*s(p3,p4)**2 + 2._dp/(s(p2,p3))*
     &    s(p1,p4)*s(p3,p4)**2 - 2._dp/(s(p2,p3))*s(p3,p4)**3 )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDp1*nDp2 * (  - 16._dp/(s(
     &    p1,p2))*s(p1,p3) - 16._dp/(s(p1,p2))*s(p1,p4) - 16._dp/(s(p1,p2
     &    ))*s(p2,p4) - 24._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,p2))/(
     &    s(p1,p2))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p3)*s(p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p3,
     &    p4) + 64._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p2,p4) + 32._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p3,p4) + 32._dp/(s(p1,p2))/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) + 32._dp/(s(p1,p2))/(
     &    s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 - 8._dp/(s(p1,p2))/(
     &    s(p1,p4))*s(p2,p4)*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p2,p3))*s(
     &    p1,p4)*s(p2,p4) + 24._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)*s(p3,
     &    p4) - 24._dp/(s(p1,p2))/(s(p2,p3))*s(p2,p4)*s(p3,p4) - 8._dp/(
     &    s(p1,p2))/(s(p2,p3))*s(p3,p4)**2 + 32._dp/(s(p1,p2))/(s(p2,p3)
     &    )/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 - 8._dp/(s(p1,p4))*s(p3,p4)
     &     - 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4)*s(p3,p4) - 8._dp/(s(p1,
     &    p4))/(
     & s(p2,p3))*s(p3,p4)**2 + 8._dp/(s(p2,p3))*s(p1,p2) - 16._dp/(s(p2,
     &    p3))*s(p2,p4) - 16._dp/(s(p2,p3))*s(p3,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDp1*nDp3 * ( 16._dp/(s(p1,
     &    p2))*s(p1,p3) - 16._dp/(s(p1,p2))*s(p1,p4) + 40._dp/(s(p1,p2))*
     &    s(p2,p4) + 32._dp/(s(p1,p2))*s(p3,p4) + 32._dp/(s(p1,p2))/(s(p1
     &    ,p2))*s(p1,p3)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p2
     &    ,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,
     &    p4)*s(p2,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3
     &    )*s(p2,p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3)**2 + 16.D
     &    0/(s(p1,p2))/(s(p1,p4))*s(p1,p3)*s(p3,p4) + 8._dp/(s(p1,p2))/(
     &    s(p1,p4))*s(p2,p4)**2 + 16._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)
     &    **2 - 8._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)*s(p2,p4) - 16._dp/(
     &    s(p1,p2))/(s(p2,p3))*s(p1,p4)*s(p3,p4) + 24._dp/(s(p1,p2))/(s(
     &    p2,p3))*s(p2,p4)**2 + 8._dp/(s(p1,p2))/(s(p2,p3))*s(p2,p4)*s(
     &    p3,p4) + 16._dp/(s(p1,p2))/(s(p2,p3))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p4)**2 - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p2,p4) - 16.D
     &    0/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4) + 16._dp/(s(p1,p2))
     &    /(
     & s(p3,p4))*s(p2,p4)**2 + 16._dp/(s(p1,p4))*s(p1,p3) + 8._dp/(s(p1,
     &    p4))*s(p2,p4) + 16._dp/(s(p1,p4))*s(p3,p4) + 8._dp/(s(p1,p4))/(
     &    s(p2,p3))*s(p2,p4)**2 + 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4)*
     &    s(p3,p4) + 16._dp/(s(p2,p3))*s(p1,p2) + 16._dp/(s(p2,p3))*s(p1,
     &    p4) + 16._dp/(s(p2,p3))*s(p2,p4) + 8._dp/(s(p2,p3))/(s(p3,p4))*
     &    s(p1,p2)*s(p2,p4) + 16._dp/(s(p2,p3))/(s(p3,p4))*s(p2,p4)**2 )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDp1**2 * ( 16._dp/(s(p1,p2)
     &    )*s(p2,p4) + 32._dp/(s(p1,p2))*s(p3,p4) + 32._dp/(s(p1,p2))/(s(
     &    p1,p2))*s(p1,p3)*s(p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,
     &    p3)*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p3,p4)
     &     - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p2,p4)**2 - 16._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p2,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(
     &    s(p3,p4))*s(p1,p4)*s(p2,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))
     &    /(s(p3,p4))*s(p2,p4)**3 + 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3
     &    )*s(p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3)*s(p3,p4) -
     &    8._dp/(s(p1,p2))/(s(p1,p4))*s(p2,p4)*s(p3,p4) - 16._dp/(s(p1,p2
     &    ))/(s(p2,p3))*s(p2,p4)**2 - 8._dp/(s(p1,p2))/(s(p2,p3))*s(p2,
     &    p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p2,p3))/(s(p3,p4))*s(p2,p4
     &    )**3 + 8._dp/(s(p1,p4))*s(p3,p4) - 8._dp/(s(p1,p4))/(s(p2,p3))*
     &    s(p1,p2)*s(p2,p4) - 8._dp/(s(p1,p4))/(s(p2,p3))*s(p1,p2)*s(p3,
     &    p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDp2*nDp3 * ( 8._dp + 32._dp
     &    /(s(p1,p2))*s(p1,p3) + 24._dp/(s(p1,p2))*s(p1,p4) - 8._dp/(s(p1
     &    ,p2))*s(p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2 + 32
     &    ._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2 + 32._dp/(s(p1,p2)
     &    )/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)*s(p2,p4) + 8._dp/(s(
     &    p1,p2))/(s(p2,p3))*s(p1,p4)**2 - 24._dp/(s(p1,p2))/(s(p2,p3))*
     &    s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)*s(p3,
     &    p4) - 16._dp/(s(p1,p2))/(s(p2,p3))*s(p2,p4)*s(p3,p4) - 16._dp/(
     &    s(p1,p2))/(s(p2,p3))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) + 16._dp
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 16._dp/(s(p1,p2))/(
     &    s(p3,p4))*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*
     &    s(p2,p4) + 16._dp/(s(p2,p3))*s(p1,p4) - 8._dp/(s(p2,p3))*s(p2,
     &    p4) + 8._dp/(s(p2,p3))*s(p3,p4) - 8._dp/(s(p2,p3))/(s(p3,p4))*
     &    s(p1,p2)*s(p1,p4) - 16._dp/(s(p2,p3))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDp2**2 * ( 16._dp/(s(p1,p2)
     &    )*s(p1,p4) + 8._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,p2))/(s(
     &    p1,p2))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,
     &    p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**3
     &     - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4
     &    ) + 24._dp/(s(p1,p2))/(s(p2,p3))*s(p1,p4)*s(p3,p4) + 8._dp/(s(
     &    p1,p2))/(s(p2,p3))*s(p3,p4)**2 - 16._dp/(s(p1,p2))/(s(p2,p3))
     &    /(s(p3,p4))*s(p1,p4)**2*s(p2,p4) + 16._dp/(s(p2,p3))*s(p1,p4)
     &     + 8._dp/(s(p2,p3))*s(p3,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDp3**2 * (  - 16._dp - 16._dp
     &    /(s(p1,p2))*s(p1,p3) - 32._dp/(s(p1,p2))*s(p1,p4) - 16._dp/(s(
     &    p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p1,p4) - 16._dp/(
     &    s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p2,p4) + 16._dp
     &    /(s(p1,p2))/(s(p2,p3))*s(p1,p4)**2 - 8._dp/(s(p1,p2))/(s(p2,p3
     &    ))*s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p2,p3))*s(p2,p4)**2
     &     - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p2,p4) - 32._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4) - 8._dp/(s(p2,p3))*s(p2,
     &    p4) - 8._dp/(s(p2,p3))/(s(p3,p4))*s(p1,p2)*s(p1,p4) - 16._dp/(
     &    s(p2,p3))/(s(p3,p4))*s(p1,p4)*s(p2,p4) - 8._dp/(s(p3,p4))*s(p1
     &    ,p4) - 8._dp/(s(p3,p4))*s(p2,p4) )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDn * ( 6._dp*s(p1,p3) + 14.D
     &    0*s(p1,p4) + 16._dp*s(p2,p4) + 16._dp*s(p3,p4) + 14._dp/(s(p1,p2
     &    ))*s(p1,p3)*s(p1,p4) - 4._dp/(s(p1,p2))*s(p1,p3)*s(p2,p4) - 8.D
     &    0/(s(p1,p2))*s(p1,p3)*s(p3,p4) - 6._dp/(s(p1,p2))*s(p1,p4)**2
     &     - 12._dp/(s(p1,p2))*s(p1,p4)*s(p2,p4) - 10._dp/(s(p1,p2))*s(p1
     &    ,p4)*s(p3,p4) - 6._dp/(s(p1,p2))*s(p2,p4)**2 - 10._dp/(s(p1,p2)
     &    )*s(p2,p4)*s(p3,p4) - 10._dp/(s(p1,p2))*s(p3,p4)**2 + 16._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4) + 16._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p1,p3)**2*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))
     &    *s(p1,p3)**2*s(p3,p4) - 8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*
     &    s(p1,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(
     &    p2,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,
     &    p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p1,
     &    p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*
     &    s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p3)**2*s(p2,p4)**2 )
      qqgghn_ba = qqgghn_ba + s123**(-1)*xn*nDn * (  - 2._dp/(s(p1,p2))
     &    /(s(p1,p4))*s(p1,p3)*s(p2,p4)**2 + 2._dp/(s(p1,p2))/(s(p1,p4))
     &    *s(p1,p3)*s(p3,p4)**2 + 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4)**2 + 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)*s(
     &    p2,p4) - 2._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**3 - 6._dp/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) - 6._dp/(s(p1,p2))/(s(p3,
     &    p4))*s(p1,p4)*s(p2,p4)**2 - 2._dp/(s(p1,p2))/(s(p3,p4))*s(p2,
     &    p4)**3 + 2._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3) + 2._dp/(s(p1,p4))*
     &    s(p1,p2)*s(p3,p4) - 2._dp/(s(p1,p4))*s(p1,p3)*s(p2,p4) - 4._dp
     &    /(s(p1,p4))*s(p1,p3)*s(p3,p4) + 2._dp/(s(p1,p4))*s(p2,p4)**2
     &     + 4._dp/(s(p1,p4))*s(p2,p4)*s(p3,p4) + 2._dp/(s(p1,p4))/(s(p2,
     &    p3))*s(p1,p2)*s(p2,p4)**2 + 4._dp/(s(p1,p4))/(s(p2,p3))*s(p1,
     &    p2)*s(p2,p4)*s(p3,p4) + 2._dp/(s(p1,p4))/(s(p2,p3))*s(p1,p2)*
     &    s(p3,p4)**2 - 2._dp/(s(p2,p3))*s(p1,p2)**2 + 4._dp/(s(p2,p3))*
     &    s(p1,p2)*s(p1,p4) + 12._dp/(s(p2,p3))*s(p1,p2)*s(p2,p4) + 12._dp
     &    /(
     & s(p2,p3))*s(p1,p2)*s(p3,p4) - 6._dp/(s(p2,p3))*s(p1,p4)*s(p2,p4)
     &     - 2._dp/(s(p2,p3))*s(p1,p4)*s(p3,p4) - 6._dp/(s(p2,p3))*s(p2,
     &    p4)**2 - 8._dp/(s(p2,p3))*s(p2,p4)*s(p3,p4) - 6._dp/(s(p2,p3))*
     &    s(p3,p4)**2 + 4._dp/(s(p2,p3))/(s(p3,p4))*s(p1,p2)*s(p1,p4)*s(
     &    p2,p4) + 4._dp/(s(p2,p3))/(s(p3,p4))*s(p1,p2)*s(p2,p4)**2 - 2.D
     &    0/(s(p2,p3))/(s(p3,p4))*s(p1,p4)**2*s(p2,p4) - 4._dp/(s(p2,p3)
     &    )/(s(p3,p4))*s(p1,p4)*s(p2,p4)**2 - 2._dp/(s(p2,p3))/(s(p3,p4)
     &    )*s(p2,p4)**3 + 4._dp/(s(p3,p4))*s(p1,p4)**2 + 8._dp/(s(p3,p4))
     &    *s(p1,p4)*s(p2,p4) + 4._dp/(s(p3,p4))*s(p2,p4)**2 )
      qqgghn_ba = qqgghn_ba + s124**(-2)*xn*nDn * ( 6._dp*s(p1,p3)**2 +
     &    12._dp*s(p1,p3)*s(p2,p3) + 12._dp*s(p1,p3)*s(p3,p4) + 6._dp*s(p2
     &    ,p3)**2 + 12._dp*s(p2,p3)*s(p3,p4) + 6._dp*s(p3,p4)**2 + 8._dp/(
     &    s(p1,p2))*s(p1,p3)**2*s(p1,p4) + 16._dp/(s(p1,p2))*s(p1,p3)*s(
     &    p1,p4)*s(p2,p3) + 16._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4)*s(p3,p4)
     &     + 8._dp/(s(p1,p2))*s(p1,p4)*s(p2,p3)**2 + 16._dp/(s(p1,p2))*s(
     &    p1,p4)*s(p2,p3)*s(p3,p4) + 8._dp/(s(p1,p2))*s(p1,p4)*s(p3,p4)
     &    **2 + 4._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4)**2 + 8.D
     &    0/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)**2*s(p2,p3) + 8._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)**2*s(p3,p4) + 4._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p4)**2*s(p2,p3)**2 + 8._dp/(s(p1,p2)
     &    )/(s(p1,p2))*s(p1,p4)**2*s(p2,p3)*s(p3,p4) + 4._dp/(s(p1,p2))
     &    /(s(p1,p2))*s(p1,p4)**2*s(p3,p4)**2 + 2._dp/(s(p1,p4))*s(p1,p2
     &    )*s(p1,p3)**2 + 4._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3)*s(p2,p3) +
     &    4._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3)*s(p3,p4) + 2._dp/(s(p1,p4))*
     &    s(p1,p2)*s(p2,p3)**2 )
      qqgghn_ba = qqgghn_ba + s124**(-2)*xn*nDn * ( 4._dp/(s(p1,p4))*s(
     &    p1,p2)*s(p2,p3)*s(p3,p4) + 2._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4)
     &    **2 )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp1*nDp2 * (  - 16._dp - 48
     &    ._dp/(s(p1,p2))*s(p1,p3) - 80._dp/(s(p1,p2))*s(p1,p4) + 16._dp/(
     &    s(p1,p2))*s(p2,p3) - 16._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2
     &    ))*s(p1,p3)*s(p3,p4) - 64._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**
     &    2 + 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)*s(p2,p3) + 32._dp/(s(
     &    p1,p2))/(s(p1,p2))*s(p1,p4)*s(p3,p4) - 32._dp/(s(p1,p2))/(s(p1
     &    ,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p3) + 8._dp/(s(p1,p2))/(s(
     &    p1,p4))*s(p1,p3)**2 + 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3)*s(
     &    p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p4))*s(p2,p3)**2 + 16._dp/(s(p1
     &    ,p2))/(s(p1,p4))*s(p2,p3)*s(p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p4
     &    ))*s(p3,p4)**2 - 64._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4
     &    ) - 32._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p3) - 16._dp/(s(
     &    p1,p4))*s(p1,p3) - 16._dp/(s(p1,p4))*s(p3,p4) - 32._dp/(s(p3,p4
     &    ))*s(p1,p3) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp1*nDp3 * (  - 64._dp - 16
     &    ._dp/(s(p1,p2))*s(p1,p3) - 80._dp/(s(p1,p2))*s(p1,p4) + 48._dp/(
     &    s(p1,p2))*s(p2,p3) - 8._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p4)**2 + 32._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p4)*s(p2,p3) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1
     &    ,p3)*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p4)**2*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2
     &     - 80._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 16._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p4)**2 - 48._dp/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p4)*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p3)**2 -
     &    16._dp/(s(p1,p4))*s(p1,p2) - 16._dp/(s(p1,p4))*s(p1,p3) + 16._dp
     &    /(s(p1,p4))*s(p2,p3) - 8._dp/(s(p1,p4))*s(p3,p4) - 16._dp/(s(p1
     &    ,p4))/(s(p3,p4))*s(p1,p2)*s(p1,p3) - 8._dp/(s(p1,p4))/(s(p3,p4
     &    ))*s(p1,p3)**2 - 8._dp/(s(p1,p4))/(s(p3,p4))*s(p2,p3)**2 + 16.D
     &    0/(s(p3,p4))*s(p1,p2) - 64._dp/(s(p3,p4))*s(p1,p3) + 32._dp/(s(
     &    p3,p4))*s(p1,p4) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp1*nDp3 * (  - 16._dp/(s(
     &    p3,p4))*s(p2,p3) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp1**2 * (  - 64._dp - 80._dp
     &    /(s(p1,p2))*s(p1,p4) + 48._dp/(s(p1,p2))*s(p2,p3) + 48._dp/(s(
     &    p1,p2))*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p3,
     &    p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2 + 32._dp/(s(p1,
     &    p2))/(s(p1,p2))*s(p1,p4)*s(p2,p3) + 32._dp/(s(p1,p2))/(s(p1,p2
     &    ))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))*s(p2,p3)*
     &    s(p3,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(
     &    p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**
     &    2*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3)**2 + 32._dp/(
     &    s(p1,p2))/(s(p1,p4))*s(p1,p3)*s(p3,p4) + 8._dp/(s(p1,p2))/(s(
     &    p1,p4))*s(p2,p3)**2 - 16._dp/(s(p1,p2))/(s(p1,p4))*s(p2,p3)*s(
     &    p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p4))*s(p3,p4)**2 - 48._dp/(s(p1
     &    ,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p3,
     &    p4))*s(p1,p4)*s(p2,p3) - 16._dp/(s(p1,p4))*s(p1,p2) + 16._dp/(
     &    s(p1,p4))*s(p2,p3) + 16._dp/(s(p1,p4))*s(p3,p4) + 8._dp/(s(p1,
     &    p4))/(
     & s(p1,p4))*s(p1,p3)**2 + 16._dp/(s(p1,p4))/(s(p1,p4))*s(p1,p3)*s(
     &    p3,p4) + 8._dp/(s(p1,p4))/(s(p1,p4))*s(p2,p3)**2 + 8._dp/(s(p1,
     &    p4))/(s(p1,p4))*s(p3,p4)**2 - 16._dp/(s(p1,p4))/(s(p3,p4))*s(
     &    p1,p2)*s(p1,p3) - 48._dp/(s(p3,p4))*s(p1,p3) - 16._dp/(s(p3,p4)
     &    )*s(p2,p3) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp2*nDp3 * ( 16._dp - 16._dp
     &    /(s(p1,p2))*s(p1,p3) + 16._dp/(s(p1,p2))*s(p1,p4) - 16._dp/(s(
     &    p1,p2))*s(p2,p3) - 8._dp/(s(p1,p2))*s(p3,p4) - 32._dp/(s(p1,p2)
     &    )/(s(p1,p2))*s(p1,p3)*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)
     &    *s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4
     &    )**2*s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2 - 48._dp
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 16._dp/(s(p1,p2))/(
     &    s(p3,p4))*s(p1,p4)**2 - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*
     &    s(p2,p3) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p3)**2 - 16._dp/(s(
     &    p3,p4))*s(p1,p3) + 16._dp/(s(p3,p4))*s(p1,p4) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp2**2 * (  - 16._dp/(s(p1,
     &    p2))*s(p1,p3) - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)
     &     - 32._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2 - 16._dp/(s(p1,p2))
     &    /(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2 - 16._dp/(s(p1,p2)
     &    )/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p3) - 16._dp/(s(p1,p2
     &    ))/(s(p3,p4))*s(p1,p3)*s(p1,p4) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDp3**2 * (  - 16._dp/(s(p1,
     &    p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)**2 - 16._dp/(s(p1
     &    ,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2*s(p2,p3) - 16._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 32._dp/(s(p1,p2))/(s(p3
     &    ,p4))*s(p1,p4)**2 - 32._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p3) + 16._dp/(s(p3,p4))*s(p1,p2) + 48._dp/(s(p3,p4))*s(p1,p4)
     &     - 16._dp/(s(p3,p4))*s(p2,p3) )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDn * (  - 6._dp*s(p1,p2) +
     &    32._dp*s(p1,p3) + 16._dp*s(p2,p3) + 20._dp*s(p3,p4) - 12._dp/(s(
     &    p1,p2))*s(p1,p3)**2 + 32._dp/(s(p1,p2))*s(p1,p3)*s(p1,p4) - 16
     &    ._dp/(s(p1,p2))*s(p1,p3)*s(p2,p3) - 26._dp/(s(p1,p2))*s(p1,p3)*
     &    s(p3,p4) + 16._dp/(s(p1,p2))*s(p1,p4)*s(p2,p3) + 8._dp/(s(p1,p2
     &    ))*s(p1,p4)*s(p3,p4) - 4._dp/(s(p1,p2))*s(p2,p3)**2 - 10._dp/(
     &    s(p1,p2))*s(p2,p3)*s(p3,p4) - 10._dp/(s(p1,p2))*s(p3,p4)**2 -
     &    8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2*s(p1,p4) + 16._dp/(s(p1
     &    ,p2))/(s(p1,p2))*s(p1,p3)*s(p1,p4)**2 - 8._dp/(s(p1,p2))/(s(p1
     &    ,p2))*s(p1,p3)*s(p1,p4)*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p2)
     &    )*s(p1,p3)*s(p1,p4)*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))*s(
     &    p1,p4)**2*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)**2*
     &    s(p3,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*
     &    s(p1,p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)
     &    *s(p1,p4)**2*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))
     &    *s(p1,p4)**2*s(p2,p3)**2 )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDn * (  - 2._dp/(s(p1,p2))
     &    /(s(p2,p3))*s(p1,p4)*s(p3,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4)
     &    )*s(p1,p3)**2*s(p1,p4) + 24._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)
     &    *s(p1,p4)*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p3)**2 - 2._dp/(s(p1,p4))*s(p1,p2)**2 + 12._dp/(s(p1,p4))*s(p1
     &    ,p2)*s(p1,p3) + 4._dp/(s(p1,p4))*s(p1,p2)*s(p2,p3) + 8._dp/(s(
     &    p1,p4))*s(p1,p2)*s(p3,p4) - 6._dp/(s(p1,p4))*s(p1,p3)**2 - 8._dp
     &    /(s(p1,p4))*s(p1,p3)*s(p2,p3) - 10._dp/(s(p1,p4))*s(p1,p3)*s(
     &    p3,p4) - 2._dp/(s(p1,p4))*s(p2,p3)**2 - 4._dp/(s(p1,p4))*s(p2,
     &    p3)*s(p3,p4) - 4._dp/(s(p1,p4))*s(p3,p4)**2 - 2._dp/(s(p1,p4))
     &    /(s(p3,p4))*s(p1,p2)**2*s(p1,p3) + 4._dp/(s(p1,p4))/(s(p3,p4))
     &    *s(p1,p2)*s(p1,p3)**2 + 4._dp/(s(p1,p4))/(s(p3,p4))*s(p1,p2)*
     &    s(p1,p3)*s(p2,p3) - 2._dp/(s(p2,p3))*s(p1,p2)**2 - 2._dp/(s(p2,
     &    p3))*s(p1,p2)*s(p1,p4) + 4._dp/(s(p2,p3))*s(p1,p2)*s(p3,p4) +
     &    4._dp/(s(p2,p3))*s(p1,p4)*s(p3,p4) - 2._dp/(s(p2,p3))*s(p3,p4)
     &    **2 )
      qqgghn_ba = qqgghn_ba + s124**(-1)*xn*nDn * (  - 2._dp/(s(p3,p4))*
     &    s(p1,p2)*s(p1,p3) - 2._dp/(s(p3,p4))*s(p1,p2)*s(p2,p3) + 12._dp
     &    /(s(p3,p4))*s(p1,p3)**2 + 16._dp/(s(p3,p4))*s(p1,p3)*s(p2,p3)
     &     + 4._dp/(s(p3,p4))*s(p2,p3)**2 )
      qqgghn_ba = qqgghn_ba + xn*nDp1*nDp2 * ( 24._dp/(s(p1,p2)) + 32._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p3) + 96._dp/(s(p1,p2))/(s(p1,p2))
     &    *s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(
     &    p2,p3) + 64._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p4) + 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3) + 16._dp/(s(p1,p2)
     &    )/(s(p1,p4))*s(p2,p3) + 8._dp/(s(p1,p2))/(s(p1,p4))*s(p2,p4)
     &     + 16._dp/(s(p1,p2))/(s(p1,p4))*s(p3,p4) + 32._dp/(s(p1,p2))/(
     &    s(p3,p4))*s(p1,p3) - 8._dp/(s(p1,p4))/(s(p2,p3))*s(p3,p4) - 8.D
     &    0/(s(p2,p3)) )
      qqgghn_ba = qqgghn_ba + xn*nDp1*nDp3 * ( 48._dp/(s(p1,p2)) - 32._dp
     &    /(s(p1,p2))/(s(p1,p2))*s(p1,p3) + 32._dp/(s(p1,p2))/(s(p1,p2))
     &    *s(p1,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(
     &    p2,p4) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2
     &    ,p3) + 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,
     &    p4) - 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3) - 16._dp/(s(p1,p2))
     &    /(s(p1,p4))*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p4))*s(p2,p4)
     &     - 8._dp/(s(p1,p2))/(s(p1,p4))*s(p3,p4) - 8._dp/(s(p1,p2))/(s(
     &    p1,p4))/(s(p3,p4))*s(p1,p3)**2 - 16._dp/(s(p1,p2))/(s(p1,p4))
     &    /(s(p3,p4))*s(p1,p3)*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p1,p4))/(
     &    s(p3,p4))*s(p2,p3)**2 - 16._dp/(s(p1,p2))/(s(p1,p4))/(s(p3,p4)
     &    )*s(p2,p3)*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p1,p4))/(s(p3,p4))*
     &    s(p2,p4)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3) - 32._dp/(
     &    s(p1,p2))/(s(p3,p4))*s(p1,p4) + 16._dp/(s(p1,p2))/(s(p3,p4))*
     &    s(p2,p3) )
      qqgghn_ba = qqgghn_ba + xn*nDp1*nDp3 * ( 24._dp/(s(p1,p2))/(s(p3,
     &    p4))*s(p2,p4) + 16._dp/(s(p1,p4)) + 8._dp/(s(p1,p4))/(s(p2,p3))
     &    *s(p2,p4) + 16._dp/(s(p1,p4))/(s(p3,p4))*s(p1,p3) - 8._dp/(s(p2
     &    ,p3))/(s(p3,p4))*s(p2,p4) - 16._dp/(s(p3,p4)) )
      qqgghn_ba = qqgghn_ba + xn*nDp1**2 * ( 56._dp/(s(p1,p2)) + 32._dp/(
     &    s(p1,p2))/(s(p1,p2))*s(p1,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))*
     &    s(p2,p3) - 64._dp/(s(p1,p2))/(s(p1,p2))*s(p2,p4) - 32._dp/(s(p1
     &    ,p2))/(s(p1,p2))*s(p3,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3
     &    ,p4))*s(p1,p3)*s(p1,p4) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,
     &    p4))*s(p1,p3)*s(p2,p4) + 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4
     &    ))*s(p1,p4)*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))
     &    *s(p2,p3)*s(p2,p4) - 32._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*
     &    s(p2,p4)**2 + 16._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3) - 16._dp/(
     &    s(p1,p2))/(s(p1,p4))*s(p2,p3) - 24._dp/(s(p1,p2))/(s(p1,p4))*
     &    s(p2,p4) - 16._dp/(s(p1,p2))/(s(p1,p4))*s(p3,p4) - 16._dp/(s(p1
     &    ,p2))/(s(p1,p4))/(s(p3,p4))*s(p1,p3)*s(p2,p4) + 32._dp/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p3) + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p2,
     &    p3) + 16._dp/(s(p1,p4)) + 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4)
     &     + 16._dp/(s(p1,p4))/(s(p3,p4))*s(p1,p3) )
      qqgghn_ba = qqgghn_ba + xn*nDp2*nDp3 * (  - 16._dp/(s(p1,p2)) + 32
     &    ._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 32._dp
     &    /(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p4)**2 + 16._dp/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p3) - 24._dp/(s(p1,p2))/(s(p3,p4))*s(p1,
     &    p4) - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p3) - 8._dp/(s(p1,p2))
     &    /(s(p3,p4))*s(p2,p4) + 8._dp/(s(p2,p3)) + 8._dp/(s(p2,p3))/(s(
     &    p3,p4))*s(p1,p4) )
      qqgghn_ba = qqgghn_ba + xn*nDp2**2 * (  - 32._dp/(s(p1,p2))/(s(p1,
     &    p2))/(s(p3,p4))*s(p1,p4)**2 )
      qqgghn_ba = qqgghn_ba + xn*nDp3**2 * ( 8._dp/(s(p1,p2)) + 32._dp/(
     &    s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) + 16._dp/(s(
     &    p1,p2))/(s(p3,p4))*s(p1,p3) - 32._dp/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p4) + 16._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p3) + 16._dp/(s(p1,
     &    p2))/(s(p3,p4))*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,
     &    p4))*s(p1,p3)**2 + 16._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(
     &    p1,p3)*s(p1,p4) + 8._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p1,
     &    p4)**2 + 8._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p2,p3)**2 +
     &    16._dp/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p2,p3)*s(p2,p4) + 8.D
     &    0/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p2,p4)**2 - 8._dp/(s(p2,
     &    p3))/(s(p3,p4))*s(p2,p4) - 16._dp/(s(p3,p4)) )
      qqgghn_ba = qqgghn_ba + xn*nDn * ( 8._dp - 32._dp/(s(p1,p2))*s(p1,
     &    p3) - 22._dp/(s(p1,p2))*s(p1,p4) - 20._dp/(s(p1,p2))*s(p2,p3)
     &     - 22._dp/(s(p1,p2))*s(p2,p4) - 24._dp/(s(p1,p2))*s(p3,p4) + 4.D
     &    0/(s(p1,p2))/(s(p1,p2))*s(p1,p3)**2 - 32._dp/(s(p1,p2))/(s(p1,
     &    p2))*s(p1,p3)*s(p1,p4) + 4._dp/(s(p1,p2))/(s(p1,p2))*s(p1,p4)
     &    **2 - 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p1,
     &    p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2*s(p2,
     &    p4) - 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4)
     &    **2 - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4
     &    )*s(p2,p3) - 16._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*
     &    s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))*s(
     &    p1,p4)**2*s(p2,p3) + 4._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))/(
     &    s(p3,p4))*s(p1,p3)**2*s(p2,p4)**2 - 8._dp/(s(p1,p2))/(s(p1,p2)
     &    )/(s(p3,p4))/(s(p3,p4))*s(p1,p3)*s(p1,p4)*s(p2,p3)*s(p2,p4)
     &     + 4._dp/(s(p1,p2))/(s(p1,p2))/(s(p3,p4))/(s(p3,p4))*s(p1,p4)
     &    **2*s(p2,p3)**2 )
      qqgghn_ba = qqgghn_ba + xn*nDn * ( 4._dp/(s(p1,p2))/(s(p1,p4))*s(
     &    p1,p3)**2 + 4._dp/(s(p1,p2))/(s(p1,p4))*s(p1,p3)*s(p2,p3) + 4.D
     &    0/(s(p1,p2))/(s(p1,p4))*s(p1,p3)*s(p2,p4) + 2._dp/(s(p1,p2))/(
     &    s(p1,p4))*s(p1,p3)*s(p3,p4) + 4._dp/(s(p1,p2))/(s(p1,p4))/(s(
     &    p3,p4))*s(p1,p3)**2*s(p2,p4) + 4._dp/(s(p1,p2))/(s(p1,p4))/(s(
     &    p3,p4))*s(p1,p3)*s(p2,p3)*s(p2,p4) + 2._dp/(s(p1,p2))/(s(p1,p4
     &    ))/(s(p3,p4))*s(p1,p3)*s(p2,p4)**2 - 2._dp/(s(p1,p2))/(s(p2,p3
     &    ))*s(p1,p4)*s(p2,p4) - 2._dp/(s(p1,p2))/(s(p2,p3))*s(p2,p4)**2
     &     + 2._dp/(s(p1,p2))/(s(p2,p3))*s(p2,p4)*s(p3,p4) + 2._dp/(s(p1,
     &    p2))/(s(p2,p3))*s(p3,p4)**2 - 2._dp/(s(p1,p2))/(s(p2,p3))/(s(
     &    p3,p4))*s(p1,p4)**2*s(p2,p4) - 4._dp/(s(p1,p2))/(s(p2,p3))/(s(
     &    p3,p4))*s(p1,p4)*s(p2,p4)**2 - 2._dp/(s(p1,p2))/(s(p2,p3))/(s(
     &    p3,p4))*s(p2,p4)**3 - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)**2
     &     - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p3)*s(p1,p4) - 12._dp/(s(p1
     &    ,p2))/(s(p3,p4))*s(p1,p3)*s(p2,p3) - 6._dp/(s(p1,p2))/(s(p3,p4
     &    ))*s(p1,p3)*s(p2,p4) )
      qqgghn_ba = qqgghn_ba + xn*nDn * (  - 8._dp/(s(p1,p2))/(s(p3,p4))*
     &    s(p1,p4)**2 - 10._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p3)
     &     - 16._dp/(s(p1,p2))/(s(p3,p4))*s(p1,p4)*s(p2,p4) - 4._dp/(s(p1
     &    ,p2))/(s(p3,p4))*s(p2,p3)**2 - 10._dp/(s(p1,p2))/(s(p3,p4))*s(
     &    p2,p3)*s(p2,p4) - 8._dp/(s(p1,p2))/(s(p3,p4))*s(p2,p4)**2 + 2.D
     &    0/(s(p1,p4))*s(p1,p2) - 12._dp/(s(p1,p4))*s(p1,p3) - 4._dp/(s(
     &    p1,p4))*s(p2,p3) - 2._dp/(s(p1,p4))*s(p2,p4) - 4._dp/(s(p1,p4))
     &    *s(p3,p4) - 2._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4)**2 - 2._dp/(s(
     &    p1,p4))/(s(p2,p3))*s(p2,p4)*s(p3,p4) + 2._dp/(s(p1,p4))/(s(p3,
     &    p4))*s(p1,p2)*s(p1,p3) - 4._dp/(s(p1,p4))/(s(p3,p4))*s(p1,p3)
     &    **2 - 4._dp/(s(p1,p4))/(s(p3,p4))*s(p1,p3)*s(p2,p3) - 2._dp/(s(
     &    p1,p4))/(s(p3,p4))*s(p1,p3)*s(p2,p4) + 2._dp/(s(p2,p3))*s(p1,
     &    p2) - 8._dp/(s(p2,p3))*s(p2,p4) - 4._dp/(s(p2,p3))*s(p3,p4) - 4
     &    ._dp/(s(p2,p3))/(s(p3,p4))*s(p1,p4)*s(p2,p4) - 4._dp/(s(p2,p3))
     &    /(s(p3,p4))*s(p2,p4)**2 + 2._dp/(s(p3,p4))*s(p2,p3) )

      qqgghn_sym=  + s123**(-2)*xn**(-1)*nDn * (  - 4._dp*s(p1,p4)**2 -
     &    8._dp*s(p1,p4)*s(p2,p4) - 8._dp*s(p1,p4)*s(p3,p4) - 4._dp*s(p2,
     &    p4)**2 - 8._dp*s(p2,p4)*s(p3,p4) - 4._dp*s(p3,p4)**2 - 2._dp/(s(
     &    p1,p3))*s(p1,p2)*s(p1,p4)**2 - 4._dp/(s(p1,p3))*s(p1,p2)*s(p1,
     &    p4)*s(p2,p4) - 4._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4)*s(p3,p4) - 2
     &    ._dp/(s(p1,p3))*s(p1,p2)*s(p2,p4)**2 - 4._dp/(s(p1,p3))*s(p1,p2)
     &    *s(p2,p4)*s(p3,p4) - 2._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4)**2 - 2
     &    ._dp/(s(p2,p3))*s(p1,p2)*s(p1,p4)**2 - 4._dp/(s(p2,p3))*s(p1,p2)
     &    *s(p1,p4)*s(p2,p4) - 4._dp/(s(p2,p3))*s(p1,p2)*s(p1,p4)*s(p3,
     &    p4) - 2._dp/(s(p2,p3))*s(p1,p2)*s(p2,p4)**2 - 4._dp/(s(p2,p3))*
     &    s(p1,p2)*s(p2,p4)*s(p3,p4) - 2._dp/(s(p2,p3))*s(p1,p2)*s(p3,p4
     &    )**2 )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDp1*
     & nDp2 * (  - 64._dp*s(p3,p4) - 16._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4)
     &     - 16._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p3))*s(p3
     &    ,p4)**2 - 16._dp/(s(p1,p3))/(s(p1,p4))*s(p1,p2)*s(p3,p4)**2 -
     &    16._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4) - 16._dp/(s(p1,p4))*s(p1,p3
     &    )*s(p3,p4) - 16._dp/(s(p1,p4))*s(p3,p4)**2 + 16._dp/(s(p2,p3))*
     &    s(p1,p4)*s(p3,p4) - 16._dp/(s(p2,p3))*s(p3,p4)**2 - 16._dp/(s(
     &    p2,p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4)**2 + 16._dp/(s(p2,p4))*s(
     &    p1,p3)*s(p3,p4) - 16._dp/(s(p2,p4))*s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDp1*
     & nDp3 * (  - 16._dp*s(p1,p2) + 16._dp/(s(p1,p3))*s(p1,p2)**2 + 32._dp
     &    /(s(p1,p3))*s(p1,p2)*s(p1,p4) - 32._dp/(s(p1,p3))*s(p1,p2)*s(
     &    p3,p4) + 16._dp/(s(p1,p3))*s(p1,p4)**2 - 16._dp/(s(p1,p3))*s(p1
     &    ,p4)*s(p3,p4) - 16._dp/(s(p1,p3))/(s(p1,p4))*s(p1,p2)**2*s(p3,
     &    p4) - 16._dp/(s(p1,p4))*s(p1,p2)**2 - 16._dp/(s(p1,p4))*s(p1,p2
     &    )*s(p1,p3) - 16._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4) + 16._dp/(s(p2
     &    ,p3))*s(p1,p2)*s(p1,p4) + 16._dp/(s(p2,p3))*s(p1,p4)**2 - 16._dp
     &    /(s(p2,p3))*s(p1,p4)*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDp1**2
     &  * (  - 32._dp*s(p3,p4) - 32._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4) - 16
     &    ._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p1,p3))/(s(p1,p4))
     &    *s(p1,p2)**2*s(p3,p4) - 32._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4) -
     &    16._dp/(s(p1,p4))*s(p1,p3)*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDp2*
     & nDp3 * (  - 16._dp*s(p1,p2) + 16._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4)
     &     + 16._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4) + 16._dp/(s(p1,p3))*s(p1
     &    ,p4)**2 + 16._dp/(s(p1,p3))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p2,p3
     &    ))*s(p1,p2)*s(p3,p4) + 16._dp/(s(p2,p3))*s(p1,p4)**2 + 16._dp/(
     &    s(p2,p3))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p2,p3))/(s(p2,p4))*s(
     &    p1,p2)**2*s(p3,p4) + 16._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3) - 16.D
     &    0/(s(p2,p4))*s(p1,p2)*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDp2**2
     &  * (  - 32._dp*s(p3,p4) - 16._dp/(s(p2,p3))*s(p1,p2)*s(p3,p4) + 16
     &    ._dp/(s(p2,p3))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p2,p3))/(s(p2,p4))
     &    *s(p1,p2)**2*s(p3,p4) - 16._dp/(s(p2,p4))*s(p1,p2)*s(p3,p4) +
     &    16._dp/(s(p2,p4))*s(p1,p3)*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDp3**2
     &  * ( 16._dp*s(p1,p2) + 16._dp/(s(p1,p3))*s(p1,p2)**2 + 32._dp/(s(p1
     &    ,p3))*s(p1,p2)*s(p1,p4) + 16._dp/(s(p1,p3))*s(p1,p4)**2 + 16._dp
     &    /(s(p2,p3))*s(p1,p4)**2 )
      qqgghn_sym = qqgghn_sym + s123**(-1)*s124**(-1)*xn**(-1)*nDn * (
     &     - 8._dp*s(p1,p2)**2 + 16._dp*s(p1,p2)*s(p3,p4) - 8._dp*s(p3,p4)
     &    **2 - 2._dp/(s(p1,p3))*s(p1,p2)**3 + 4._dp/(s(p1,p3))*s(p1,p2)
     &    **2*s(p3,p4) - 2._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4)**2 - 2._dp/(
     &    s(p1,p3))/(s(p1,p4))*s(p1,p2)**3*s(p3,p4) + 4._dp/(s(p1,p3))/(
     &    s(p1,p4))*s(p1,p2)**2*s(p3,p4)**2 - 2._dp/(s(p1,p3))/(s(p1,p4)
     &    )*s(p1,p2)*s(p3,p4)**3 - 2._dp/(s(p1,p4))*s(p1,p2)**3 + 4._dp/(
     &    s(p1,p4))*s(p1,p2)**2*s(p3,p4) - 2._dp/(s(p1,p4))*s(p1,p2)*s(
     &    p3,p4)**2 - 2._dp/(s(p2,p3))*s(p1,p2)**3 + 4._dp/(s(p2,p3))*s(
     &    p1,p2)**2*s(p3,p4) - 2._dp/(s(p2,p3))*s(p1,p2)*s(p3,p4)**2 - 2
     &    ._dp/(s(p2,p3))/(s(p2,p4))*s(p1,p2)**3*s(p3,p4) + 4._dp/(s(p2,p3
     &    ))/(s(p2,p4))*s(p1,p2)**2*s(p3,p4)**2 - 2._dp/(s(p2,p3))/(s(p2
     &    ,p4))*s(p1,p2)*s(p3,p4)**3 - 2._dp/(s(p2,p4))*s(p1,p2)**3 + 4.D
     &    0/(s(p2,p4))*s(p1,p2)**2*s(p3,p4) - 2._dp/(s(p2,p4))*s(p1,p2)*
     &    s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDp1*nDp2 * (  - 32
     &    ._dp - 16._dp/(s(p1,p3))*s(p1,p2) + 8._dp/(s(p1,p3))*s(p3,p4) - 8
     &    ._dp/(s(p1,p3))/(s(p1,p4))*s(p1,p2)*s(p3,p4) + 8._dp/(s(p1,p3))
     &    /(s(p1,p4))*s(p3,p4)**2 - 32._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4
     &    )*s(p2,p4) - 16._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4)*s(p3,p4) -
     &    16._dp/(s(p1,p3))/(s(p2,p3))*s(p2,p4)*s(p3,p4) + 8._dp/(s(p1,p3
     &    ))/(s(p2,p4))*s(p1,p4)*s(p3,p4) + 8._dp/(s(p1,p3))/(s(p2,p4))*
     &    s(p3,p4)**2 + 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4)*s(p3,p4) +
     &    8._dp/(s(p1,p4))/(s(p2,p3))*s(p3,p4)**2 - 16._dp/(s(p2,p3))*s(
     &    p1,p2) - 8._dp/(s(p2,p3))*s(p3,p4) - 8._dp/(s(p2,p3))/(s(p2,p4)
     &    )*s(p1,p2)*s(p3,p4) + 8._dp/(s(p2,p3))/(s(p2,p4))*s(p3,p4)**2
     &     )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDp1*nDp3 * (  - 16
     &    ._dp/(s(p1,p3))*s(p1,p2) - 8._dp/(s(p1,p3))*s(p1,p4) + 8._dp/(s(
     &    p1,p3))*s(p2,p4) + 8._dp/(s(p1,p3))*s(p3,p4) + 8._dp/(s(p1,p3))
     &    /(s(p1,p4))*s(p1,p2)*s(p2,p4) + 16._dp/(s(p1,p3))/(s(p1,p4))*
     &    s(p1,p2)*s(p3,p4) - 8._dp/(s(p1,p3))/(s(p1,p4))*s(p2,p4)*s(p3,
     &    p4) - 16._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4)*s(p2,p4) - 16._dp/(
     &    s(p1,p3))/(s(p2,p3))*s(p1,p4)*s(p3,p4) + 16._dp/(s(p1,p3))/(s(
     &    p2,p3))*s(p2,p4)**2 + 16._dp/(s(p1,p4))*s(p1,p2) + 16._dp/(s(p1
     &    ,p4))*s(p1,p3) + 16._dp/(s(p1,p4))*s(p3,p4) - 8._dp/(s(p1,p4))
     &    /(s(p2,p3))*s(p2,p4)**2 - 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4)
     &    *s(p3,p4) - 8._dp/(s(p2,p3))*s(p1,p2) - 16._dp/(s(p2,p3))*s(p1,
     &    p4) + 8._dp/(s(p2,p3))*s(p2,p4) + 8._dp/(s(p2,p3))*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDp1**2 * ( 8._dp/(
     &    s(p1,p3))*s(p3,p4) + 8._dp/(s(p1,p3))/(s(p1,p4))*s(p1,p2)*s(p2
     &    ,p4) + 16._dp/(s(p1,p3))/(s(p1,p4))*s(p1,p2)*s(p3,p4) - 8._dp/(
     &    s(p1,p3))/(s(p1,p4))*s(p2,p4)*s(p3,p4) + 16._dp/(s(p1,p3))/(s(
     &    p2,p3))*s(p2,p4)**2 + 16._dp/(s(p1,p3))/(s(p2,p3))*s(p2,p4)*s(
     &    p3,p4) + 8._dp/(s(p1,p3))/(s(p2,p3))*s(p3,p4)**2 + 16._dp/(s(p1
     &    ,p4))*s(p2,p4) + 32._dp/(s(p1,p4))*s(p3,p4) + 8._dp/(s(p1,p4))
     &    /(s(p2,p3))*s(p1,p2)*s(p2,p4) + 8._dp/(s(p1,p4))/(s(p2,p3))*s(
     &    p1,p2)*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDp2*nDp3 * (  - 8.D
     &    0/(s(p1,p3))*s(p1,p2) - 8._dp/(s(p1,p3))*s(p1,p4) - 8._dp/(s(p1
     &    ,p3))*s(p3,p4) + 16._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4)**2 - 16
     &    ._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4)*s(p2,p4) - 16._dp/(s(p1,p3))
     &    /(s(p2,p3))*s(p2,p4)*s(p3,p4) - 8._dp/(s(p1,p3))/(s(p2,p4))*s(
     &    p1,p4)**2 - 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p4)*s(p3,p4) - 8.D
     &    0/(s(p2,p3))*s(p1,p4) + 8._dp/(s(p2,p3))*s(p2,p4) - 8._dp/(s(p2
     &    ,p3))*s(p3,p4) + 8._dp/(s(p2,p3))/(s(p2,p4))*s(p1,p2)*s(p1,p4)
     &     + 16._dp/(s(p2,p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4) - 8._dp/(s(p2
     &    ,p3))/(s(p2,p4))*s(p1,p4)*s(p3,p4) - 16._dp/(s(p2,p4))*s(p1,p3
     &    ) + 16._dp/(s(p2,p4))*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDp2**2 * ( 16._dp/(
     &    s(p1,p3))/(s(p2,p3))*s(p1,p4)**2 + 16._dp/(s(p1,p3))/(s(p2,p3)
     &    )*s(p1,p4)*s(p3,p4) + 8._dp/(s(p1,p3))/(s(p2,p3))*s(p3,p4)**2
     &     + 8._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p1,p4) + 8._dp/(s(p1,
     &    p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4) - 8._dp/(s(p2,p3))*s(p3,p4)
     &     + 8._dp/(s(p2,p3))/(s(p2,p4))*s(p1,p2)*s(p1,p4) + 16._dp/(s(p2
     &    ,p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4) - 8._dp/(s(p2,p3))/(s(p2,p4
     &    ))*s(p1,p4)*s(p3,p4) + 16._dp/(s(p2,p4))*s(p1,p4) + 32._dp/(s(
     &    p2,p4))*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDp3**2 * (  - 16._dp
     &     - 16._dp/(s(p1,p3))*s(p1,p2) - 8._dp/(s(p1,p3))*s(p1,p4) + 8._dp
     &    /(s(p1,p3))*s(p2,p4) + 8._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4)**2
     &     + 8._dp/(s(p1,p3))/(s(p2,p3))*s(p2,p4)**2 - 8._dp/(s(p2,p3))*
     &    s(p1,p4) + 8._dp/(s(p2,p3))*s(p2,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDn * ( 4._dp*s(p1,
     &    p2) - 6._dp*s(p1,p4) - 6._dp*s(p2,p4) - 8._dp*s(p3,p4) + 2._dp/(
     &    s(p1,p3))*s(p1,p2)**2 - 2._dp/(s(p1,p3))*s(p1,p2)*s(p1,p4) - 4
     &    ._dp/(s(p1,p3))*s(p1,p2)*s(p3,p4) - 2._dp/(s(p1,p3))*s(p1,p4)**2
     &     - 4._dp/(s(p1,p3))*s(p1,p4)*s(p2,p4) - 4._dp/(s(p1,p3))*s(p1,
     &    p4)*s(p3,p4) - 2._dp/(s(p1,p3))*s(p2,p4)**2 - 4._dp/(s(p1,p3))*
     &    s(p2,p4)*s(p3,p4) - 2._dp/(s(p1,p3))*s(p3,p4)**2 + 2._dp/(s(p1,
     &    p3))/(s(p1,p4))*s(p1,p2)**2*s(p3,p4) - 2._dp/(s(p1,p3))/(s(p1,
     &    p4))*s(p1,p2)*s(p3,p4)**2 - 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1,
     &    p2)*s(p1,p4)**2 - 4._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p1,p4
     &    )*s(p3,p4) - 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p2)*s(p3,p4)**2
     &     + 2._dp/(s(p1,p4))*s(p1,p2)**2 - 2._dp/(s(p1,p4))*s(p1,p2)*s(
     &    p2,p4) - 4._dp/(s(p1,p4))*s(p1,p2)*s(p3,p4) - 2._dp/(s(p1,p4))*
     &    s(p2,p4)**2 + 2._dp/(s(p1,p4))*s(p3,p4)**2 - 2._dp/(s(p1,p4))/(
     &    s(p2,p3))*s(p1,p2)*s(p2,p4)**2 - 4._dp/(s(p1,p4))/(s(p2,p3))*
     &    s(p1,p2)*s(p2,p4)*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s123**(-1)*xn**(-1)*nDn * (  - 2._dp/(s(
     &    p1,p4))/(s(p2,p3))*s(p1,p2)*s(p3,p4)**2 + 2._dp/(s(p2,p3))*s(
     &    p1,p2)**2 - 2._dp/(s(p2,p3))*s(p1,p2)*s(p2,p4) - 4._dp/(s(p2,p3
     &    ))*s(p1,p2)*s(p3,p4) - 2._dp/(s(p2,p3))*s(p1,p4)**2 - 4._dp/(s(
     &    p2,p3))*s(p1,p4)*s(p2,p4) - 4._dp/(s(p2,p3))*s(p1,p4)*s(p3,p4)
     &     - 2._dp/(s(p2,p3))*s(p2,p4)**2 - 4._dp/(s(p2,p3))*s(p2,p4)*s(
     &    p3,p4) - 2._dp/(s(p2,p3))*s(p3,p4)**2 + 2._dp/(s(p2,p3))/(s(p2,
     &    p4))*s(p1,p2)**2*s(p3,p4) - 2._dp/(s(p2,p3))/(s(p2,p4))*s(p1,
     &    p2)*s(p3,p4)**2 + 2._dp/(s(p2,p4))*s(p1,p2)**2 - 2._dp/(s(p2,p4
     &    ))*s(p1,p2)*s(p1,p4) - 4._dp/(s(p2,p4))*s(p1,p2)*s(p3,p4) - 2.D
     &    0/(s(p2,p4))*s(p1,p4)**2 + 2._dp/(s(p2,p4))*s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + s124**(-2)*xn**(-1)*nDn * (  - 4._dp*s(
     &    p1,p3)**2 - 8._dp*s(p1,p3)*s(p2,p3) - 8._dp*s(p1,p3)*s(p3,p4)
     &     - 4._dp*s(p2,p3)**2 - 8._dp*s(p2,p3)*s(p3,p4) - 4._dp*s(p3,p4)
     &    **2 - 2._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3)**2 - 4._dp/(s(p1,p4))*
     &    s(p1,p2)*s(p1,p3)*s(p2,p3) - 4._dp/(s(p1,p4))*s(p1,p2)*s(p1,p3
     &    )*s(p3,p4) - 2._dp/(s(p1,p4))*s(p1,p2)*s(p2,p3)**2 - 4._dp/(s(
     &    p1,p4))*s(p1,p2)*s(p2,p3)*s(p3,p4) - 2._dp/(s(p1,p4))*s(p1,p2)
     &    *s(p3,p4)**2 - 2._dp/(s(p2,p4))*s(p1,p2)*s(p1,p3)**2 - 4._dp/(
     &    s(p2,p4))*s(p1,p2)*s(p1,p3)*s(p2,p3) - 4._dp/(s(p2,p4))*s(p1,
     &    p2)*s(p1,p3)*s(p3,p4) - 2._dp/(s(p2,p4))*s(p1,p2)*s(p2,p3)**2
     &     - 4._dp/(s(p2,p4))*s(p1,p2)*s(p2,p3)*s(p3,p4) - 2._dp/(s(p2,p4
     &    ))*s(p1,p2)*s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDp1*nDp2 * (  - 32
     &    ._dp + 16._dp/(s(p1,p4))/(s(p2,p4))*s(p1,p3)**2 + 16._dp/(s(p1,p4
     &    ))/(s(p2,p4))*s(p1,p3)*s(p3,p4) + 16._dp/(s(p1,p4))/(s(p2,p4))
     &    *s(p2,p3)**2 + 16._dp/(s(p1,p4))/(s(p2,p4))*s(p2,p3)*s(p3,p4)
     &     - 16._dp/(s(p2,p4))*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDp1*nDp3 * (  - 16
     &    ._dp + 16._dp/(s(p1,p4))*s(p1,p3) - 16._dp/(s(p1,p4))*s(p2,p3) +
     &    16._dp/(s(p1,p4))*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDp1**2 * (  - 16._dp
     &     - 16._dp/(s(p1,p4))*s(p1,p2) + 16._dp/(s(p1,p4))*s(p3,p4) - 8.D
     &    0/(s(p1,p4))/(s(p1,p4))*s(p1,p3)**2 - 16._dp/(s(p1,p4))/(s(p1,
     &    p4))*s(p1,p3)*s(p3,p4) - 8._dp/(s(p1,p4))/(s(p1,p4))*s(p2,p3)
     &    **2 - 8._dp/(s(p1,p4))/(s(p1,p4))*s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDp2*nDp3 * (  - 16
     &    ._dp - 16._dp/(s(p2,p4))*s(p1,p2) - 16._dp/(s(p2,p4))*s(p1,p3) +
     &    16._dp/(s(p2,p4))*s(p2,p3) + 16._dp/(s(p2,p4))*s(p3,p4) )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDp2**2 * (  - 16._dp
     &     - 16._dp/(s(p2,p4))*s(p1,p2) - 8._dp/(s(p2,p4))/(s(p2,p4))*s(
     &    p1,p3)**2 - 8._dp/(s(p2,p4))/(s(p2,p4))*s(p2,p3)**2 - 16._dp/(
     &    s(p2,p4))/(s(p2,p4))*s(p2,p3)*s(p3,p4) - 8._dp/(s(p2,p4))/(s(
     &    p2,p4))*s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDp3**2 * (  - 16._dp
     &     )
      qqgghn_sym = qqgghn_sym + s124**(-1)*xn**(-1)*nDn * ( 8._dp*s(p1,
     &    p2) - 8._dp*s(p3,p4) + 2._dp/(s(p1,p3))*s(p1,p2)**2 - 4._dp/(s(
     &    p1,p3))*s(p1,p2)*s(p3,p4) + 2._dp/(s(p1,p3))*s(p3,p4)**2 + 2._dp
     &    /(s(p1,p4))*s(p1,p3)**2 + 4._dp/(s(p1,p4))*s(p1,p3)*s(p2,p3)
     &     + 2._dp/(s(p1,p4))*s(p2,p3)**2 - 2._dp/(s(p1,p4))*s(p3,p4)**2
     &     + 2._dp/(s(p2,p3))*s(p1,p2)**2 - 4._dp/(s(p2,p3))*s(p1,p2)*s(
     &    p3,p4) + 2._dp/(s(p2,p3))*s(p3,p4)**2 + 2._dp/(s(p2,p4))*s(p1,
     &    p3)**2 + 4._dp/(s(p2,p4))*s(p1,p3)*s(p2,p3) + 2._dp/(s(p2,p4))*
     &    s(p2,p3)**2 - 2._dp/(s(p2,p4))*s(p3,p4)**2 )
      qqgghn_sym = qqgghn_sym + xn**(-1)*nDp1*nDp2 * ( 16._dp/(s(p1,p3))
     &     + 8._dp/(s(p1,p3))/(s(p2,p4))*s(p3,p4) + 8._dp/(s(p1,p4))/(s(
     &    p2,p3))*s(p3,p4) + 16._dp/(s(p2,p3)) )
      qqgghn_sym = qqgghn_sym + xn**(-1)*nDp1*nDp3 * ( 8._dp/(s(p1,p3))
     &     - 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4) )
      qqgghn_sym = qqgghn_sym + xn**(-1)*nDp1**2 * (  - 8._dp/(s(p1,p3))
     &    /(s(p1,p4))*s(p2,p4) - 8._dp/(s(p1,p3))/(s(p1,p4))*s(p3,p4) +
     &    16._dp/(s(p1,p4)) - 8._dp/(s(p1,p4))/(s(p2,p3))*s(p2,p4) )
      qqgghn_sym = qqgghn_sym + xn**(-1)*nDp2*nDp3 * (  - 8._dp/(s(p1,p3
     &    ))/(s(p2,p4))*s(p1,p4) + 8._dp/(s(p2,p3)) + 16._dp/(s(p2,p4)) )
      qqgghn_sym = qqgghn_sym + xn**(-1)*nDp2**2 * (  - 8._dp/(s(p1,p3))
     &    /(s(p2,p4))*s(p1,p4) - 8._dp/(s(p2,p3))/(s(p2,p4))*s(p1,p4) -
     &    8._dp/(s(p2,p3))/(s(p2,p4))*s(p3,p4) + 16._dp/(s(p2,p4)) )
      qqgghn_sym = qqgghn_sym + xn**(-1)*nDn * (  - 8._dp - 2._dp/(s(p1,
     &    p3))*s(p1,p2) + 2._dp/(s(p1,p3))*s(p1,p4) + 2._dp/(s(p1,p3))*s(
     &    p3,p4) + 4._dp/(s(p1,p3))/(s(p2,p3))*s(p1,p4)**2 + 8._dp/(s(p1,
     &    p3))/(s(p2,p3))*s(p1,p4)*s(p2,p4) + 8._dp/(s(p1,p3))/(s(p2,p3)
     &    )*s(p1,p4)*s(p3,p4) + 4._dp/(s(p1,p3))/(s(p2,p3))*s(p2,p4)**2
     &     + 8._dp/(s(p1,p3))/(s(p2,p3))*s(p2,p4)*s(p3,p4) + 4._dp/(s(p1,
     &    p3))/(s(p2,p3))*s(p3,p4)**2 + 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1
     &    ,p4)**2 + 2._dp/(s(p1,p3))/(s(p2,p4))*s(p1,p4)*s(p3,p4) + 2._dp
     &    /(s(p1,p4))/(s(p2,p3))*s(p2,p4)**2 + 2._dp/(s(p1,p4))/(s(p2,p3
     &    ))*s(p2,p4)*s(p3,p4) + 4._dp/(s(p1,p4))/(s(p2,p4))*s(p3,p4)**2
     &     - 2._dp/(s(p2,p3))*s(p1,p2) + 2._dp/(s(p2,p3))*s(p2,p4) + 2._dp
     &    /(s(p2,p3))*s(p3,p4) )

      return
      end
