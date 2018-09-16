      function gmgmjetn(zanb,j1,j2,p3,p4,p5)
      implicit none
      include 'types.f'
      real(dp):: gmgmjetn

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: i,j1,j2,p1,p2,p3,p4,p5
      complex(dp):: zanb(mxpart,mxpart),amp_pp,amp_mp
c--- amplitude squared for
c---    0 --> q(p1) + qb(p2) + gluon(p3) + photon(p4) + photon(p5)
c---    with gluon(p3) contracted with arbitrary vector n instead
c---    of polarization vector

      gmgmjetn=0._dp

c--- the loop accounts for the opposite helicities
      do i=1,2
      if (i == 1) then
        p1=j1
        p2=j2
      else
        p1=j2
        p2=j1
      endif

c--- amplitude for p4+, p5+
       amp_pp=
     &+ zanb(p1,p3)*za(p1,p2)*za(p1,p3)
     &   /za(p1,p4)/za(p1,p5)/za(p2,p4)/za(p2,p5)

c--- amplitude for p4-, p5+
       amp_mp=
     &- zanb(p1,p1)*za(p1,p4)*zb(p2,p5)
     &   /za(p1,p3)/za(p2,p5)/zb(p1,p3)/zb(p2,p4)
     &+ zanb(p1,p2)*za(p1,p2)*zb(p1,p2)
     &   /za(p1,p5)/za(p2,p5)/zb(p1,p4)/zb(p2,p4)
     &- zanb(p1,p3)*za(p3,p4)*zb(p2,p5)
     &   /za(p1,p3)/za(p2,p5)/zb(p1,p3)/zb(p2,p4)
     &+ zanb(p1,p5)*zb(p1,p2)
     &   /za(p2,p5)/zb(p1,p4)/zb(p2,p4)
     &- zanb(p2,p2)*za(p1,p4)*zb(p2,p5)
     &   /za(p1,p5)/za(p2,p3)/zb(p1,p4)/zb(p2,p3)
     &- zanb(p3,p2)*za(p1,p4)*zb(p3,p5)
     &   /za(p1,p5)/za(p2,p3)/zb(p1,p4)/zb(p2,p3)
     &- zanb(p4,p2)*za(p1,p2)/za(p1,p5)/za(p2,p5)/zb(p1,p4)
     &- zanb(p4,p5)/za(p2,p5)/zb(p1,p4)

      gmgmjetn=gmgmjetn+abs(amp_pp)**2+abs(amp_mp)**2
      enddo

      return
      end

