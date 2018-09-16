      subroutine ampmidWpWp(p1,p2,p3,p4,p5,p6,p7,p8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cmplxmass.f'
      include 'masses.f'
      include 'runstring.f'
      include 'WWbits.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: sqwmass,amp,
     & propw34,propw56,propw28,propw17,proph1347,proph1567
      real(dp):: t4,s34,s56,s17,s28,s1347,s1567
C-----Begin statement functions
      t4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
      include 'cplx.h'
C-----end statement functions

c--- special fix for Madgraph check
      if (index(runstring,'mad') > 0) then
        sqwmass=cplx2(wmass**2,zip)
      else
        sqwmass=cwmass2
      endif

      s1347=t4(p1,p3,p4,p7)
      s1567=t4(p1,p5,p6,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s17=s(p1,p7)
      s28=s(p2,p8)

      proph1347=cplx2(s1347-hmass**2,hmass*hwidth)
      proph1567=cplx2(s1567-hmass**2,hmass*hwidth)
      propw34=s34-cwmass2
      propw56=s56-cwmass2
      propw17=s17-cwmass2
      propw28=s28-cwmass2

      amp= + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1) * ( 2._dp*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)
     &     - za(p3,p7)*za(p5,p8)*zb(p1,p4)*zb(p2,p6) - za(p3,p8)*za(p5,
     &    p7)*zb(p1,p6)*zb(p2,p4) )
      amp = amp + Hbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1) * (  - za(p3,p7)*za(p5,p8)*zb(p1,p4)*zb(p2,p6)*
     &    sqwmass*proph1347**(-1) - za(p3,p8)*za(p5,p7)*zb(p1,p6)*zb(p2
     &    ,p4)*sqwmass*proph1567**(-1) )
      amp=amp/cxw**3
      return
      end
