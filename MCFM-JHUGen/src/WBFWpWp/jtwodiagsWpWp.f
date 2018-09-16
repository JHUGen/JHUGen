      subroutine amptwodiagsWpWp(p1,p2,p3,p4,p5,p6,p7,p8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: amp,propW34,propW56,propw28,propw17,zba2
      real(dp):: t3,s34,s56,s17,s28,
     & s134,s234,s156,s256,s347,s348,s567,s568
C-----Begin statement functions
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
C-----end statement functions

      s17=s(p1,p7)
      s28=s(p2,p8)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s234=t3(p2,p3,p4)
      s134=t3(p1,p3,p4)
      s156=t3(p1,p5,p6)
      s256=t3(p2,p5,p6)
      s347=t3(p3,p4,p7)
      s348=t3(p3,p4,p8)
      s567=t3(p5,p6,p7)
      s568=t3(p5,p6,p8)
      propw28=s28-cwmass2
      propw17=s17-cwmass2
      propW34=s34-cwmass2
      propW56=s56-cwmass2
      amp=za(p8,p3)*zb(p2,p6)*zba2(p4,p3,p8,p7)*zba2(p1,p2,p6,p5)*
     & propw56**(-1)*propw34**(-1)*s256**(-1)*s348**(-1)*propw17**(-1)
     &  + za(p8,p5)*zb(p2,p4)*zba2(p1,p2,p4,p3)*zba2(p6,p5,p8,p7)*
     & propw56**(-1)*propw34**(-1)*propw17**(-1)*s234**(-1)*s568**(-1)
     &  - za(p3,p7)*zb(p1,p6)*zba2(p2,p1,p6,p5)*zba2(p4,p3,p7,p8)*
     & propw56**(-1)*propw34**(-1)*propw28**(-1)*s156**(-1)*s347**(-1)
     &  - za(p7,p5)*zb(p4,p1)*zba2(p2,p1,p4,p3)*zba2(p6,p5,p7,p8)*
     & propw56**(-1)*propw34**(-1)*s134**(-1)*s567**(-1)*propw28**(-1)
      amp=amp/cxw**3
      return
      end
