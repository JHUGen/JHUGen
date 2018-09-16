      subroutine ampBdiags(p1,p2,p3,p4,p5,p6,p7,p8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zba2,amp,
     & propw34,propw56,propw28,propw17,rxw
      real(dp):: t3,s17,s28,s34,s56,
     & s137,s147,s157,s167,s238,s248,s258,s268
C-----Begin statement functions
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      t3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s17=s(p1,p7)
      s28=s(p2,p8)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s137=t3(p1,p3,p7)
      s147=t3(p1,p4,p7)
      s157=t3(p1,p5,p7)
      s167=t3(p1,p6,p7)
      s238=t3(p2,p3,p8)
      s248=t3(p2,p4,p8)
      s258=t3(p2,p5,p8)
      s268=t3(p2,p6,p8)
      propw34=s34-cwmass2
      propw56=s56-cwmass2
      propw17=s17-cwmass2
      propw28=s28-cwmass2

      amp= + propw56**(-1)*propw28**(-1)*propw17**(-1) * ( za(p7,p3)*
     &    zb(p2,p4)*zba2(p1,p3,p7,p5)*zba2(p6,p2,p4,p8)*s248**(-1)*
     &    s137**(-1) + za(p8,p3)*zb(p1,p4)*zba2(p2,p3,p8,p5)*zba2(p6,p1
     &    ,p4,p7)*s238**(-1)*s147**(-1) )
      amp = amp + propw34**(-1)*propw28**(-1)*propw17**(-1) * ( za(p7,
     &    p5)*zb(p2,p6)*zba2(p1,p5,p7,p3)*zba2(p4,p2,p6,p8)*s268**(-1)*
     &    s157**(-1) + za(p8,p5)*zb(p1,p6)*zba2(p2,p5,p8,p3)*zba2(p4,p1
     &    ,p6,p7)*s258**(-1)*s167**(-1) )
      amp=amp/cxw**3
      return
      end
