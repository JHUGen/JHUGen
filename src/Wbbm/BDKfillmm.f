      subroutine BDKfillmm(k1,k2,k3,k4,k5,k6,mq,za,zb,ASLmm,Aslpp)
************************************************************************
*     Author: R.K. Ellis                                               *
*     August, 2010.                                                    *
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      integer k1,k2,k3,k4,k5,k6,e
      double complex Lnrat,xl14,ASLmm(-2:0),ASLpp(-2:0),
     & A(-1:0),C0(-1:0),B0p23(-1:0),B0p3(-1:0),VV(-2:0),zba2
      double complex BDK1211mm,loreduced,lomm,lopp,qlI1,qlI2,qlI3,
     & BDKfinite
      double precision mq,mqsq,p2Dp3,s123,s234
      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)

      s123=+s(k1,k2)+s(k2,k3)+s(k3,k1)
      s234=+s(k2,k3)+s(k3,k4)+s(k4,k2)
      loreduced=2d0/(s(k5,k6)*s(k2,k3))
     & *(+za(k1,k6)/s234*zba2(k5,k3,k4,k2)*zb(k2,k4)
     &   -zb(k5,k4)/s123*za(k1,k2)*zba2(k2,k1,k3,k6))
     
      lomm=mq/zb(k2,k3)*loreduced
      lopp=mq/za(k2,k3)*loreduced
      xl14=Lnrat(musq,-s(k1,k4))

      VV(-2)=-cone
      VV(-1)=-(dcmplx(1.5d0)+xl14)
      VV(0)=-(dcmplx(4d0)+1.5d0*xl14+0.5d0*xl14**2)


      mqsq=mq**2
      do e=-1,0
      C0(e)=qlI3(mqsq,s(k2,k3),mqsq,zip,mqsq,mqsq,musq,e)
      B0p23(e)=qlI2(s(k2,k3),mqsq,mqsq,musq,e)
      B0p3(e)=qlI2(mqsq,zip,mqsq,musq,e)
      A(e)=qlI1(mqsq,musq,e)
      enddo
      p2Dp3=0.5d0*(s(k2,k3)-2d0*mqsq)
    
C---Contribution of Final vertex function (massive  case)
C   *(-2*p2.p3*C0(p2,p3,m,0,m)-B0(p23,m,m)+B0(p3,0,m)+1/2*A0(m)/m^2);

      VV(-1)=VV(-1)
     & +(-2d0*p2Dp3*C0(-1)+dcmplx(half))
      VV(0)=VV(0)
     & +(-2d0*p2Dp3*C0(0)-B0p23(0)+B0p3(0)+0.5d0*A(0)/mqsq)

      BDKfinite=BDK1211mm(k1,k2,k3,k4,k5,k6,za,zb)

      ASLmm(-2)=VV(-2)*lomm
      ASLmm(-1)=VV(-1)*lomm
      ASLmm(0)=VV(0)*lomm+BDKfinite*mq/zb(k2,k3)

      ASLpp(-2)=VV(-2)*lopp
      ASLpp(-1)=VV(-1)*lopp
      ASLpp(0)=VV(0)*lopp+BDKfinite*mq/za(k2,k3)
      return
      end
 
 
 
 

