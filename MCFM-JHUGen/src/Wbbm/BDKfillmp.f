      subroutine BDKfillmp(k1,k2,k3,k4,k5,k6,mq,za,zb,ASLmp)
************************************************************************
*     Author: R.K. Ellis                                               *
*     August 2010.                                                     *
************************************************************************
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      integer k1,k2,k3,k4,k5,k6,e
      double complex Lnrat,xl14,ASLmp(-2:0),C0(-1:0),
     & B0p23(-1:0),B0p3(-1:0),VV(-2:0)
      double complex BDK1211mp,atree,lo,qlI2,qlI3
      double precision mqsq,p2Dp3,mq

      lo=atree('sl',k3,k2,k1,k4,k5,k6,zb,za)

      xl14=Lnrat(musq,-s(k1,k4))
      VV(-2)=-cone
      VV(-1)=-(dcmplx(1.5d0)+xl14)
      VV(0)=-(dcmplx(4d0)+1.5d0*xl14+0.5d0*xl14**2)

      mqsq=mq**2
c      write(6,*) 'BDKFill:musq',musq
      do e=-1,0
      C0(e)=qlI3(mqsq,s(k2,k3),mqsq,zip,mqsq,mqsq,musq,e)
      B0p23(e)=qlI2(s(k2,k3),mqsq,mqsq,musq,e)
      B0p3(e)=qlI2(mqsq,zip,mqsq,musq,e)
      enddo
      p2Dp3=0.5d0*(s(k2,k3)-2d0*mqsq)
    
C---Contribution of Final vertex function (massless case)
c       ASLmp(-2)=ASLmp(-2)-lo
c      ASLmp(-1)=ASLmp(-1)-(1.5d0+xl23)*lo
c      ASL(0)=ASL(0)-(3.5d0+1.5d0*xl23+0.5d0*xl23**2)*lo

C---Contribution of Final vertex function (massive  case)
C   (-2*p2.p3*C0(p2,p3,m,0,m)-3/2*B0(p23,m,m)+2*B0(p3,0,m)-1/2)

      VV(-1)=VV(-1)
     & +(-2d0*p2Dp3*C0(-1)-1.5d0*B0p23(-1)+2d0*B0p3(-1))
      VV(0)=VV(0)
     & +(-2d0*p2Dp3*C0(0)-1.5d0*B0p23(0)+2d0*B0p3(0)-dcmplx(0.5d0))

      ASLmp(-2)=VV(-2)*lo
      ASLmp(-1)=VV(-1)*lo
      ASLmp(0)=VV(0)*lo+BDK1211mp(k3,k2,k1,k4,k5,k6,zb,za)
       
      end
 
 
 

