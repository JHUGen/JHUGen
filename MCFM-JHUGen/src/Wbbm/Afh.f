      subroutine Afh(k1,k2,k3,k4,k5,k6,mq,Afmm,Afmp,Afpm,Afpp,
     & Ahmm,Ahmp,Ahpm,Ahpp,mhloopsq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'scale.f'
      integer:: k1,k2,k3,k4,k5,k6
      complex(dp):: Afmm(-2:0),Afmp(-2:0),Afpm(-2:0),Afpp(-2:0),
     & Ahmm(-2:0),Ahmp(-2:0),Ahpm(-2:0),Ahpp(-2:0),
     & LOmm,LOmp,LOpm,LOpp,B0f,B0fm,A0fm,qlI2,qlI1
      real(dp):: s23,mhloopsq,mq

      s23=s(k2,k3)


      call a6treemass(k1,k2,k3,k4,k5,k6,mq,LOmm,LOmp,LOpm,LOpp)

c--- compute finite parts of scalar integrals
      B0f=qlI2(s23,zip,zip,musq,0)
      B0fm=qlI2(s23,mhloopsq,mhloopsq,musq,0)
      A0fm=qlI1(mhloopsq,musq,0)

c--- fill coefficients according to:
c        tagheavy -> -(2/3*(1/3-1/e-B0f(p23,mb,mb))
c                   -4/3*mb^2*iza(k2,k3)*izb(k3,k2)*(B0f(p23,mb,mb)-A0f(mb)+1)
c                    )*ampLO;
c        taglight -> (-1/3+1/e+B0f(p23,0,0))*2/3*ampLO;

c--- mp
      Afmp(-2)=czip
      Afmp(-1)=2d0/3d0*LOmp
      Afmp( 0)=(-1d0/3d0+B0f)*2d0/3d0*LOmp

      Ahmp(-2)=czip
      Ahmp(-1)=2d0/3d0*LOmp
      Ahmp( 0)=-(2d0/3d0*(1d0/3d0-B0fm)
     &         -4d0/3d0*mhloopsq/s23*(B0fm-A0fm/mhloopsq+1d0))*LOmp

c--- pm
      Afpm(-2)=czip
      Afpm(-1)=2d0/3d0*LOpm
      Afpm( 0)=(-1d0/3d0+B0f)*2d0/3d0*LOpm

      Ahpm(-2)=czip
      Ahpm(-1)=2d0/3d0*LOpm
      Ahpm( 0)=-(2d0/3d0*(1d0/3d0-B0fm)
     &         -4d0/3d0*mhloopsq/s23*(B0fm-A0fm/mhloopsq+1d0))*LOpm

c--- mm
      Afmm(-2)=czip
      Afmm(-1)=2d0/3d0*LOmm
      Afmm( 0)=(-1d0/3d0+B0f)*2d0/3d0*LOmm

      Ahmm(-2)=czip
      Ahmm(-1)=2d0/3d0*LOmm
      Ahmm( 0)=-(2d0/3d0*(1d0/3d0-B0fm)
     &         -4d0/3d0*mhloopsq/s23*(B0fm-A0fm/mhloopsq+1d0))*LOmm

c--- pp
      Afpp(-2)=czip
      Afpp(-1)=2d0/3d0*LOpp
      Afpp( 0)=(-1d0/3d0+B0f)*2d0/3d0*LOpp

      Ahpp(-2)=czip
      Ahpp(-1)=2d0/3d0*LOpp
      Ahpp( 0)=-(2d0/3d0*(1d0/3d0-B0fm)
     &         -4d0/3d0*mhloopsq/s23*(B0fm-A0fm/mhloopsq+1d0))*LOpp

      return
      end

