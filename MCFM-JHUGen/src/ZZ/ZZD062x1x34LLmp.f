      subroutine ZZD062x1x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, October 2013
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      integer:: k1,k2,k3,k4,k5,k6
      integer:: h3,h5,j1,j2,j3,j4,j5,j6
      real(dp):: mt,s134
      complex(dp):: zab2,amp,Xmp(2,2),Xpm(2,2)

C---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
C---end statement functions

      s134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      k1=j1
      k2=j2

      do h3=1,2
      do h5=1,2
      if (h3 == 1) then
        k3=j3
        k4=j4
      elseif (h3 == 2) then
        k3=j4
        k4=j3
      endif
      if (h5 == 1) then
        k5=j5
        k6=j6
      elseif (h5 == 2) then
        k5=j6
        k6=j5
      endif

      amp=
     & -zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)**3/zb(k3,k4)/za(k5,k6)
     & *((zab2(k2,k1,k3,k4)**2*zab2(k5,k3,k4,k1)**2)/s134
     &   +s134*za(k2,k5)**2*zb(k1,k4)**2)

      Xmp(h3,h5)=amp

      if (docheck) call ggZZcapture('d62x1x34',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              amp,czip,czip)

      enddo
      enddo

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5)=conjg(Xmp(3-h3,3-h5))
      enddo
      enddo

      return
      end
