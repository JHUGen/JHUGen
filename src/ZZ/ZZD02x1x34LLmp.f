      subroutine ZZD02x1x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      implicit none
c--- Author: J. M. Campbell, October 2013
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      integer k1,k2,k3,k4,k5,k6
      integer h3,h5,j1,j2,j3,j4,j5,j6,itot,irat
      double precision mt,mtsq,s134
      double complex zab2,iza,izb,amp2,amp4,Xmp(2,2,2),Xpm(2,2,2)
      parameter(itot=1,irat=2)

C---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      iza(k1,k2)=cone/za(k1,k2)
      izb(k1,k2)=cone/zb(k1,k2)
C---end statement functions

      s134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      k1=j1
      k2=j2
      mtsq=mt**2
      
      do h3=1,2
      do h5=1,2
      if (h3 .eq. 1) then
        k3=j3
        k4=j4
      elseif (h3 .eq. 2) then
        k3=j4
        k4=j3
      endif
      if (h5 .eq. 1) then
        k5=j5
        k6=j6
      elseif (h5 .eq. 2) then
        k5=j6
        k6=j5
      endif
      
      amp2=(
     & +za(k2,k5)*zb(k1,k4)/zab2(k2,k3,k4,k1)**2*(
     &  -0.5d0*s134**2*(za(k1,k3)*zb(k2,k6)+za(k1,k5)*za(k2,k3)
     &      *zb(k1,k6)*zb(k2,k4)*iza(k2,k5)*izb(k1,k4))
     &  -0.5d0*s134*(
     &                   +za(k1,k5)*za(k3,k4)*zb(k2,k4)*zb(k5,k6)
     &                   -za(k3,k5)*zb(k5,k6)*zab2(k1,k3,k4,k2)
     &                   +2*za(k3,k4)*zb(k4,k6)*zab2(k1,k3,k4,k2)
     &                   +zab2(k1,k3,k4,k6)*zab2(k3,k1,k4,k2))
     &  +0.5d0*za(k3,k4)*zb(k2,k6)*zab2(k2,k5,k6,k4)*zab2(k1,k3,k4,k2)
     &  -2d0*za(k3,k4)*zb(k5,k6)*zab2(k1,k3,k4,k2)*zab2(k5,k1,k3,k4))

     & +1d0/zab2(k2,k3,k4,k1)*(
     &  -0.5d0*s134*za(k1,k3)*zb(k2,k6)*zab2(k5,k1,k3,k4)
     &  -2d0/s134*za(k3,k4)*zb(k5,k6)
     &      *zab2(k1,k3,k4,k2)*zab2(k5,k1,k3,k4)**2
     &  -0.5d0*za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k5,k6)*zab2(k5,k3,k4,k1)
     &  -0.5d0*za(k1,k5)*za(k3,k4)*zb(k2,k4)*zb(k2,k6)*zab2(k2,k5,k6,k4)
     &  +0.5d0*zab2(k1,k3,k4,k2)*zab2(k5,k1,k3,k4)
     &   *(za(k2,k3)*zb(k2,k6)+za(k3,k5)*zb(k5,k6)
     &     -2d0*za(k3,k4)*zb(k4,k6))) )/s(k3,k4)/s(k5,k6)
      
      amp4=2d0/(s(k1,k2)*s(k3,k4)*s(k5,k6)*zab2(k2,k3,k4,k1)**2)
     & *(za(k2,k5)*zb(k1,k4)*zab2(k1,k3,k4,k2)
     &  -za(k1,k5)*zb(k2,k4)*zab2(k2,k3,k4,k1))
     & *(za(k2,k3)*zb(k1,k6)*zab2(k1,k3,k4,k2)
     &  -za(k1,k3)*zb(k2,k6)*zab2(k2,k3,k4,k1))

      Xmp(h3,h5,itot)=amp2*mtsq+amp4*mtsq**2
      Xmp(h3,h5,irat)=amp4*(-1d0/6d0)

      if (docheck) call ggZZcapture('2x1x34',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              czip,amp2,amp4)

      enddo
      enddo

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5,:)=dconjg(Xmp(3-h3,3-h5,:))
      enddo
      enddo

      return
      end
