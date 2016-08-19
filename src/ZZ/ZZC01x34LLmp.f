      subroutine ZZC01x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      implicit none
c--- Author: J. M. Campbell, October 2013
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      integer k1,k2,k3,k4,k5,k6
      integer h3,h5,j1,j2,j3,j4,j5,j6,itot,irat
      double precision mt,mtsq,t,s134,s234
      double complex zab2,amp2,Xmp(2,2,2),Xpm(2,2,2)
      parameter(itot=1,irat=2)

C---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
C---end statement functions

      s134=t(j1,j3,j4)
      s234=t(j2,j3,j4)
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

c---- NB: amp0=0      
c      amp0=1d0/zb(k5,k6)/zab2(k2,k3,k4,k1)**2*(
c     & -0.5d0*(s234+s(k3,k4))*za(k2,k3)**2*zb(k2,k6)**2/za(k3,k4)
c     & -s234*za(k2,k3)*zb(k2,k4)*zb(k2,k6)*zab2(k2,k3,k4,k6)
c     &  /zab2(k2,k3,k4,k2)
c     & +za(k1,k3)*zb(k2,k4)/za(k1,k2)*zab2(k2,k3,k4,k6)**2
c     & -0.5d0*s(k3,k4)*za(k2,k3)*zb(k2,k6)/za(k3,k4)*zab2(k3,k2,k4,k6)
c     & +0.5d0*za(k1,k3)*zb(k4,k6)/za(k1,k2)
c     &  *zab2(k2,k3,k4,k2)*zab2(k2,k3,k4,k6)
c     & -0.5d0*za(k2,k3)*zb(k2,k6)/za(k3,k4)
c     &  *zab2(k2,k3,k4,k2)*zab2(k3,k2,k4,k6)  
c     & -za(k2,k3)*zb(k2,k6)*za(k2,k3)*zb(k3,k4)/za(k1,k2)
c     &  *zab2(k2,k3,k4,k6)*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k2)   
c     & +0.5d0*za(k2,k3)*zb(k4,k6)/za(k1,k2)
c     &  *zab2(k1,k3,k4,k6)*zab2(k2,k3,k4,k2)
c     & -0.5d0*zb(k2,k4)*zab2(k2,k3,k4,k6)*zab2(k3,k2,k4,k6))

      amp2=1d0/s(k1,k2)/s(k3,k4)/s(k5,k6)*(
      
     & -2d0*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)**3*s134
     &  *zb(k1,k4)**2*zb(k1,k6)*za(k2,k5)*za(k4,k3)*za(k2,k1)

     & +zb(k1,k4)/zab2(k2,k3,k4,k1)**2*(
     &  +zab2(k2,k5,k6,k2)*za(k1,k3)*zab2(k5,k3,k4,k1)*zab2(k1,k3,k4,k6)
     &  +zab2(k2,k5,k6,k2)*za(k1,k3)*s(k3,k4)*zb(k1,k6)*za(k1,k5)
     &  +zb(k1,k4)*zb(k2,k6)*za(k1,k2)**2*za(k3,k4)*zab2(k5,k3,k4,k1)
     &   *zab2(k1,k3,k4,k2)/zab2(k1,k3,k4,k1)
     &  +zb(k1,k4)*zb(k1,k6)*za(k2,k5)/zb(k1,k2)*za(k3,k4)
     &   *zab2(k1,k3,k4,k2)**2
     &  +2d0*za(k3,k4)*zb(k1,k4)*zb(k1,k6)*za(k1,k5)*za(k1,k2)
     &   *zab2(k1,k3,k4,k2)*s134/zab2(k1,k3,k4,k1))
   
     & +1d0/zab2(k2,k3,k4,k1)*zb(k1,k4)*(
     &  +(s234+s(k3,k4))*zb(k2,k6)*za(k1,k3)*za(k1,k5)
     &  +zb(k1,k6)*za(k1,k3)*za(k1,k5)*zab2(k1,k3,k4,k2)
     &  -zb(k1,k4)*zb(k1,k6)*za(k1,k5)*za(k3,k4)*zab2(k1,k3,k4,k2)**2
     &   /zb(k1,k2)/zab2(k1,k3,k4,k1)
     &  -2d0*za(k3,k4)*zb(k3,k4)*zb(k1,k2)*za(k1,k3)*za(k1,k5)**2
     &   *zb(k6,k5)/zab2(k1,k3,k4,k1))
   
     & +zb(k1,k2)*zb(k3,k4)*za(k1,k3)**2*za(k1,k5)**2
     &   /za(k1,k2)*zb(k6,k5)/zab2(k1,k3,k4,k1) )
  
      Xmp(h3,h5,itot)=amp2*mtsq
      Xmp(h3,h5,irat)=amp2*(+0.5d0)

      if (docheck) call ggZZcapture('1x34',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              czip,amp2,czip)

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
