      subroutine ZZtri1_34LL(j1,j2,j3,j4,j5,j6,za,zb,mt,Xpp,Xmp,Xpm,Xmm,
     & Xrat)
      implicit none
      include 'types.f'
C-----Author: R.K. Ellis (September 2013)
C-----Trianglecoefficient for LL coupling
C-----Triangle C0(p1,p2,mt,mt,mt)
C-----Xpp and Xmp refer to the initial state gluon polarizations

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      include 'ggZZcomputemp.f'
      real(dp):: mt,mtsq,t
      integer:: j1,j2,j3,j4,j5,j6,k1,k2,k3,k4,k5,k6,h3,h5
      complex(dp):: Xpp(2,2),Xmp(2,2),Xmm(2,2),Xpm(2,2),zab2,
     & Funcpp_2,Funcpp_0,Funcmp_0,Funcmp_2,app0,app2,amp0,amp2,
     & Xrat(2,2,2,2)
c--- statement functions
      t(k1,k2,k3)=s(k1,k2)+s(k2,k3)+s(k1,k3)
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)

C----Functions for the pp amplitude
      Funcpp_0(k1,k2,k3,k4,k5,k6)=
     & -((za(k1,k5)**2*za(k2,k3)**2+za(k1,k3)**2*za(k2,k5)**2)
     & *(s(k1,k3)+s(k1,k4))/(two*za(k1,k2)**4*za(k3,k4)*za(k5,k6)))


      Funcpp_2(k1,k2,k3,k4,k5,k6)=(
     & +((za(k1,k3)*za(k2,k5)*zab2(k1,k3,k4,k6)*zb(k4,k1))
     & /(za(k1,k2)**3))

     & -(za(k1,k3)**2*za(k2,k5)*zab2(k1,k3,k4,k6)*zb(k2,k1)*zb(k4,k3))
     & /(za(k1,k2)**3*zab2(k1,k3,k4,k2))

     & +(za(k1,k3)*zab2(k2,k1,k3,k4)*zab2(k5,k3,k4,k1)*zb(k6,k1))
     & /(za(k1,k2)**2*zab2(k2,k3,k4,k1))

     & -(za(k1,k3)**2*za(k1,k5)*zb(k2,k1)*zb(k4,k3)*zb(k6,k1))
     & /(za(k1,k2)**2*zab2(k1,k3,k4,k2))

     & +(2*za(k1,k3)*za(k1,k5)*za(k3,k4)*zb(k4,k1)*zb(k4,k3)*zb(k6,k1))
     & /(za(k1,k2)**2*zab2(k1,k3,k4,k1))

     & +(2*za(k1,k3)*zab2(k5,k3,k4,k1)*zb(k4,k1)*zb(k6,k2))
     & /(za(k1,k2)*zab2(k1,k3,k4,k1))

     & -(za(k1,k3)*za(k3,k5)*zab2(k1,k3,k4,k1)*zb(k4,k3)*zb(k6,k2))
     & /(za(k1,k2)**2*zab2(k1,k3,k4,k2))

     & -(za(k1,k3)**2*za(k1,k5)*zab2(k1,k3,k4,k1)
     &  *zab2(k2,k5,k6,k2)*zb(k4,k3)*zb(k6,k2))
     & /(za(k1,k2)**3*zab2(k1,k3,k4,k2)**2)

     & +(za(k1,k3)*za(k1,k5)*za(k2,k3)*zb(k2,k1)*zb(k4,k3)*zb(k6,k2))
     & /(za(k1,k2)**2*zab2(k1,k3,k4,k2))

     & +(za(k1,k3)*za(k1,k5)*zab2(k1,k3,k4,k1)
     & *zab2(k2,k3,k4,k2)*zb(k6,k4))
     & /(za(k1,k2)**3*zab2(k1,k3,k4,k2))

     & +(za(k2,k3)*zab2(k5,k3,k4,k1)**2*zb(k4,k1)*zb(k6,k5))/
     & (za(k1,k2)*zab2(k2,k3,k4,k1)**2))
     & /(s(k3,k4)*s(k5,k6))

C----Functions for the mp amplitude
      Funcmp_0(k1,k2,k3,k4,k5,k6)=
     & +half*(za(k1,k5)**2*za(k3,k4)*zb(k4,k1)**2)
     & /(za(k5,k6)*zab2(k2,k3,k4,k1)**2)

     & -(t(k1,k3,k4)*za(k1,k2)*za(k1,k3)*zb(k4,k1)
     & *zab2(k5,k3,k4,k1)**2)
     & /(za(k5,k6)*(s(k1,k3)+s(k1,k4))*zab2(k2,k3,k4,k1)**3)

     &+half*(za(k3,k5)*(s(k1,k3)+s(k1,k4))*zab2(k5,k3,k4,k1)*zb(k4,k2))
     & /(za(k5,k6)*zb(k2,k1)*zab2(k2,k3,k4,k1)**2)

     & +(za(k1,k3)*zab2(k5,k3,k4,k1)**2*zb(k4,k2))
     & /(za(k5,k6)*zb(k2,k1)*zab2(k2,k3,k4,k1)**2)

     & -(za(k1,k5)*za(k3,k4)*zb(k4,k1)**2
     & *zab2(k5,k3,k4,k1)*zab2(k1,k3,k4,k2))
     & /(za(k5,k6)*zb(k2,k1)*(s(k1,k3)+s(k1,k4))*zab2(k2,k3,k4,k1)**2)

     &+half*(za(k3,k5)*zb(k4,k1)*(s(k1,k3)+s(k1,k4))*zab2(k5,k3,k4,k2))
     & /(za(k5,k6)*zb(k2,k1)*zab2(k2,k3,k4,k1)**2)

     & +half*(t(k1,k3,k4)*za(k1,k5)**2*zb(k4,k1)**2)
     & /(za(k5,k6)*zab2(k2,k3,k4,k1)**2*zb(k4,k3))

     & -(t(k1,k3,k4)**2*za(k1,k2)*za(k2,k5)*zb(k4,k1)**2
     & *zab2(k5,k3,k4,k1))
     & /(za(k5,k6)*zab2(k2,k3,k4,k1)**4*zb(k4,k3))

      Funcmp_2(k1,k2,k3,k4,k5,k6)=
     &  ((za(k1,k3)*za(k1,k5)*zb(k2,k1)*zb(k4,k1)*zab2(k1,k3,k4,k6))
     & /(zab2(k2,k3,k4,k1))

     & +(za(k1,k3)*zb(k4,k1)*zab2(k5,k3,k4,k1)
     & *(s(k1,k2)+s(k2,k3)+s(k2,k4))*zab2(k1,k3,k4,k6))
     & /(zab2(k2,k3,k4,k1)**2)

     & -(za(k1,k3)*za(k1,k5)*zb(k4,k1)
     & *(s(k1,k2)+s(k2,k3)+s(k2,k4))
     & *zb(k6,k1)*s(k3,k4))/(zab2(k2,k3,k4,k1)**2)

     & -(za(k1,k2)*za(k3,k4)*zb(k4,k1)**2*zab2(k5,k3,k4,k1)
     & *zab2(k1,k3,k4,k2)
     & *(t(k1,k3,k4)*s(k1,k2)
     & +(s(k1,k3)+s(k1,k4))*(s(k1,k2)+s(k2,k3)+s(k2,k4)))
     & *zb(k6,k1))
     & /(zb(k2,k1)*(s(k1,k3)+s(k1,k4))*zab2(k2,k3,k4,k1)**3)

     & +(t(k1,k3,k4)*za(k1,k3)*za(k1,k5)*zb(k4,k1)*zb(k6,k2))
     & /(zab2(k2,k3,k4,k1))

     & +(t(k2,k3,k4)*za(k1,k3)*za(k1,k5)*zb(k4,k1)*zb(k6,k2))
     & /(zab2(k2,k3,k4,k1))

     & -(za(k1,k2)**2*za(k3,k4)*zb(k4,k1)**2
     & *zab2(k5,k3,k4,k1)*zab2(k1,k3,k4,k2)*zb(k6,k2))
     & /((s(k1,k3)+s(k1,k4))*zab2(k2,k3,k4,k1)**2)

     & +(two*za(k1,k3)*za(k1,k5)**2
     & *zb(k2,k1)*zb(k4,k1)*zb(k6,k5)*s(k3,k4))
     & /((s(k1,k3)+s(k1,k4))*zab2(k2,k3,k4,k1))

     & +(za(k1,k3)**2*za(k1,k5)**2*zb(k2,k1)*zb(k4,k3)*zb(k6,k5))
     & /(za(k1,k2)*(s(k1,k3)+s(k1,k4))))
     & /(s(k1,k2)*s(k3,k4)*s(k5,k6))

c--- end statement functions

      mtsq=mt**2

      k1=j1
      k2=j2

      do h3=1,2
         if (h3 ==1) then
             k3=j3
             k4=j4
         elseif (h3 ==2) then
             k3=j4
             k4=j3
         endif
      do h5=1,2
         if (h5 ==1) then
             k5=j5
             k6=j6
         elseif (h5 ==2) then
             k5=j6
             k6=j5
         endif

         app0=Funcpp_0(k1,k2,k3,k4,k5,k6)
         app2=Funcpp_2(k1,k2,k3,k4,k5,k6)
c         write(6,*) 'app0',app0
c         write(6,*) 'app2',app2
         Xpp(h3,h5)=app0+mtsq*app2
c--- contribution to rational part (coefficient of -mt^2)
         Xrat(2,2,h3,h5)=-app2

         if (computemp) then
           amp0=Funcmp_0(k1,k2,k3,k4,k5,k6)
           amp2=Funcmp_2(k1,k2,k3,k4,k5,k6)
c           write(6,*) 'amp0',amp0
c           write(6,*) 'amp2',amp2
           Xmp(h3,h5)=amp0+mtsq*amp2
c--- contribution to rational part (coefficient of -mt^2)
           Xrat(1,2,h3,h5)=-amp2
         else
           Xmp(h3,h5)=czip
           Xrat(1,2,h3,h5)=czip
         endif

      if (docheck) call ggZZcapture('1x34pp',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              app0,app2,czip)

      enddo
      enddo

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5)=conjg(Xmp(3-h3,3-h5))
      Xmm(h3,h5)=conjg(Xpp(3-h3,3-h5))
      Xrat(2,1,h3,h5)=conjg(Xrat(1,2,3-h3,3-h5))
      Xrat(1,1,h3,h5)=conjg(Xrat(2,2,3-h3,3-h5))
      enddo
      enddo

      return
      end

