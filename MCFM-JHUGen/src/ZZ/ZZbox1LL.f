      subroutine ZZbox1LL(j1,j2,j3,j4,j5,j6,za,zb,mt,Xpp,Xmp,Xpm,Xmm,
     & Xrat)
      implicit none
      include 'types.f'
C-----Author: R.K. Ellis (September 2013)
C-----Box coefficient for LL coupling 
C-----Box D0(p1,p34,p2,mt,mt,mt,mt)
C-----Opposite off-shellnesses: so called easy box
C-----Xpp and Xmp refer to the initial state gluon polarizations
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      real(dp):: mt,mtsq,s134,s234
      integer:: j1,j2,j3,j4,j5,j6,k1,k2,k3,k4,k5,k6,h3,h5
      complex(dp):: zab2,app0,app2,app4,funcpp2,amp2,amp4,funcmp2,
     & Xpp(2,2),Xmp(2,2),Xpm(2,2),Xmm(2,2),Xrat(2,2,2,2)
c      complex(dp):: funcpp4
c--- statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      funcpp2(k1,k2,k3,k4,k5,k6)=-half*(
     & +za(k1,k3)*za(k2,k5)*za(k3,k4)*zb(k4,k2)*zb(k4,k3)*zb(k6,k1)
     & +za(k1,k3)*za(k1,k5)*za(k2,k3)*zb(k2,k1)*zb(k4,k3)*zb(k6,k1)
     & -za(k1,k2)*za(k3,k4)*za(k3,k5)*zb(k4,k2)*zb(k4,k3)*zb(k6,k1)
     & -za(k1,k5)*za(k2,k3)*zab2(k5,k3,k4,k2)*zb(k4,k1)*zb(k6,k5)
     & -4._dp*(za(k1,k5)*za(k2,k3))**2*zb(k2,k1)*zb(k4,k3)*zb(k6,k5)
     & /za(k1,k2)
     & -za(k1,k3)*za(k2,k5)*za(k3,k5)*zb(k2,k1)*zb(k4,k3)*zb(k6,k5)
     & +za(k1,k2)*za(k3,k5)**2*zb(k2,k1)*zb(k4,k3)*zb(k6,k5))
     & /(za(k1,k2)**2*s(k3,k4)*s(k5,k6))
c      funcpp4(k1,k2,k3,k4,k5,k6)=
c     & -two*za(k2,k3)*zb(k6,k1)/(s(k3,k4)*s(k5,k6)*za(k1,k2)**2)
c     & *(za(k2,k5)*zb(k4,k1)*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)
c     & -za(k1,k5)*zb(k4,k2))
      funcmp2(k1,k2,k3,k4,k5,k6)=half*(-
     & za(k5,k6)*zab2(k1,k3,k4,k2)*zab2(k3,k2,k4,k1)
     & *zb(k4,k1)*zb(k6,k2)**2/(zab2(k2,k3,k4,k1)*zb(k2,k1)**2)
     & +za(k1,k5)**2*za(k2,k3)
     & *zab2(k1,k3,k4,k2)*zab2(k2,k1,k3,k4)*zb(k6,k5)
     & /(za(k1,k2)**2*zab2(k2,k3,k4,k1)))/(s(k3,k4)*s(k5,k6))

c--- end statement function

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
      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      s234=s(k2,k3)+s(k2,k4)+s(k3,k4)

C----Helicity 1^+ 2^+
      app0=half*(two*za(k3,k2)*za(k3,k1)*za(k2,k5)*za(k1,k5)
     & +(za(k3,k5)*za(k1,k2))**2)*(s134*s234-s(k3,k4)*s(k5,k6))
     & /(za(k1,k2)**4*za(k3,k4)*za(k5,k6))
      app2=funcpp2(k1,k2,k3,k4,k5,k6)+funcpp2(k2,k1,k3,k4,k5,k6)

c      app4=funcpp4(k1,k2,k3,k4,k5,k6)+funcpp4(k2,k1,k3,k4,k5,k6)

      app4=two/(s(k3,k4)*s(k5,k6)*za(k1,k2)**2)
     & *(zb(k4,k1)*zab2(k1,k3,k4,k2)*za(k2,k5)
     &  -zb(k4,k2)*zab2(k2,k3,k4,k1)*za(k1,k5))
     & *(za(k1,k3)*zb(k6,k2)/zab2(k1,k3,k4,k2)
     &  -za(k2,k3)*zb(k6,k1)/zab2(k2,k3,k4,k1))




      Xpp(h3,h5)=app0+mtsq*app2+mtsq**2*app4
c--- contribution to rational part (coefficient of mt^4)
      Xrat(2,2,h3,h5)=app4      

C----Helicity 1^- 2^+
      amp2=funcmp2(k1,k2,k3,k4,k5,k6)+funcmp2(k1,k2,k5,k6,k3,k4)


      amp4=two/(s(j1,j2)*s(k3,k4)*s(k5,k6)*zab2(k2,k3,k4,k1)**2)
     & *(za(k2,k5)*zb(k1,k4)*zab2(k1,k3,k4,k2)
     &  -za(k1,k5)*zb(k2,k4)*zab2(k2,k3,k4,k1))
     & *(za(k2,k3)*zb(k1,k6)*zab2(k1,k3,k4,k2)
     &  -za(k1,k3)*zb(k2,k6)*zab2(k2,k3,k4,k1))

      Xmp(h3,h5)=mtsq*amp2+mtsq**2*amp4
c--- contribution to rational part (coefficient of mt^4)
      Xrat(1,2,h3,h5)=amp4      

      if (docheck) call ggZZcapture('1x34x2mp',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              czip,amp2,amp4)
      if (docheck) call ggZZcapture('1x34x2pp',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              app0,app2,app4)

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

      if (docheck) then
        call ggZZparsecheck('LL',-1,+1,Xmp,
     &   '( s56,   0, s34,   0,s156,s134;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,Xpp,
     &   '( s56,   0, s34,   0,s156,s134;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,Xpm,
     &   '( s56,   0, s34,   0,s156,s134;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,-1,Xmm,
     &   '( s56,   0, s34,   0,s156,s134;mtsq,mtsq,mtsq,mtsq)')
      endif
      
      return
      end
      
