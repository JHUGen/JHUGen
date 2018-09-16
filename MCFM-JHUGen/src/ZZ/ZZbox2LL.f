      subroutine ZZbox2LL(j1,j2,j3,j4,j5,j6,za,zb,mt,Xpp,Xmp,Xpm,Xmm,
     & Xrat)
      implicit none
      include 'types.f'
C-----Author: R.K. Ellis (September 2013)
C-----Box coefficient for a coupling
C-----Box D0(p2,p1,p34,mt,mt,mt,mt)
C-----Adjacent off-shellnesses: so-called hard box
C-----Xpp and Xmp refer to the initial state gluon polarizations
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      include 'ggZZcomputemp.f'
      real(dp):: mt,mtsq,s134,s34,s56
      integer:: k1,k2,k3,k4,k5,k6,j1,j2,j3,j4,j5,j6,h3,h5
      complex(dp):: zab2,funcapp2,
     & amp2,amp0,amp4,app2,app4,Xswap(2,2),
     & Xpp(2,2),Xmp(2,2),Xpm(2,2),Xmm(2,2),Xrat(2,2,2,2)
c----- statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      funcapp2(k1,k2,k3,k4,k5,k6)=zb(k2,k1)*zb(k4,k3)
     & *(za(k1,k5)*zab2(k1,k3,k4,k6)*zab2(k3,k5,k6,k2)**2
     &  -za(k1,k3)*za(k5,k6)*zb(k6,k2)**2
     & *(zab2(k1,k3,k4,k5)*za(k5,k3)+zab2(k1,k3,k4,k6)*za(k6,k3)))
     & /(two*s34*s56*za(k1,k2)*zab2(k1,k3,k4,k2)**2)
c----- end statement functions
      mtsq=mt**2
      k1=j1
      k2=j2
      s134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      s34=s(j3,j4)
      s56=s(j5,j6)
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

C----Helicity 1^+ 2^+ (of order mtsq)
c      app0=czip
c      app2=zb(k2,k1)*zb(k4,k3)
c     & *(za(k1,k5)*zab2(k1,k3,k4,k6)*zab2(k3,k5,k6,k2)**2
c     &  -za(k1,k3)*za(k5,k6)*zb(k6,k2)**2
c     & *(zab2(k1,k3,k4,k5)*za(k5,k3)+zab2(k1,k3,k4,k6)*za(k6,k3)))
c     & /(two*s34*s56*za(k1,k2)*zab2(k1,k3,k4,k2)**2)
c     & +zb(k1,k2)*zb(k6,k5)
c     & *(za(k2,k3)*zab2(k2,k5,k6,k4)*zab2(k5,k3,k4,k1)**2
c     &  +za(k2,k5)*za(k3,k4)*zb(k4,k1)**2
c     & *(zab2(k2,k1,k3,k4)*za(k4,k5)+zab2(k2,k1,k4,k3)*za(k3,k5)))
c     & /(two*s34*s56*za(k2,k1)*zab2(k2,k3,k4,k1)**2)
c      write(6,*) 'app2_old',app2

      app2=funcapp2(k1,k2,k3,k4,k5,k6)+funcapp2(k2,k1,k5,k6,k3,k4)
      app4=two/(s34*s56*za(k1,k2)**2)
     & *(zb(k4,k1)*zab2(k1,k3,k4,k2)*za(k2,k5)
     &  -zb(k4,k2)*zab2(k2,k3,k4,k1)*za(k1,k5))
     & *(za(k1,k3)*zb(k6,k2)/zab2(k1,k3,k4,k2)
     &  -za(k2,k3)*zb(k6,k1)/zab2(k2,k3,k4,k1))

  

      Xpp(h3,h5)=mtsq*app2+mtsq**2*app4
c--- contribution to rational part (coefficient of mt^4)
      Xrat(2,2,h3,h5)=app4      


      if (computemp) then
C----Helicity 1^- 2^+
      amp0=half*s(k1,k2)*s134
     & *((zab2(k2,k1,k3,k4)*zab2(k5,k3,k4,k1))**2
     &  +(s134*za(k2,k5)*zb(k1,k4))**2)
     & /(zb(k3,k4)*za(k5,k6)*zab2(k2,k3,k4,k1)**4)

!===== old formula
!      amp2=(za(k1,k3)*za(k1,k5)*zb(k2,k4)*zb(k2,k6) 
!     & +( 
!     & zab2(k1,k3,k4,k2)
!     & *(za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6)
!     &  -za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)
!     &  +za(k1,k3)*za(k3,k5)*zb(k1,k6)*zb(k3,k4)
!     &  -za(k2,k3)*za(k3,k5)*zb(k2,k6)*zb(k3,k4)
!     &  -za(k1,k5)*za(k2,k3)*zb(k1,k4)*zb(k2,k6))

!     & +half*s134*
!     & (two*za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4)
!     &     +za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k2,k6)
!     &     +za(k1,k3)*za(k3,k5)*zb(k2,k6)*zb(k3,k4)
!     &     -za(k1,k5)*za(k2,k3)*zb(k2,k4)*zb(k2,k6)
!     &     -za(k1,k5)*za(k3,k4)*zb(k2,k4)*zb(k4,k6))
!     & )*izab2(k2,k3,k4,k1)

!     & +s134*(  
!     & -za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)*s134
!     & -za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*s134!
!
!     & +3._dp/two*zab2(k1,k3,k4,k2)*(
!     & +za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6)
!     & +za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(k4,k6)
!     & -za(k1,k5)*za(k2,k3)*zb(k1,k4)*zb(k1,k6)
!     & -za(k2,k3)*za(k3,k5)*zb(k1,k6)*zb(k3,k4)
!     & -two*za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6))
!     & )*izab2(k2,k3,k4,k1)**2!

!     & +4._dp*za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)
!     & *zab2(k1,k3,k4,k2)*s134**2
!     & *izab2(k2,k3,k4,k1)**3)
!     & /(s34*s56)

      amp2=(s134*za(k2,k5)*zab2(k1,k3,k4,k6)*
     -     zab2(k3,k1,k4,k2)*zb(k4,k1))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)**2) + 
     -  (za(k3,k4)*zab2(k1,k3,k4,k2)*
     -     zab2(k2,k1,k3,k4)**2*zab2(k5,k3,k4,k1)*
     -     zb(k6,k1))/(s34*s56*zab2(k2,k3,k4,k1)**3) - 
     -  (s134**2*za(k1,k5)*za(k2,k3)*zb(k4,k2)*
     -     zb(k6,k1))/(2.*s34*s56*zab2(k2,k3,k4,k1)**2)
     -    - (za(k2,k3)*zab2(k1,k3,k4,k2)*
     -     zab2(k5,k1,k3,k4)*zb(k6,k2))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)) + 
     -  (s134*za(k1,k3)*zab2(k2,k1,k3,k4)*
     -     zab2(k5,k3,k4,k1)*zb(k6,k2))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)**2) - 
     -  (za(k2,k5)*za(k3,k4)*zab2(k1,k3,k4,k2)*
     -     zab2(k2,k1,k3,k4)*zb(k4,k1)*zb(k6,k2))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)**2) + 
     -  (za(k1,k5)*za(k3,k4)*zab2(k2,k1,k3,k4)*
     -     zb(k4,k2)*zb(k6,k2))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)) + 
     -  (za(k2,k3)*zab2(k1,k3,k4,k2)*zab2(k2,k1,k3,k4)*
     -     zab2(k5,k3,k4,k1)**2*zb(k6,k5))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)**3) + 
     -  (za(k2,k5)*za(k3,k4)*zab2(k1,k3,k4,k2)*
     -     zab2(k2,k1,k3,k4)*zab2(k5,k3,k4,k1)*
     -     zb(k4,k1)*zb(k6,k5))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)**3) - 
     -  (2*s134*za(k2,k5)**2*za(k3,k4)*
     -     zab2(k1,k3,k4,k2)*zb(k4,k1)**2*zb(k6,k5))/
     -   (s34*s56*zab2(k2,k3,k4,k1)**3) - 
     -  (za(k1,k3)*za(k1,k5)*zab2(k5,k3,k4,k1)*
     -     zb(k4,k2)*zb(k6,k5))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)) + 
     -  (s134*za(k1,k5)*za(k2,k5)*za(k3,k4)*zb(k4,k1)*
     -     zb(k4,k2)*zb(k6,k5))/
     -   (2.*s34*s56*zab2(k2,k3,k4,k1)**2)

      amp4=two/(s(k1,k2)*s(k3,k4)*s(k5,k6)*zab2(k2,k3,k4,k1)**2)
     & *(za(k2,k5)*zb(k1,k4)*zab2(k1,k3,k4,k2)
     &  -za(k1,k5)*zb(k2,k4)*zab2(k2,k3,k4,k1))
     & *(za(k2,k3)*zb(k1,k6)*zab2(k1,k3,k4,k2)
     &  -za(k1,k3)*zb(k2,k6)*zab2(k2,k3,k4,k1))

        Xmp(h3,h5)=amp0+mtsq*amp2+mtsq**2*amp4
c--- contribution to rational part (coefficient of mt^4)
        Xrat(1,2,h3,h5)=amp4  

      else
        Xmp(h3,h5)=czip
        Xrat(1,2,h3,h5)=czip
      endif
      
      if (docheck) call ggZZcapture('2x1x34pp',h3,h5,j1,j2,j3,j4,j5,j6,
     &                              czip,app2,app4)

      enddo
      enddo

c--- swapped case: need to interchange helicity labels
      if (j3 == 5) then
        Xswap(:,:)=Xmp(:,:)
        Xmp(1,2)=Xswap(2,1)
        Xmp(2,1)=Xswap(1,2)
        Xswap(:,:)=Xpp(:,:)
        Xpp(1,2)=Xswap(2,1)
        Xpp(2,1)=Xswap(1,2)        
        Xswap(:,:)=Xrat(1,2,:,:)
        Xrat(1,2,1,2)=Xswap(2,1)
        Xrat(1,2,2,1)=Xswap(1,2)
        Xswap(:,:)=Xrat(2,2,:,:)
        Xrat(2,2,1,2)=Xswap(2,1)
        Xrat(2,2,2,1)=Xswap(1,2)
      endif
      
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
        if (j3 == 3) then ! normal case: 1,2,3,4,5,6
        call ggZZparsecheck('LL',-1,+1,Xmp,
     &   '( s34,   0,   0, s56,s134, s12;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,Xpp,
     &   '( s34,   0,   0, s56,s134, s12;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,Xpm,
     &   '( s34,   0,   0, s56,s134, s12;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,-1,Xmm,
     &   '( s34,   0,   0, s56,s134, s12;mtsq,mtsq,mtsq,mtsq)')
        endif
        if (j3 == 5) then ! swapped case: 1,2,5,6,3,4
        call ggZZparsecheck('LL',-1,+1,Xmp,
     &   '( s56,   0,   0, s34,s156, s12;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,Xpp,
     &   '( s56,   0,   0, s34,s156, s12;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,Xpm,
     &   '( s56,   0,   0, s34,s156, s12;mtsq,mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,-1,Xmm,
     &   '( s56,   0,   0, s34,s156, s12;mtsq,mtsq,mtsq,mtsq)')
        endif
      endif
      
      return
      end
