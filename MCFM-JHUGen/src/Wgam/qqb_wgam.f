      subroutine qqb_wgam(p,msq)
      implicit none
      include 'types.f'
C-----Author Keith Ellis, September 2002
C----- updated: John Campbell, August 2011 (anomalous couplings)
c----Matrix element for W gam production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5)
C For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'nwz.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qbq,qqb,fac
      complex(dp):: agamtree

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(5,p,za,zb)
      fac=aveqq*2._dp*xn*gwsq**2*esq

      if (nwz == -1) then 
         qbq=fac*(abs(agamtree(1,2,3,4,5,za,zb,-1))**2
     &           +abs(agamtree(1,2,3,4,5,za,zb,+1))**2)
         qqb=fac*(abs(agamtree(2,1,3,4,5,za,zb,-1))**2
     &           +abs(agamtree(2,1,3,4,5,za,zb,+1))**2)
      elseif (nwz == +1) then 
         qbq=fac*(abs(agamtree(2,1,4,3,5,zb,za,-1))**2
     &           +abs(agamtree(2,1,4,3,5,zb,za,+1))**2)
         qqb=fac*(abs(agamtree(1,2,4,3,5,zb,za,-1))**2
     &           +abs(agamtree(1,2,4,3,5,zb,za,+1))**2)
      endif

      do j=-nf,nf
      do k=-nf,nf 
         if ((j > 0) .and. (k < 0)) then            
            msq(j,k)=Vsq(j,k)*qqb
       
         elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*qbq
           
         endif
         
      enddo
      enddo
            
      return
      end

      function agamtree(p1,p2,p3,p4,p5,za,zb,hgamma)
      implicit none
      include 'types.f'
      complex(dp):: agamtree
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      real(dp):: xfac
      complex(dp):: prp12,prp34
      integer:: p1,p2,p3,p4,p5,hgamma
      prp34=s(p3,p4)/cplx2((s(p3,p4)-wmass**2),wmass*wwidth)

c--- apply a dipole form factor to anomalous couplings, with power two (only if tevscale > 0)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s(p1,p2)/(tevscale*1d3)**2)**2
      else
        xfac=1._dp
      endif
      xdelk_g=xfac*delk_g
      xlambda_g=xfac*lambda_g
      
      if (zerowidth) then
c--- zerowidth: no final state radiation, so we can set prp12 to zero
        prp12=czip
      else
c--- otherwise, usual Breit-Wigner form
        prp12=s(p1,p2)/cplx2((s(p1,p2)-wmass**2),wmass*wwidth)
      endif


c---  c.f. Eqs.(4.4),(4.5) of hep-ph/9803250 (multiplied by -i)
c---       for the terms proportional to prp34
      if (hgamma == -1) then
        agamtree=-zb(p2,p4)**2/(s(p1,p2)-s(p3,p4))*(
     &    Qu*(za(p2,p5)/(zb(p4,p3)*zb(p1,p5))*prp34
     &       -za(p4,p5)/(zb(p1,p2)*zb(p3,p5))*prp12)
     &   +Qd*(za(p1,p5)/(zb(p4,p3)*zb(p2,p5))*prp34
     &       +za(p4,p5)/(zb(p1,p2)*zb(p3,p5))*prp12))
      elseif (hgamma == +1) then
        agamtree=-za(p1,p3)**2/(s(p1,p2)-s(p3,p4))*(
     &    Qd*(zb(p1,p5)/(za(p3,p4)*za(p2,p5))*prp34
     &       +zb(p4,p5)/(za(p2,p1)*za(p3,p5))*prp12)
     &   +Qu*(zb(p2,p5)/(za(p3,p4)*za(p1,p5))*prp34
     &       -zb(p4,p5)/(za(p2,p1)*za(p3,p5))*prp12))
      endif
 
c--- additional anomalous coupling component, Eqs. (7)-(9) of hep-ph/0002138
c--- note partial fractioning of W propagators to go beyond zerowidth approx.
      if     (hgamma == -1) then
       agamtree=agamtree+(prp34-prp12)*(Qu-Qd)
     &  *za(p5,p3)/(2._dp*s(p3,p4)*(s(p1,p2)-s(p3,p4))*za(p4,p3))
     &  *((xdelk_g+xlambda_g)*zb(p2,p4)*za(p1,p5)*za(p4,p3)
     &    +xlambda_g*zb(p2,p5)*za(p5,p1)*za(p3,p5))
      elseif (hgamma == +1) then
       agamtree=agamtree+(prp34-prp12)*(Qd-Qu)
     &  *zb(p5,p4)/(2._dp*s(p3,p4)*(s(p1,p2)-s(p3,p4))*zb(p3,p4))
     &  *((xdelk_g+xlambda_g)*za(p1,p3)*zb(p2,p5)*zb(p3,p4)
     &    +xlambda_g*za(p1,p5)*zb(p5,p2)*zb(p4,p5))
       endif
       
      return
      end

c      function ubdmsq(i1,i2,i3,i4,i5)
c      implicit none
c      include 'types.f'
c      real(dp):: ubdmsq
c      
Cc     Matrix element for 
Cc     ub(-p1)+d(-p2)=e-(p3)+nu~(p4)+gamma(p5)
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'ewcouple.f'
c      include 'ewcharge.f'
c      include 'ckm.f'
c      include 'sprods_com.f'
c      include 'zprods_com.f'
c      integer:: i1,i2,i3,i4,i5
c      real(dp):: fac,propsq12,propsq34,propsq1234 
cc     ans=
cc      +2*xn*gw^4*ee^2*(s24^2+s13^2)*[2*p1.p5+2*p2.p5]^-2*( 
cc       s12*[Qu-Qd]^2*s45/s35/[s12-mw^2]^2
c
cc       +s34*(Qu/s15+Qd/s25)^2*s15*s25/[s34-mw^2]^2
c
Cc       +((s23*s45-s24*s35)*s15-(s13*s45-s14*s35)*s25)
Cc       *[Qu-Qd]*(Qu/s15+Qd/s25)/s35/[s12-mw^2]/[s34-mw^2])
c      propsq12=(s(i1,i2)-wmass**2)**2+(wmass*wwidth)**2
c      propsq34=(s(i3,i4)-wmass**2)**2+(wmass*wwidth)**2
c      propsq1234=
c     & (s(i1,i2)-wmass**2)*(s(i3,i4)-wmass**2)+(wmass*wwidth)**2
c      ubdmsq=
c     & +s(i1,i2)*(Qu-Qd)**2*s(i4,i5)/(s(i3,i5)*propsq12)
c     &     +s(i3,i4)*(Qu/s(i1,i5)+Qd/s(i2,i5))**2
c     &      *s(i1,i5)*s(i2,i5)/propsq34
c     & +((s(i2,i3)*s(i4,i5)-s(i2,i4)*s(i3,i5))*s(i1,i5)
c     &  -(s(i1,i3)*s(i4,i5)-s(i1,i4)*s(i3,i5))*s(i2,i5))
c     & *(Qu-Qd)*(Qu/s(i1,i5)+Qd/s(i2,i5))
c     & /(s(i3,i5)*propsq1234)
CCCC Epsilon piece missing! 

c      ubdmsq=(s(i2,i4)**2+s(i1,i3)**2)/(s(i1,i5)+s(i2,i5))**2*ubdmsq
c      fac=aveqq*2._dp*xn*gwsq**2*esq
c      ubdmsq=fac*ubdmsq
c      return
c      end
