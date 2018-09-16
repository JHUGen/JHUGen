      function asGamma1(mt,besq,omsq)
      implicit none
      include 'types.f'
      real(dp):: asGamma1
C--   Author: John M. Campbell and R.K. Ellis, January 2012  
C-----Taken from formula (5) of
C-----Fermilab-PUB-12-078-T
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: mt,P0,P3,PP,PM,W0,WP,wm,YW,z,omsq,om,be,
     & P0b,P3b,Pmb,PPb,Ywb,Wmb,f,besq,ddilog,GammaInfty,term4,
     & term7,term9

c     Statement functions.
      P0(z)=0.5_dp*(1._dp-omsq+z)
      P3(z)=0.5_dp*sqrt(1._dp+omsq**2+z**2-2._dp*(omsq+z+omsq*z))
      PP(z)=P0(z)+P3(z)
      PM(z)=P0(z)-P3(z)
      W0(z)=0.5_dp*(1._dp+omsq-z)
      WP(z)=W0(z)+P3(z)
      WM(z)=W0(z)-P3(z)
      YW(z)=0.5_dp*log(WP(z)/WM(z))
c      YP(z)=0.5_dp*log(PP(z)/PM(z))
c     End statement functions.
     
      f=(1._dp-besq)**2+omsq*(1._dp+besq)-2._dp*omsq**2
      om=sqrt(omsq)
      be=sqrt(besq)
      P0b=P0(besq)
      P3b=P3(besq)
      Pmb=PM(besq)
      Wmb=WM(besq)
c      WPb=WP(besq)
      PPb=PP(besq)
c      Ypb=YP(besq)
      Ywb=YW(besq)

      GammaInfty=Gf*mt**3/8._dp/pi/rt2      
      if (besq > 0._dp) then
c      term4=Ypb*log(4._dp*P3b**2/Ppb**2/Wpb)
      term4=(log(PPb)-log(be))*log(4._dp*P3b**2*Wmb/(omsq*PPb**2))
      term7=
     & +(3._dp-besq+11._dp*besq**2-besq**3+omsq*(6._dp-12._dp*besq+2._dp*besq**2)
     & -omsq**2*(21._dp+5._dp*besq)+12._dp*omsq**3)*log(Ppb)
     & -(-besq+11._dp*besq**2-besq**3+omsq*(-12._dp*besq+2._dp*besq**2)
     & -omsq**2*(5._dp*besq))*log(be)
c     & -(3._dp+6._dp*omsq-omsq**2*21._dp+12._dp*omsq**3)*log(be)
      term9= 
     & +6._dp*(1._dp-4._dp*besq+3._dp*besq**2+omsq*(3._dp+besq)-4._dp*omsq**2)
     & *(P3b-0.5_dp*(1._dp-omsq))*log(be)
     & +3._dp*(1._dp-omsq)*(-4._dp*besq+3._dp*besq**2+omsq*(besq))*log(be)
c     & +3._dp*(1._dp-omsq)*(1._dp+3._dp*omsq-4._dp*omsq**2)*log(be)
      else
      term4=log(PPb)*log(4._dp*P3b**2*Wmb/(omsq*PPb**2))
      term7=
     & +(3._dp-besq+11._dp*besq**2-besq**3+omsq*(6._dp-12._dp*besq+2._dp*besq**2)
     & -omsq**2*(21._dp+5._dp*besq)+12._dp*omsq**3)*log(Ppb)
      term9=0._dp
      endif


c--- equation for alphas*Gamma1      
      asGamma1=GammaInfty*ason2pi*Cf*(
     & 8._dp*f*P0b*(ddilog(1._dp-Pmb)-ddilog(1._dp-Ppb)
     &  -2._dp*ddilog(1._dp-Pmb/Ppb)+term4
     &  +Ywb*log(Ppb))
     & +4._dp*(1._dp-besq)*((1._dp-besq)**2+omsq*(1._dp+besq)-4._dp*omsq**2)*Ywb
     & +term7
     & +8._dp*f*P3b*log(om/4._dp/P3b**2)+term9
     & +(5._dp-22._dp*besq+5._dp*besq**2+9._dp*omsq*(1._dp+besq)-6._dp*omsq**2)*P3b)

      return
      end
      
