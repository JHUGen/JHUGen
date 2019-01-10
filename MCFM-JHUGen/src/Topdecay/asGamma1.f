      double precision function asGamma1(mt,besq,omsq)
C--   Author: John M. Campbell and R.K. Ellis, January 2012  
C-----Taken from formula (5) of
C-----Fermilab-PUB-12-078-T
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      double precision mt,P0,P3,PP,PM,W0,WP,wm,YW,z,omsq,om,be,
     & P0b,P3b,Pmb,PPb,Ywb,Wmb,f,besq,ddilog,GammaInfty,term4,
     & term7,term9

c     Statement functions.
      P0(z)=0.5d0*(1d0-omsq+z)
      P3(z)=0.5d0*dsqrt(1d0+omsq**2+z**2-2d0*(omsq+z+omsq*z))
      PP(z)=P0(z)+P3(z)
      PM(z)=P0(z)-P3(z)
      W0(z)=0.5d0*(1d0+omsq-z)
      WP(z)=W0(z)+P3(z)
      WM(z)=W0(z)-P3(z)
      YW(z)=0.5d0*dlog(WP(z)/WM(z))
c      YP(z)=0.5d0*dlog(PP(z)/PM(z))
c     End statement functions.
     
      f=(1d0-besq)**2+omsq*(1d0+besq)-2d0*omsq**2
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

      GammaInfty=Gf*mt**3/8d0/pi/rt2      
      if (besq .gt. 0d0) then
c      term4=Ypb*dlog(4d0*P3b**2/Ppb**2/Wpb)
      term4=(log(PPb)-log(be))*dlog(4d0*P3b**2*Wmb/(omsq*PPb**2))
      term7=
     & +(3d0-besq+11d0*besq**2-besq**3+omsq*(6d0-12d0*besq+2d0*besq**2)
     & -omsq**2*(21d0+5d0*besq)+12d0*omsq**3)*log(Ppb)
     & -(-besq+11d0*besq**2-besq**3+omsq*(-12d0*besq+2d0*besq**2)
     & -omsq**2*(5d0*besq))*log(be)
c     & -(3d0+6d0*omsq-omsq**2*21d0+12d0*omsq**3)*log(be)
      term9= 
     & +6d0*(1d0-4d0*besq+3d0*besq**2+omsq*(3d0+besq)-4d0*omsq**2)
     & *(P3b-0.5d0*(1d0-omsq))*dlog(be)
     & +3d0*(1d0-omsq)*(-4d0*besq+3d0*besq**2+omsq*(besq))*dlog(be)
c     & +3d0*(1d0-omsq)*(1d0+3d0*omsq-4d0*omsq**2)*dlog(be)
      else
      term4=log(PPb)*dlog(4d0*P3b**2*Wmb/(omsq*PPb**2))
      term7=
     & +(3d0-besq+11d0*besq**2-besq**3+omsq*(6d0-12d0*besq+2d0*besq**2)
     & -omsq**2*(21d0+5d0*besq)+12d0*omsq**3)*log(Ppb)
      term9=0d0
      endif


c--- equation for alphas*Gamma1      
      asGamma1=GammaInfty*ason2pi*Cf*(
     & 8d0*f*P0b*(ddilog(1d0-Pmb)-ddilog(1d0-Ppb)
     &  -2d0*ddilog(1d0-Pmb/Ppb)+term4
     &  +Ywb*dlog(Ppb))
     & +4d0*(1d0-besq)*((1d0-besq)**2+omsq*(1d0+besq)-4d0*omsq**2)*Ywb
     & +term7
     & +8d0*f*P3b*dlog(om/4d0/P3b**2)+term9
     & +(5d0-22d0*besq+5d0*besq**2+9d0*omsq*(1d0+besq)-6d0*omsq**2)*P3b)

      return
      end
      
