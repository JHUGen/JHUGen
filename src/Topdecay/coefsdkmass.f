      subroutine coefsdkmass(s12,mt,mb,ct,c0L,c0R,c1L,C1R)
      implicit none
C     Authors: John Campbell and R.Keith Ellis, April 2012
C     This gives the coefficients of the 
C     virtual corrections, C0L,C0R,C1L,C1R, 
C     and if includect is set to .true. 
C     the integrated counterterm for top semi-leptonic decay, ct
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'includect.f'
      double precision s12,mb,mt,ct,
     & rsq,omrsq,eta,wlog,rlog,mulog,lambda,
     & rb,rbsq,rwsq,Zt,Zb,mtsq,mbsq,lnrbsq,ddilog
      double complex B0f1,B0f2,C0(-2:0),c0L,c0R,c1L,C1R,s12log,lnrat,
     & qlI3,qlI2
      double precision P0,P3,PP,PM,W0,WP,WM,YP,YW,z,om,omsq
      double precision P0b,P3b,PMb,PPb,Ypb,Ywb,be,besq
c     Statement functions.
      P0(z)=0.5d0*(1d0-om**2+z)
      P3(z)=0.5d0*sqrt(1d0+om**4+z**2-2d0*(om**2+z+om**2*z))
      PP(z)=P0(z)+P3(z)
      PM(z)=P0(z)-P3(z)
      YP(z)=0.5d0*log(PP(z)/PM(z))
      W0(z)=0.5d0*(1d0+om**2-z)
      WP(z)=W0(z)+P3(z)
      WM(z)=W0(z)-P3(z)
      YW(z)=0.5d0*log(WP(z)/WM(z))
c     End statement functions.

      omsq=s12/mt**2
      om=sqrt(omsq)
      be=mb/mt
      besq=be**2
C      zm=(1d0-om)**2


      if (scheme .eq.'dred') then
C------        eta=0 4d-hel
         eta=0d0
      elseif (scheme .eq. 'tH-V') then
C------       eta=1 t'Hooft Veltman
         eta=1d0
      endif

      if ((mb .eq. 0d0) .and. (mt .eq. 0d0)) then 
           s12log=lnrat(musq,-s12)
C----Again these are corrections to vertex in units of as/4/pi*CF
           C0L=-2d0*(epinv**2+epinv*s12log+0.5d0*s12log**2)
     &     -3d0*(epinv+s12log)-7d0-eta
           C0R=czip
           C1L=czip
           C1R=czip
           return
      endif

      rb=mb/mt
      rbsq=rb**2
      mbsq=mb**2
      mtsq=mt**2
      mulog=log(musq/mtsq)
      C0(-2)=qlI3(mbsq,s12,mtsq,0d0,mbsq,mtsq,musq,-2)
      C0(-1)=qlI3(mbsq,s12,mtsq,0d0,mbsq,mtsq,musq,-1)
      C0( 0)=qlI3(mbsq,s12,mtsq,0d0,mbsq,mtsq,musq, 0)
      Zt=-0.5d0*(3d0*(epinv+mulog)+5d0-eta)

      if (abs(mb) .lt. 1d-6) then 
      Zb=0d0
      lnrbsq=0d0
      B0f1=-qlI2(mtsq,0d0,mtsq,musq,-1)*epinv-qlI2(mtsq,0d0,mtsq,musq,0)
      B0f2=qlI2(s12,mbsq,mtsq,musq,0)-qlI2(mtsq,0d0,mtsq,musq,0)
      else
      lnrbsq=log(rbsq)
      Zb=-0.5d0*(3d0*(epinv+mulog-lnrbsq)+5d0-eta)
      B0f1=qlI2(mbsq,0d0,mbsq,musq,0)-qlI2(mtsq,0d0,mtsq,musq,0)
      B0f2=qlI2(s12,mbsq,mtsq,musq,0)-qlI2(mtsq,0d0,mtsq,musq,0)
      endif

      rwsq=s12/mtsq
      lambda=(1d0-rwsq-rbsq)**2-4d0*rbsq*rwsq
      
      if (includect) then
          if (abs(mb) .lt. 1d-6) then
c---          normal massless result
              rsq=s12/mtsq
              omrsq=1d0-rsq
              wlog=dlog(omrsq)
              rlog=dlog(rsq)
c---   epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967

               Ct=(epinv2*epinv+epinv*mulog+0.5d0*mulog**2)
     &          +(epinv+mulog)*(2.5d0-2d0*wlog)
     &          +25d0/4d0+0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog
     &          +0.5d0/omrsq+2d0*ddilog(omrsq)-5d0*pisqo6
     &          -5d0*wlog+2d0*wlog**2+eta/2d0
     &          -2d0*dlog(aff)**2-(3.5d0-4d0*aff+aff**2/2d0)*dlog(aff)
     &          +2d0*(1d0-aff)*rsq/omrsq*dlog(rsq/(1d0-aff+rsq*aff))
          else

c--- massive result: counter-term from cleaned-up Czarnecki, cf. Eqs. (25)
c--- (modified to absorb an overall factor of log(2*P3b))
c--- [in units of 2*gsq*Cg*CF -> ason2pi*Cf]
c      Cz1=2d0*(-0.5d0*(epinv+mulog)-1d0+dlog(4d0*P3b**2/(be*om))
c     & -(1d0-omsq)*Ypb/P3b-0.5d0*(1d0-omsq-besq)*Ywb/P3b)
c      Cz2=-2d0*P0b/P3b*(
c     &         2d0*Ypb*(-0.5d0*(epinv+mulog)+dlog(4d0*P3b**2*(zm-besq)/be))
c     & +ddilog(1d0-Pmb)-ddilog(1d0-Ppb)-3d0*ddilog(1d0-Pmb/Ppb)
c     & -3d0*Ypb**2+2d0*dlog(Pmb)*dlog(1d0-om-Pmb)
c     & -2d0*dlog(Ppb)*dlog(Ppb-1d0+om))
c      Cz3=2d0*(-0.5d0*(epinv+mulog)-1d0+dlog(4d0*P3b**2/(be*om))
c     & -P0b*Ypb/P3b-W0b*Ywb/P3b)     
c--- note incude a minus sign for integrated subtraction    
c      Ct=-(Cz1+Cz3+Cz2)
c--- massive result: formula from paper [in units of 2*gsq*Cg*CF -> ason2pi*Cf]

               P0b=P0(besq)
               P3b=P3(besq)
               PMb=PM(besq)
               PPb=PP(besq)
               Ypb=YP(besq)
               Ywb=YW(besq)

               Ct=(2d0*(epinv+mulog)-4d0*dlog(4d0*P3b**2/om/be))
     &          *(1d0-P0b/P3b*Ypb)
     &          +4d0+2d0/P3b*((1d0-omsq)*Ypb+(1d0-besq)*Ywb)
     &          +P0b/P3b*(2d0*Ypb-6d0*Ypb**2+4d0*Ywb*dlog(be)
     &          -6d0*ddilog(1d0-Pmb/Ppb)-2d0*ddilog(1d0-Ppb)
     &         +2d0*ddilog(1d0-Pmb))
          endif
      else
          Ct=0d0
      endif

      C0L=Zt+Zb+epinv+mulog+2d0-eta
     & +2d0*(1d0+rbsq-rwsq)*mt**2*(C0(-2)*epinv**2+C0(-1)*epinv+C0(0))
     & +(2d0*rwsq*(1d0+rbsq-rwsq/2d0)-(1d0-rbsq)**2
     & +2d0*B0f1*(1d0-rbsq-rwsq*(2d0+rbsq-rwsq))
     & -B0f2*((1d0-rbsq)**2-4d0*rwsq*(1d0+rbsq-0.75d0*rwsq)))/lambda
      C0R=2d0*rb/lambda*(B0f1*(1d0-rbsq-rwsq)+2d0*B0f2*rwsq)
      C1L=-2d0*rb/lambda*((lnrbsq+2d0*B0f1)*(1d0+rbsq-rwsq)
     &     +B0f2*(1d0-rbsq+rwsq))
      C1R=2d0/lambda*(2d0*lnrbsq*rbsq
     & +4d0*B0f1*rbsq+B0f2*(1d0-rbsq-rwsq))

      return
      end

