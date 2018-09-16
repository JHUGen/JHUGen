      subroutine coefsdkmass(s12,mt,mb,ct,c0L,c0R,c1L,C1R)
      implicit none
      include 'types.f'
      
C     Authors: John Campbell and R.Keith Ellis, April 2012
C     This gives the coefficients of the 
C     virtual corrections, C0L,C0R,C1L,C1R, 
C     and if includect is set to .true. 
C     the integrated counterterm for top semi-leptonic decay, ct
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'includect.f'
      real(dp)::s12,mb,mt,ct,
     & rsq,omrsq,eta,wlog,rlog,mulog,lambda,
     & rb,rbsq,rwsq,Zt,Zb,mtsq,mbsq,lnrbsq,ddilog
      complex(dp)::B0f1,B0f2,C0(-2:0),c0L,c0R,c1L,C1R,s12log,lnrat,
     & qlI3,qlI2
      real(dp)::P0,P3,PP,PM,W0,WP,WM,YP,YW,z,om,omsq
      real(dp)::P0b,P3b,PMb,PPb,Ypb,Ywb,be,besq
c     Statement functions.
      P0(z)=0.5_dp*(1._dp-om**2+z)
      P3(z)=0.5_dp*sqrt(1._dp+om**4+z**2-2._dp*(om**2+z+om**2*z))
      PP(z)=P0(z)+P3(z)
      PM(z)=P0(z)-P3(z)
      YP(z)=0.5_dp*log(PP(z)/PM(z))
      W0(z)=0.5_dp*(1._dp+om**2-z)
      WP(z)=W0(z)+P3(z)
      WM(z)=W0(z)-P3(z)
      YW(z)=0.5_dp*log(WP(z)/WM(z))
c     End statement functions.

      omsq=s12/mt**2
      om=sqrt(omsq)
      be=mb/mt
      besq=be**2
C      zm=(1._dp-om)**2


      if (scheme =='dred') then
C------        eta=0 4.e-_dphel
         eta=0._dp
      elseif (scheme == 'tH-V') then
C------       eta=1 t'Hooft Veltman
         eta=1._dp
      endif

      if ((mb == 0._dp) .and. (mt == 0._dp)) then 
           s12log=lnrat(musq,-s12)
C----Again these are corrections to vertex in units of as/4/pi*CF
           C0L=-2._dp*(epinv**2+epinv*s12log+0.5_dp*s12log**2)
     &     -3._dp*(epinv+s12log)-7._dp-eta
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
      C0(-2)=qlI3(mbsq,s12,mtsq,0._dp,mbsq,mtsq,musq,-2)
      C0(-1)=qlI3(mbsq,s12,mtsq,0._dp,mbsq,mtsq,musq,-1)
      C0( 0)=qlI3(mbsq,s12,mtsq,0._dp,mbsq,mtsq,musq, 0)
      Zt=-0.5_dp*(3._dp*(epinv+mulog)+5._dp-eta)

      if (abs(mb) < 1.e-6_dp) then 
      Zb=0._dp
      lnrbsq=0._dp
      B0f1=-qlI2(mtsq,0._dp,mtsq,musq,-1)*epinv-qlI2(mtsq,0._dp,mtsq,musq,0)
      B0f2=qlI2(s12,mbsq,mtsq,musq,0)-qlI2(mtsq,0._dp,mtsq,musq,0)
      else
      lnrbsq=log(rbsq)
      Zb=-0.5_dp*(3._dp*(epinv+mulog-lnrbsq)+5._dp-eta)
      B0f1=qlI2(mbsq,0._dp,mbsq,musq,0)-qlI2(mtsq,0._dp,mtsq,musq,0)
      B0f2=qlI2(s12,mbsq,mtsq,musq,0)-qlI2(mtsq,0._dp,mtsq,musq,0)
      endif

      rwsq=s12/mtsq
      lambda=(1._dp-rwsq-rbsq)**2-4._dp*rbsq*rwsq
      
      if (includect) then
          if (abs(mb) < 1.e-6_dp) then
c---          normal massless result
              rsq=s12/mtsq
              omrsq=1._dp-rsq
              wlog=log(omrsq)
              rlog=log(rsq)
c---   epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967

               Ct=(epinv2*epinv+epinv*mulog+0.5_dp*mulog**2)
     &          +(epinv+mulog)*(2.5_dp-2._dp*wlog)
     &          +25._dp/4._dp+0.5_dp*(1._dp/omrsq**2-8._dp/omrsq+7._dp)*rlog
     &          +0.5_dp/omrsq+2._dp*ddilog(omrsq)-5._dp*pisqo6
     &          -5._dp*wlog+2._dp*wlog**2+eta/2._dp
     &          -2._dp*log(aff)**2-(3.5_dp-4._dp*aff+aff**2/2._dp)*log(aff)
     &          +2._dp*(1._dp-aff)*rsq/omrsq*log(rsq/(1._dp-aff+rsq*aff))
          else

c--- massive result: counter-term from cleane.e-_dpup Czarnecki, cf. Eqs. (25)
c--- (modified to absorb an overall factor of log(2*P3b))
c--- [in units of 2*gsq*Cg*CF -> ason2pi*Cf]
c      Cz1=2._dp*(-0.5_dp*(epinv+mulog)-1._dp+log(4._dp*P3b**2/(be*om))
c     & -(1._dp-omsq)*Ypb/P3b-0.5_dp*(1._dp-omsq-besq)*Ywb/P3b)
c      Cz2=-2._dp*P0b/P3b*(
c     &         2._dp*Ypb*(-0.5_dp*(epinv+mulog)+log(4._dp*P3b**2*(zm-besq)/be))
c     & +ddilog(1._dp-Pmb)-ddilog(1._dp-Ppb)-3._dp*ddilog(1._dp-Pmb/Ppb)
c     & -3._dp*Ypb**2+2._dp*log(Pmb)*log(1._dp-om-Pmb)
c     & -2._dp*log(Ppb)*log(Ppb-1._dp+om))
c      Cz3=2._dp*(-0.5_dp*(epinv+mulog)-1._dp+log(4._dp*P3b**2/(be*om))
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

               Ct=(2._dp*(epinv+mulog)-4._dp*log(4._dp*P3b**2/om/be))
     &          *(1._dp-P0b/P3b*Ypb)
     &          +4._dp+2._dp/P3b*((1._dp-omsq)*Ypb+(1._dp-besq)*Ywb)
     &          +P0b/P3b*(2._dp*Ypb-6._dp*Ypb**2+4._dp*Ywb*log(be)
     &          -6._dp*ddilog(1._dp-Pmb/Ppb)-2._dp*ddilog(1._dp-Ppb)
     &         +2._dp*ddilog(1._dp-Pmb))
          endif
      else
          Ct=0._dp
      endif

      C0L=Zt+Zb+epinv+mulog+2._dp-eta
     & +2._dp*(1._dp+rbsq-rwsq)*mt**2*(C0(-2)*epinv**2+C0(-1)*epinv+C0(0))
     & +(2._dp*rwsq*(1._dp+rbsq-rwsq/2._dp)-(1._dp-rbsq)**2
     & +2._dp*B0f1*(1._dp-rbsq-rwsq*(2._dp+rbsq-rwsq))
     & -B0f2*((1._dp-rbsq)**2-4._dp*rwsq*(1._dp+rbsq-0.75_dp*rwsq)))/lambda
      C0R=2._dp*rb/lambda*(B0f1*(1._dp-rbsq-rwsq)+2._dp*B0f2*rwsq)
      C1L=-2._dp*rb/lambda*((lnrbsq+2._dp*B0f1)*(1._dp+rbsq-rwsq)
     &     +B0f2*(1._dp-rbsq+rwsq))
      C1R=2._dp/lambda*(2._dp*lnrbsq*rbsq
     & +4._dp*B0f1*rbsq+B0f2*(1._dp-rbsq-rwsq))

      return
      end

