      subroutine qg_tbqdk_v(p,msq)
************************************************************************
*     Virtual t-channel single top, with explicit b-quark              *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *      
*                                                                      *
*     Originally: R. Frederix and F. Tramontano, February 2008         *
*        Adapted: J. Campbell, February 27, 2008                       *
*    added decay: J. Campbell, May 2011                                *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'scheme.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'stopscales.f'
      include 'sck.f'
      double precision p(mxpart,4),fac
c--- needed for pole check
c     . ,xs,xsn,xsd
      double precision msq(-nf:nf,-nf:nf),msq_qg,msq_gq,msq_qbarg,
     . msq_gqbar,dot,Wprop15,Wprop25,xsqV,xsqR,mq,ma,gsq_H
      double complex 
     . LOamps_qg(2,2),Virtamps_qg(2,2),
     . LOamps_qbarg(2,2),Virtamps_qbarg(2,2),
     . LOamps_gq(2,2),Virtamps_gq(2,2),
     . LOamps_gqbar(2,2),Virtamps_gqbar(2,2)
c--- needed for pole check
c     . ,lp,lRc1,lRs1
      integer hg,hs,j,k,i3,i4

      sck=1d0

      scheme='dred'

c---initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--- DEBUG: to check alpha-dependence      
c      return

c--- set mass of quark and antiquark according to nwz
      if (nwz .eq. +1) then
        mq=mt
        ma=mb
        i3=3
        i4=4
      else
        mq=mb
        ma=mt
        i3=4
        i4=3
      endif

c--- variables to pass renormalization scale in to virtual routines
      xsqV=renscale_H**2
      xsqR=renscale_H**2

c--- note: factor of ason2pi moved inside virtwrap routine compared
c---       to earlier versions, to allow for different scales
      gsq_H=fourpi*as_H
      fac=aveqg*2d0*xn**2*Cf*gsq_H*gwsq**2

c--- propagator for qg and qbarg amplitudes
      Wprop15=1d0/(2d0*dot(p,1,7)-wmass**2)

      call virtwrapdk(p,1,2,i3,i4,7,
     .                mq,ma,xsqV,xsqR,LOamps_qg,Virtamps_qg)
      call virtwrapdk(p,7,2,i3,i4,1,
     .                mq,ma,xsqV,xsqR,LOamps_qbarg,Virtamps_qbarg)

c--- propagator for gq and gqbar amplitudes
      Wprop25=1d0/(2d0*dot(p,2,7)-wmass**2)

      call virtwrapdk(p,2,1,i3,i4,7,
     .              mq,ma,xsqV,xsqR,LOamps_gq,Virtamps_gq)
      call virtwrapdk(p,7,1,i3,i4,2,
     .              mq,ma,xsqV,xsqR,LOamps_gqbar,Virtamps_gqbar)

      msq_qg=0d0
      msq_qbarg=0d0
      msq_gq=0d0
      msq_gqbar=0d0
      do hg=1,2
      do hs=1,2
c      write(6,*) hg,hc,hs,LOamps_gq(hg,hc,hs),Virtamps_gqbar(hg,hc,hs)
      msq_qg=msq_qg+Wprop15**2*dble(
     .       Virtamps_qg(hg,hs)*dconjg(LOamps_qg(hg,hs)))
      msq_qbarg=msq_qbarg+Wprop15**2*dble(
     .       Virtamps_qbarg(hg,hs)*dconjg(LOamps_qbarg(hg,hs)))
      msq_gq=msq_gq+Wprop25**2*dble(
     .       Virtamps_gq(hg,hs)*dconjg(LOamps_gq(hg,hs)))
      msq_gqbar=msq_gqbar+Wprop25**2*dble(
     .       Virtamps_gqbar(hg,hs)*dconjg(LOamps_gqbar(hg,hs)))
      enddo
      enddo

c--- fill matrix elements
      do j=1,4
         msq(+j,0)=fac*Vsum(+j)*msq_qg
         msq(-j,0)=fac*Vsum(-j)*msq_qbarg
         msq(0,+j)=fac*Vsum(+j)*msq_gq
         msq(0,-j)=fac*Vsum(-j)*msq_gqbar
      enddo

      return
      end


c--- wrapper to virtual and LO amplitude routines that allows the
c---  momenta to be permuted according to i1,i2,i5
      subroutine virtwrapdk(p,i1,i2,i3,i4,i5,
     .                    mh,ml,xsqV,xsqR,LOampsdk,Virtampsdk)
      implicit none
      include 'constants.f'
      include 'b0.f'
      include 'epinv.f'
      include 'scale.f'
      include 'zprods_com.f'
      include 'stopscales.f'
      include 'sck.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'nwz.f'
      integer i1,i2,i3,i4,i5,j,hg,hs
      double precision p(mxpart,4),q(mxpart,4),dot,mh,ml,xsqV,xsqR,
     . colA,colB,virt_massless,eta,ren,ason2pi_H,ason2pi_L,s345,fac
      double complex LOamps(2,2,2),Virtamps(2,2,2),
     .  LOampsdk(2,2),Virtampsdk(2,2),
     . Appp,Appm,Apmp,Apmm,Ampp,Ampm,Ammp,Ammm,
     . Bppp,Bppm,Bpmp,Bpmm,Bmpp,Bmpm,Bmmp,Bmmm,
     . xl15,lnrat
      sck=1d0

c--- factors of ason2pi now included in this routine
      ason2pi_H=as_H/twopi
      ason2pi_L=as_L/twopi
            
      do j=1,4
        q(1,j)=p(i1,j)
        q(2,j)=p(i2,j) ! eta
        q(i3,j)=p(3,j)+p(4,j)+p(5,j)
     &     -mt**2/2d0/(dot(p,i2,3)+dot(p,i2,4)+dot(p,i2,5))*p(i2,j) ! t1
        q(i4,j)=p(6,j)-mb**2/2d0/dot(p,i2,6)*p(i2,j)
        q(5,j)=p(i5,j)
        q(6,j)=p(i4,j) ! e
      enddo
      
c--- set up spinor products
      call spinoru(6,q,za,zb)

      colA=ca/2d0
      colB=(2d0*cf-ca)/2d0

      eta=0d0  !dred scheme
      ren=(-b0*epinv
     .     -cf*(3d0/2d0*epinv + (4d0+1d0-eta)/2d0 + 3d0*log(scale/mh))
     .     -cf*(1d0/2d0*epinv +
     .             sck*(epinv + (4d0+1d0-eta)/2d0 + 3d0*log(scale/ml)))) 

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      ren=ren+xn/6d0

c--- include finite renormalization per FT's message of 20/12/2008
c---  -2*cf*LeadingOrder*alphas/2*pi
      ren=ren-2d0*cf

c--- go from OS scheme for a massive quark to MSbar for massless one
c--- (strong-coupling correction)
      ren=ren-2d0/3d0*log(renscale_H/ml)

c--- go from OS scheme for a massive quark to MSbar for massless one
c--- (gluon PDF correction)
      ren=ren-2d0/3d0*log(ml/facscale_H)

c--- filling common block
      call stop_def(xsqV,xsqR,q,mh,ml)
c---- calling amps(hg,hc,hs)
      call bornampsN(q,mh,ml,LOamps)
      call Aamp_ppp(q,mh,ml,Appp)
      call Aamp_ppm(q,mh,ml,Appm)
      call Aamp_pmp(q,mh,ml,Apmp)
      call Aamp_pmm(q,mh,ml,Apmm)
      call Aamp_mpp(q,mh,ml,Ampp)
      call Aamp_mpm(q,mh,ml,Ampm)
      call Aamp_mmp(q,mh,ml,Ammp)
      call Aamp_mmm(q,mh,ml,Ammm)
      call Bamp_ppp(q,mh,ml,Bppp)
      call Bamp_ppm(q,mh,ml,Bppm)
      call Bamp_pmp(q,mh,ml,Bpmp)
      call Bamp_pmm(q,mh,ml,Bpmm)
      call Bamp_mpp(q,mh,ml,Bmpp)
      call Bamp_mpm(q,mh,ml,Bmpm)
      call Bamp_mmp(q,mh,ml,Bmmp)
      call Bamp_mmm(q,mh,ml,Bmmm)
 
      xl15=lnrat(-2d0*dot(p,i1,i5),renscale_L**2)
c--- correction to the massless line (cf. cv0 in qqb_tbb_v.f)
      virt_massless=-2d0*epinv*(epinv-dble(xl15))-dble(xl15**2)
     .              -3d0*(epinv-dble(xl15))-7d0
      virt_massless=virt_massless*cf

c--- apply factors of ason2pi now
      colA=colA*ason2pi_H
      colB=colB*ason2pi_H    ! corrections on heavy line
      virt_massless=virt_massless*ason2pi_L ! on light line
      ren=ren*ason2pi_H ! renormalization is like LO
 
      Virtamps(2,2,2)=Appp*colA+Bppp*colB
     .               +(ren+virt_massless)*LOamps(2,2,2)
      Virtamps(2,2,1)=Appm*colA+Bppm*colB
     .               +(ren+virt_massless)*LOamps(2,2,1)
      Virtamps(2,1,2)=Apmp*colA+Bpmp*colB
     .               +(ren+virt_massless)*LOamps(2,1,2)
      Virtamps(2,1,1)=Apmm*colA+Bpmm*colB
     .               +(ren+virt_massless)*LOamps(2,1,1)
      Virtamps(1,2,2)=Ampp*colA+Bmpp*colB
     .               +(ren+virt_massless)*LOamps(1,2,2)
      Virtamps(1,2,1)=Ampm*colA+Bmpm*colB
     .               +(ren+virt_massless)*LOamps(1,2,1)
      Virtamps(1,1,2)=Ammp*colA+Bmmp*colB
     .               +(ren+virt_massless)*LOamps(1,1,2)
      Virtamps(1,1,1)=Ammm*colA+Bmmm*colB
     .               +(ren+virt_massless)*LOamps(1,1,1)

      s345=(p(3,4)+p(4,4)+p(5,4))**2
     &    -(p(3,1)+p(4,1)+p(5,1))**2
     &    -(p(3,2)+p(4,2)+p(5,2))**2
     &    -(p(3,3)+p(4,3)+p(5,3))**2
     
c--- now dress up with appropriate factors to include the top quark decay
      fac=gwsq*dsqrt(2d0*dot(p,i3,5))
     &    /dsqrt((2d0*dot(p,3,4)-wmass**2)**2+(wmass*wwidth)**2)
     &    /dsqrt((s345-mt**2)**2+(mt*twidth)**2)
     
      do hg=1,2
      do hs=1,2
      if     (nwz .eq. +1) then
        LOampsdk(hg,hs)=fac*(
     &   +LOamps(hg,1,hs)*zb(6,i3)
     &   +LOamps(hg,2,hs)*mt*zb(6,2)/zb(i3,2))
        Virtampsdk(hg,hs)=fac*(
     &   +Virtamps(hg,1,hs)*zb(6,i3)
     &   +Virtamps(hg,2,hs)*mt*zb(6,2)/zb(i3,2))
      elseif (nwz .eq. -1) then
        LOampsdk(hg,hs)=fac*(
     &   +LOamps(hg,hs,2)*za(6,i3)
     &   +LOamps(hg,hs,1)*mt*za(6,2)/za(i3,2))
        Virtampsdk(hg,hs)=fac*(
     &   +Virtamps(hg,hs,2)*za(6,i3)
     &   +Virtamps(hg,hs,1)*mt*za(6,2)/za(i3,2))
      else
        write(6,*) 'nwz must be +1 or -1 in virtwrapdk'
      stop
      endif
      enddo
      enddo
      
      return



************************************************************************
*   CODE BELOW HERE IS FOR CHECKING POLES ONLY                         *
************************************************************************


c      call checkint(p,xsqV,xsqR,ml**2,mh**2)

c--- needed for pole check

c      xsn=(1d0-dsqrt(1d0-4d0*ml*mh/(2*dot(p,3,4)+2*ml*mh)))
c      xsd=(1d0+dsqrt(1d0-4d0*ml*mh/(2*dot(p,3,4)+2*ml*mh)))
c      xs=-xsn/xsd
c--- these are to check the qg terms
c      lRc1=lnrat(mh*sqrt(xsqR),-2*dot(p,2,3))
c      lRs1=lnrat(ml*sqrt(xsqR),-2*dot(p,2,4))
c      lp=lnrat(xsn,-xsd)

c--- these are to check the gq terms
c      lRc1=lnrat(mh*sqrt(xsqR),-2*dot(p,1,3))
c      lRs1=lnrat(ml*sqrt(xsqR),-2*dot(p,1,4))
c--- pole check
c
c      write(6,*) 'Appp/LOamp(2,2,2)=',Appp/(+LOamps(2,2,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Appm/LOamp(2,2,1)=',Appm/(+LOamps(2,2,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Apmp/LOamp(2,1,2)=',Apmp/(+LOamps(2,1,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Apmm/LOamp(2,1,1)=',Apmm/(+LOamps(2,1,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ampp/LOamp(1,2,2)=',Ampp/(+LOamps(1,2,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ampm/LOamp(1,2,1)=',Ampm/(+LOamps(1,2,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ammp/LOamp(1,1,2)=',Ammp/(+LOamps(1,1,2)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))
c      write(6,*) 'Ammm/LOamp(1,1,1)=',Ammm/(+LOamps(1,1,1)
c     . *(-2d0*epinv**2 - epinv*(-1d0 + 2d0*lRc1 + 2d0*lRs1)))

c      write(6,*) 'Bppp/LOamp(2,2,2)=',Bppp/(+LOamps(2,2,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bppm/LOamp(2,2,1)=',Bppm/(+LOamps(2,2,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bpmp/LOamp(2,1,2)=',Bpmp/(+LOamps(2,1,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bpmm/LOamp(2,1,1)=',Bpmm/(+LOamps(2,1,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmpp/LOamp(1,2,2)=',Bmpp/(+LOamps(1,2,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmpm/LOamp(1,2,1)=',Bmpm/(+LOamps(1,2,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmmp/LOamp(1,1,2)=',Bmmp/(+LOamps(1,1,2)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)
c      write(6,*) 'Bmmm/LOamp(1,1,1)=',Bmmm/(+LOamps(1,1,1)
c     .    *epinv*(-4d0*dot(p,3,4)*lp*xs + mh*ml*(-1d0 + xs**2))/
c     .  (ml*(-1d0 + xs**2))/mh)

c      pause
c      return
      
      end
      
      
