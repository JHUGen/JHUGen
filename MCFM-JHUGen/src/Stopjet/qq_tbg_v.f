      subroutine qq_tbg_v(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Virtual s-channel single top + jet                               *
*                                                                      *
*     q(p1) + q(p2) -> t(p3) + b(p4) + g(p5)                           *
*                                                                      *
*     Originally: R. Frederix and F. Tramontano, February 2008         *
*        Adapted: J. Campbell, June 20, 2008                           *
*                                                                      *
************************************************************************
*                                                                      *
*     IMPORTANT NOTE!                                                  *
*                                                                      *
*     For now, we only include radiation from heavy quark line         *
*      and virtual corrections on that line too                        *
*                                                                      *
************************************************************************
c     u + g  ->  c + s + d  (t-channel single-charm)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'scheme.f'
      include 'nwz.f'
      include 'stopscales.f'
      include 'sck.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     & dot,Wprop12,xsqV,xsqR,mq,ma,gsq_H,
     & msq_qqbar_lc,msq_qqbar_slc,msq_qqbar_tr
c--- needed for pole check
c     & ,xs,xsn,xsd
      complex(dp)::
     & LOamps_qqbar(2,2,2),
     & Virtamps_qqbar_lc(2,2,2),Virtamps_qqbar_slc(2,2,2),
     & Virtamps_qqbar_tr(2,2,2)
c--- needed for pole check
c     & ,lp,lRc1,lRs1
      integer:: hg,hc,hs,j,k,i3,i4

      sck=1._dp

      scheme='dred'

c---initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--- DEBUG: to check alpha-dependence
c      return

c--- set mass of quark and antiquark according to nwz
      if (nwz == +1) then
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
      fac=aveqq*2._dp*xn**2*Cf*gsq_H*gwsq**2

c--- propagator for qg and qbarg amplitudes
c--- note: in principle we could need a width here, if the top quark
c---       is far off-shell; for simplicity, we do not include it.
      Wprop12=1._dp/(2._dp*dot(p,1,2)-wmass**2)

      call virtwrap_heavyonly(
     &      p,1,5,i3,i4,2,mq,ma,xsqV,xsqR,LOamps_qqbar,
     &      Virtamps_qqbar_lc,Virtamps_qqbar_slc,Virtamps_qqbar_tr)
c--- Note: qbq matrix elements are later excluded by PDF choice anyway
c      call virtwrap_heavyonly(p,2,5,i3,i4,1,
c     &              mq,ma,xsqV,xsqR,LOamps_qbarq,Virtamps_qbarq)

      msq_qqbar_lc=0._dp
      msq_qqbar_slc=0._dp
      msq_qqbar_tr=0._dp
c      msq_qbarq=0._dp
      do hg=1,2
      do hc=1,2
      do hs=1,2
      msq_qqbar_lc=msq_qqbar_lc+Wprop12**2*real(
     &     Virtamps_qqbar_lc(hg,hc,hs)*conjg(LOamps_qqbar(hg,hc,hs)))
      msq_qqbar_slc=msq_qqbar_slc+Wprop12**2*real(
     &     Virtamps_qqbar_slc(hg,hc,hs)*conjg(LOamps_qqbar(hg,hc,hs)))
      msq_qqbar_tr=msq_qqbar_tr+Wprop12**2*real(
     &     Virtamps_qqbar_tr(hg,hc,hs)*conjg(LOamps_qqbar(hg,hc,hs)))
c      msq_qbarq=msq_qbarq+Wprop12**2*real(
c     &       Virtamps_qbarq(hg,hc,hs)*conjg(LOamps_qbarq(hg,hc,hs)))
      enddo
      enddo
      enddo

c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)
      msq(0,0)=fac*msq_qqbar_lc
      msq(0,1)=fac*msq_qqbar_slc
      msq(1,0)=fac*msq_qqbar_tr

cc--- fill matrix elements
c      do j=-4,4
c      do k=-4,4
c        if     ((j > 0) .and. (k < 0)) then
c        msq(j,k)=fac*msq_qqbar
cc        elseif ((j < 0) .and. (k > 0)) then
cc        msq(j,k)=fac*Vsq(j,k)*msq_qbarq
c      endif
c      enddo
c      enddo

      return
      end


c--- wrapper to virtual and LO amplitude routines that allows the
c---  momenta to be permuted according to i1,i2,i5
c--- NOTE: this is an exact copy of "virtwrap", but the corrections
c---       to the massless line are set to zero
      subroutine virtwrap_heavyonly(p,i1,i2,i3,i4,i5,
     &                    mh,ml,xsqV,xsqR,LOamps,
     &                    Virtamps_lc,Virtamps_slc,Virtamps_tr)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'nflav.f'
      include 'scale.f'
      include 'zprods_com.f'
      include 'stopscales.f'
      include 'sck.f'
      integer:: i1,i2,i3,i4,i5,j
      real(dp):: p(mxpart,4),q(mxpart,4),dot,mh,ml,xsqV,xsqR,
     & colA,colB,ren_lc,ren_slc,ren_tr,eta,
     & ason2pi_H
      complex(dp):: LOamps(2,2,2),
     & Virtamps_lc(2,2,2),Virtamps_slc(2,2,2),Virtamps_tr(2,2,2),
     & Appp,Appm,Apmp,Apmm,Ampp,Ampm,Ammp,Ammm,
     & Bppp,Bppm,Bpmp,Bpmm,Bmpp,Bmpm,Bmmp,Bmmm
      sck=1._dp

c--- factors of ason2pi now included in this routine
      ason2pi_H=as_H/twopi
c      ason2pi_L=as_L/twopi

      do j=1,4
        q(1,j)=p(i1,j)
        q(2,j)=p(i2,j)
        q(3,j)=p(i3,j)-mh**2/2._dp/dot(p,i2,i3)*p(i2,j)
        q(4,j)=p(i4,j)-ml**2/2._dp/dot(p,i2,i4)*p(i2,j)
        q(5,j)=p(i5,j)
      enddo

c--- set up spinor products
      call spinoru(5,q,za,zb)

      colA=ca/2._dp
      colB=(2._dp*cf-ca)/2._dp
      eta=0._dp  !dred scheme
c      ren=(-b0*epinv
c     &     -cf*(3._dp/2._dp*epinv + (4._dp+1._dp-eta)/2._dp + 3._dp*log(scale/mh))
c     &     -cf*(1._dp/2._dp*epinv +
c     &             sck*(epinv + (4._dp+1._dp-eta)/2._dp + 3._dp*log(scale/ml))))
      ren_tr=+2._dp*tr*real(nflav,dp)/3._dp*epinv

      ren_lc=(-11._dp/6._dp*ca*epinv
     &     -ca/2._dp*(3._dp/2._dp*epinv
     &              + (4._dp+1._dp-eta)/2._dp + 3._dp*log(scale/mh))
     &     -ca/2._dp*(1._dp/2._dp*epinv +
     &             sck*(epinv + (4._dp+1._dp-eta)/2._dp + 3._dp*log(scale/ml))))

      ren_slc=(
     &     -(2._dp*cf-ca)/2._dp*(3._dp/2._dp*epinv
     &              + (4._dp+1._dp-eta)/2._dp + 3._dp*log(scale/mh))
     &     -(2._dp*cf-ca)/2._dp*(1._dp/2._dp*epinv +
     &             sck*(epinv + (4._dp+1._dp-eta)/2._dp + 3._dp*log(scale/ml))))

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      ren_lc=ren_lc+ca/6._dp

c--- include finite renormalization per FT's message of 20/12/2008
c---  -2*cf*LeadingOrder*alphas/2*pi
      ren_lc=ren_lc-ca
      ren_slc=ren_slc-(2._dp*cf-ca)

c--- go from OS scheme for a massive quark to MSbar for massless one
c--- (strong-coupling correction)
c      ren_tr=ren_tr-tr*4._dp/3._dp*log(scale/mh)

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

c      xl15=lnrat(-2._dp*dot(p,i1,i5),renscale_L**2)
c--- correction to the massless line (cf. cv0 in qqb_tbb_v.f)
c      virt_massless=-2._dp*epinv*(epinv-real(xl15))-real(xl15**2)
c     &              -3._dp*(epinv-real(xl15))-7._dp
c      virt_massless=virt_massless*cf
c      virt_massless=0._dp

c--- apply factors of ason2pi now
      colA=colA*ason2pi_H
      colB=colB*ason2pi_H    ! corrections on heavy line
c      virt_massless=virt_massless*ason2pi_L ! on light line
      ren_lc=ren_lc*ason2pi_H ! renormalization is like LO
      ren_slc=ren_slc*ason2pi_H ! renormalization is like LO
      ren_tr=ren_tr*ason2pi_H ! renormalization is like LO

      Virtamps_lc(2,2,2)=Appp*colA+ren_lc*LOamps(2,2,2)
      Virtamps_lc(2,2,1)=Appm*colA+ren_lc*LOamps(2,2,1)
      Virtamps_lc(2,1,2)=Apmp*colA+ren_lc*LOamps(2,1,2)
      Virtamps_lc(2,1,1)=Apmm*colA+ren_lc*LOamps(2,1,1)
      Virtamps_lc(1,2,2)=Ampp*colA+ren_lc*LOamps(1,2,2)
      Virtamps_lc(1,2,1)=Ampm*colA+ren_lc*LOamps(1,2,1)
      Virtamps_lc(1,1,2)=Ammp*colA+ren_lc*LOamps(1,1,2)
      Virtamps_lc(1,1,1)=Ammm*colA+ren_lc*LOamps(1,1,1)

      Virtamps_slc(2,2,2)=Bppp*colB+ren_slc*LOamps(2,2,2)
      Virtamps_slc(2,2,1)=Bppm*colB+ren_slc*LOamps(2,2,1)
      Virtamps_slc(2,1,2)=Bpmp*colB+ren_slc*LOamps(2,1,2)
      Virtamps_slc(2,1,1)=Bpmm*colB+ren_slc*LOamps(2,1,1)
      Virtamps_slc(1,2,2)=Bmpp*colB+ren_slc*LOamps(1,2,2)
      Virtamps_slc(1,2,1)=Bmpm*colB+ren_slc*LOamps(1,2,1)
      Virtamps_slc(1,1,2)=Bmmp*colB+ren_slc*LOamps(1,1,2)
      Virtamps_slc(1,1,1)=Bmmm*colB+ren_slc*LOamps(1,1,1)

      Virtamps_tr(2,2,2)=ren_tr*LOamps(2,2,2)
      Virtamps_tr(2,2,1)=ren_tr*LOamps(2,2,1)
      Virtamps_tr(2,1,2)=ren_tr*LOamps(2,1,2)
      Virtamps_tr(2,1,1)=ren_tr*LOamps(2,1,1)
      Virtamps_tr(1,2,2)=ren_tr*LOamps(1,2,2)
      Virtamps_tr(1,2,1)=ren_tr*LOamps(1,2,1)
      Virtamps_tr(1,1,2)=ren_tr*LOamps(1,1,2)
      Virtamps_tr(1,1,1)=ren_tr*LOamps(1,1,1)

      return
      end
