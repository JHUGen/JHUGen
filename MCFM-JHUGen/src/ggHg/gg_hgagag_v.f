      subroutine gg_hgagag_v(p,msq)
C----Author: R.K. Ellis August 2011
c----Virtual corrections matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
C----f(p1)+f(p2) --> gamma(p3)+gamma(p4)+g(p5)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
C     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
C     Modified by overall factors
      integer iglue,j,k
      double precision p(mxpart,4),msq(fn:nf,fn:nf)
      double precision ss,tt,uu,msqgamgam,
     . virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac
      parameter(iglue=5)

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      Asq=(as/(3d0*pi))**2/vevsq

C   Deal with Higgs decay
      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)

      fac=ason2pi*Asq*gsq*hdecay
      call hjetfill(ss,tt,uu,virtgg,virtqa,virtaq,virtqg,virtgq)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      if ((j.eq.0).and.(k.eq.0)) msq(j,k)=avegg*fac*virtgg
      if ((j.gt.0).and.(k.eq.-j)) msq(j,k)=aveqq*fac*virtqa
      if ((j.lt.0).and.(k.eq.-j)) msq(j,k)=aveqq*fac*virtaq
      if ((j.eq.0).and.(k.ne.0)) msq(j,k)=aveqg*fac*virtgq
      if ((j.ne.0).and.(k.eq.0)) msq(j,k)=aveqg*fac*virtqg
      enddo
      enddo

      return
      end
