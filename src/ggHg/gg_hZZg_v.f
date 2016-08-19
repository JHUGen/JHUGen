      subroutine gg_hZZg_v(p,msq)
C----Author: R.K. Ellis Dec 2009
c----Virtual corrections matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  Z (e^-(p5)+e^(p6)) 
c                          + Z(mu^-(p3)+mu^+(p4))+g(p_iglue1=7)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
C     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
C     Modified by overall factors
      integer iglue,j,k
      double precision p(mxpart,4),msq(fn:nf,fn:nf),s3456
      double precision ss,tt,uu,
     . virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac
      parameter(iglue=7)

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      Asq=(as/(3d0*pi))**2/vevsq

      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*zmass**2*4d0*xw**2/(one-xw)*
     . ( ((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     .  +((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))
      hdecay=hdecay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s3456-hmass**2)**2+(hmass*hwidth)**2)
      
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
