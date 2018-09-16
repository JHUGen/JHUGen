      subroutine qqb_w_twdk_v(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     Virtual matrix element squared and averaged over initial         *
*      colours and spins (corrections in production)                   *
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'nwz.f'
      integer:: j,k,i3,i4,i5,i6,iq
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & q(mxpart,4),msq_gq,msq_qg,fac,r(mxpart,4),wprop
      complex(dp):: ampl0(2),amplv(2),amplp(2,2),ampld(2)
      real(dp):: twotDg,dot
      complex(dp):: spp_ft,spm_ft,
     & smm_ft,smp_ft,virt_pp,virt_pm,virt_mp,virt_mm

c---initialize
      msq_gq=zero
      msq_qg=zero
      msq(:,:)=zero


c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        i3=3
        i4=4
        i5=5
        i6=6
        iq=1 ! top quark
      elseif (nwz == +1) then
        i3=4
        i4=3
        i5=6
        i6=5
        iq=-1 ! antitop quark
      else
        write(6,*) 'Error in qqb_w_twdk_v, nwz is not +1 or -1 :   ',nwz
        stop
      endif

c--- overall factor contained in diag.prc
      fac=ason2pi*aveqg*(V/two)*two*gsq*gwsq**4

c-- calculate auxiliary momentum array - gq case
      twotDg=two*(dot(p,5,1)+dot(p,6,1)+dot(p,7,1))
      do k=1,4
      do j=1,7
      q(j,k)=p(j,k)
      enddo
      q(8,k)=p(5,k)+p(6,k)+p(7,k)-p(1,k)*mt**2/twotDg
      enddo

      do k=1,4
      do j=1,4
      r(j,k)=p(j,k)
      enddo
      r(5,k)=p(5,k)+p(6,k)+p(7,k)
      enddo

c---fill matrices of spinor products
      call spinoru(8,q,za,zb)

c--- Note: call to tree now passes top mass as a parameter
      call tree(mt,1,2,i3,i4,8,amplp)
      call wampd(mt,twidth,1,i5,i6,7,8,ampld)

      spp_ft=virt_pp(mt,1,2,i3,i4,5,r)
      spm_ft=virt_pm(mt,1,2,i3,i4,5,r)
      smm_ft=virt_mm(mt,1,2,i3,i4,5,r)
      smp_ft=virt_mp(mt,1,2,i3,i4,5,r)

      wprop=sqrt((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

c--- we have demonstrated that amplp(1,1)=lomm_ft/wprop etc.

c      write(6,*) 'amplp(1,1)',amplp(1,1)
c      write(6,*) 'amplp(1,2)',amplp(1,2)
c      write(6,*) 'amplp(2,1)',amplp(2,1)
c      write(6,*) 'amplp(2,2)',amplp(2,2)
c      write(6,*)'lomm_ft',lomm_ft/wprop
c      write(6,*)'lomp_ft',lomp_ft/wprop
c      write(6,*)'lopm_ft',lopm_ft/wprop
c      write(6,*)'lopp_ft',lopp_ft/wprop

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark, for Born and virtual

c--- lowest order amplitudes
      ampl0(1)=amplp(1,1)*ampld(1)
     &        +amplp(2,1)*ampld(2)

      ampl0(2)=amplp(1,2)*ampld(1)
     &        +amplp(2,2)*ampld(2)

c--- virtual amplitudes
      amplv(1)=(smm_ft/wprop)*ampld(1)
     &        +(spm_ft/wprop)*ampld(2)

      amplv(2)=(smp_ft/wprop)*ampld(1)
     &        +(spp_ft/wprop)*ampld(2)

      do j=1,2
      msq_gq=msq_gq+real(amplv(j)*conjg(ampl0(j)))
      enddo

c-- calculate auxiliary momentum array - qg case
      twotDg=two*(dot(p,5,2)+dot(p,6,2)+dot(p,7,2))
      do k=1,4
c      do j=1,7
c      q(j,k)=p(j,k)
c      enddo
      q(8,k)=p(5,k)+p(6,k)+p(7,k)-p(2,k)*mt**2/twotDg
      enddo

c      do k=1,4
c      do j=1,4
c      r(j,k)=p(j,k)
c      enddo
c      r(5,k)=p(5,k)+p(6,k)+p(7,k)
c      enddo

c---fill matrices of spinor products
      call spinoru(8,q,za,zb)

      call tree(mt,2,1,i3,i4,8,amplp)
      call wampd(mt,twidth,2,i5,i6,7,8,ampld)

      spp_ft=virt_pp(mt,2,1,i3,i4,5,r)
      spm_ft=virt_pm(mt,2,1,i3,i4,5,r)
      smm_ft=virt_mm(mt,2,1,i3,i4,5,r)
      smp_ft=virt_mp(mt,2,1,i3,i4,5,r)

      wprop=sqrt((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

c--- we have demonstrated that amplp(1,1)=lomm_ft/wprop etc.

c      write(6,*) 'amplp(1,1)',amplp(1,1)
c      write(6,*) 'amplp(1,2)',amplp(1,2)
c      write(6,*) 'amplp(2,1)',amplp(2,1)
c      write(6,*) 'amplp(2,2)',amplp(2,2)
c      write(6,*)'lomm_ft',lomm_ft/wprop
c      write(6,*)'lomp_ft',lomp_ft/wprop
c      write(6,*)'lopm_ft',lopm_ft/wprop
c      write(6,*)'lopp_ft',lopp_ft/wprop

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark, for Born and virtual

c--- lowest order amplitudes
      ampl0(1)=amplp(1,1)*ampld(1)
     &        +amplp(2,1)*ampld(2)

      ampl0(2)=amplp(1,2)*ampld(1)
     &        +amplp(2,2)*ampld(2)

c--- virtual amplitudes
      amplv(1)=(smm_ft/wprop)*ampld(1)
     &        +(spm_ft/wprop)*ampld(2)

      amplv(2)=(smp_ft/wprop)*ampld(1)
     &        +(spp_ft/wprop)*ampld(2)

      do j=1,2
      msq_qg=msq_qg+real(amplv(j)*conjg(ampl0(j)))
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msq(j,k)=zero
      if     ((j == +5*iq) .and. (k == 0)) then
          msq(j,k)=fac*msq_qg
      elseif ((j == 0) .and. (k == +5*iq)) then
          msq(j,k)=fac*msq_gq
      endif
      enddo
      enddo

      return
      end


