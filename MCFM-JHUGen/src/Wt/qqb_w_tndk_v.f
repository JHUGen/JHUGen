      subroutine qqb_w_tndk_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: J.M. Campbell                                            *
*     December, 2004.                                                  *
************************************************************************
c---- One-loop matrix element for W+t production, including t mass
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + tbar(p5)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ t(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      include 'qcdcouple.f'
      include 'scheme.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),wprop,fac
      real(dp):: virtqg,virtgq,virtqbg,virtgqb
      real(dp):: twotDg,q(mxpart,4),dot
      complex(dp):: amp0(2,2),spp,spm,smp,smm,
     & virt_pp,virt_pm,virt_mp,virt_mm

      scheme='dred'

      msq(:,:)=zero

      wprop=(s(3,4)-wmass**2)**2+(wmass*wwidth)**2

c--- calculate auxiliary momentum array - gq case
      twotDg=two*dot(p,5,1)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(1,k)*mt**2/twotDg
      enddo

c---fill matrices of spinor products
      call spinoru(6,q,za,zb)

      if     (nwz == -1) then
c---- basic process is g+b -> W- + t
c--- Note: call to tree now passes top mass as a parameter
        call tree(mt,1,2,3,4,6,amp0)
        spp=virt_pp(mt,1,2,3,4,5,q)/sqrt(wprop)
        spm=virt_pm(mt,1,2,3,4,5,q)/sqrt(wprop)
        smm=virt_mm(mt,1,2,3,4,5,q)/sqrt(wprop)
        smp=virt_mp(mt,1,2,3,4,5,q)/sqrt(wprop)

        virtgq=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &             +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtgq=virtsqwcg(2,1,3,4,5,p)/wprop   ! (with mc=mt)
      elseif (nwz == +1) then
c---- basic process is g+b~ -> W+ + t~
        call tree(mt,1,2,4,3,6,amp0)
        spp=virt_pp(mt,1,2,4,3,5,q)/sqrt(wprop)
        spm=virt_pm(mt,1,2,4,3,5,q)/sqrt(wprop)
        smm=virt_mm(mt,1,2,4,3,5,q)/sqrt(wprop)
        smp=virt_mp(mt,1,2,4,3,5,q)/sqrt(wprop)

        virtgqb=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &              +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtgqb=virtsqwcg(2,1,4,3,5,p)/wprop   ! (with mc=mt)
      else
        write(6,*) 'Problem with nwz in qqb_w_tndk.f: nwz=',nwz
        stop
      endif

c--- calculate auxiliary momentum array - qg case
      twotDg=two*dot(p,5,2)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(2,k)*mt**2/twotDg
      enddo

c---fill matrices of spinor products
      call spinoru(6,q,za,zb)

      if     (nwz == -1) then
c---- basic process is b+g -> W- + t
        call tree(mt,2,1,3,4,6,amp0)
        spp=virt_pp(mt,2,1,3,4,5,q)/sqrt(wprop)
        spm=virt_pm(mt,2,1,3,4,5,q)/sqrt(wprop)
        smm=virt_mm(mt,2,1,3,4,5,q)/sqrt(wprop)
        smp=virt_mp(mt,2,1,3,4,5,q)/sqrt(wprop)

        virtqg=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &             +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtqg=virtsqwcg(1,2,3,4,5,p)/wprop   ! (with mc=mt)
      elseif (nwz == +1) then
c---- basic process is b~+g -> W+ + t~
        call tree(mt,2,1,4,3,6,amp0)
        spp=virt_pp(mt,2,1,4,3,5,q)/sqrt(wprop)
        spm=virt_pm(mt,2,1,4,3,5,q)/sqrt(wprop)
        smm=virt_mm(mt,2,1,4,3,5,q)/sqrt(wprop)
        smp=virt_mp(mt,2,1,4,3,5,q)/sqrt(wprop)

        virtqbg=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &              +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

c        virtqbg=virtsqwcg(1,2,4,3,5,p)/wprop   ! (with mc=mt)
      endif

      fac=ason2pi*gsq*gwsq**2*V*aveqg

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j == +5) .and. (k == 0) .and. (nwz == -1)) then
          msq(j,k)=fac*virtqg
      elseif ((j == -5) .and. (k == 0) .and. (nwz == +1)) then
          msq(j,k)=fac*virtqbg
      elseif ((j == 0) .and. (k == +5) .and. (nwz == -1)) then
          msq(j,k)=fac*virtgq
      elseif ((j == 0) .and. (k == -5) .and. (nwz == +1)) then
          msq(j,k)=fac*virtgqb
      endif

      enddo
      enddo
      return
      end

