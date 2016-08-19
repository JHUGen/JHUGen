      subroutine qqb_w_tndk_v(p,msq)
      implicit none
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
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      include 'qcdcouple.f'
      include 'scheme.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),wprop,fac
      double precision virtqg,virtgq,virtqbg,virtgqb
      double precision twotDg,q(mxpart,4),dot
      double complex amp0(2,2),spp,spm,smp,smm,
     . virt_pp,virt_pm,virt_mp,virt_mm
      
      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      wprop=(s(3,4)-wmass**2)**2+(wmass*wwidth)**2

c--- calculate auxiliary momentum array - gq case
      twotDg=2d0*dot(p,5,1)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(1,k)*mt**2/twotDg
      enddo
      
c---fill matrices of spinor products  
      call spinoru(6,q,za,zb)

      if     (nwz .eq. -1) then
c---- basic process is g+b -> W- + t 
c--- Note: call to tree now passes top mass as a parameter
        call tree(mt,1,2,3,4,6,amp0)
        spp=virt_pp(mt,1,2,3,4,5,q)/dsqrt(wprop)
        spm=virt_pm(mt,1,2,3,4,5,q)/dsqrt(wprop)
        smm=virt_mm(mt,1,2,3,4,5,q)/dsqrt(wprop)
        smp=virt_mp(mt,1,2,3,4,5,q)/dsqrt(wprop)
     
        virtgq=dble(smm*Dconjg(amp0(1,1))+spp*Dconjg(amp0(2,2))
     .             +smp*Dconjg(amp0(1,2))+spm*Dconjg(amp0(2,1)))

c        virtgq=virtsqwcg(2,1,3,4,5,p)/wprop   ! (with mc=mt)
      elseif (nwz .eq. +1) then
c---- basic process is g+b~ -> W+ + t~ 
        call tree(mt,1,2,4,3,6,amp0)
        spp=virt_pp(mt,1,2,4,3,5,q)/dsqrt(wprop)
        spm=virt_pm(mt,1,2,4,3,5,q)/dsqrt(wprop)
        smm=virt_mm(mt,1,2,4,3,5,q)/dsqrt(wprop)
        smp=virt_mp(mt,1,2,4,3,5,q)/dsqrt(wprop)
     
        virtgqb=dble(smm*Dconjg(amp0(1,1))+spp*Dconjg(amp0(2,2))
     .              +smp*Dconjg(amp0(1,2))+spm*Dconjg(amp0(2,1)))

c        virtgqb=virtsqwcg(2,1,4,3,5,p)/wprop   ! (with mc=mt)
      else
        write(6,*) 'Problem with nwz in qqb_w_tndk.f: nwz=',nwz
        stop
      endif
      
c--- calculate auxiliary momentum array - qg case
      twotDg=2d0*dot(p,5,2)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(2,k)*mt**2/twotDg
      enddo
      
c---fill matrices of spinor products  
      call spinoru(6,q,za,zb)

      if     (nwz .eq. -1) then
c---- basic process is b+g -> W- + t 
        call tree(mt,2,1,3,4,6,amp0)
        spp=virt_pp(mt,2,1,3,4,5,q)/dsqrt(wprop)
        spm=virt_pm(mt,2,1,3,4,5,q)/dsqrt(wprop)
        smm=virt_mm(mt,2,1,3,4,5,q)/dsqrt(wprop)
        smp=virt_mp(mt,2,1,3,4,5,q)/dsqrt(wprop)
     
        virtqg=dble(smm*Dconjg(amp0(1,1))+spp*Dconjg(amp0(2,2))
     .             +smp*Dconjg(amp0(1,2))+spm*Dconjg(amp0(2,1)))

c        virtqg=virtsqwcg(1,2,3,4,5,p)/wprop   ! (with mc=mt)
      elseif (nwz .eq. +1) then
c---- basic process is b~+g -> W+ + t~ 
        call tree(mt,2,1,4,3,6,amp0)
        spp=virt_pp(mt,2,1,4,3,5,q)/dsqrt(wprop)
        spm=virt_pm(mt,2,1,4,3,5,q)/dsqrt(wprop)
        smm=virt_mm(mt,2,1,4,3,5,q)/dsqrt(wprop)
        smp=virt_mp(mt,2,1,4,3,5,q)/dsqrt(wprop)
     
        virtqbg=dble(smm*Dconjg(amp0(1,1))+spp*Dconjg(amp0(2,2))
     .              +smp*Dconjg(amp0(1,2))+spm*Dconjg(amp0(2,1)))

c        virtqbg=virtsqwcg(1,2,4,3,5,p)/wprop   ! (with mc=mt)
      endif
      
      fac=ason2pi*gsq*gwsq**2*V*aveqg
      
      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j .eq. +5) .and. (k .eq. 0) .and. (nwz .eq. -1)) then
          msq(j,k)=fac*virtqg
      elseif ((j .eq. -5) .and. (k .eq. 0) .and. (nwz .eq. +1)) then
          msq(j,k)=fac*virtqbg
      elseif ((j .eq. 0) .and. (k .eq. +5) .and. (nwz .eq. -1)) then
          msq(j,k)=fac*virtgq
      elseif ((j .eq. 0) .and. (k .eq. -5) .and. (nwz .eq. +1)) then
          msq(j,k)=fac*virtgqb
      endif

      enddo
      enddo
      return
      end
 
