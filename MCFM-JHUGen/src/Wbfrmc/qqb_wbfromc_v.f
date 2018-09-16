      subroutine qqb_wbfromc_v(p,msq)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     Nov, 2011.                                                       *
************************************************************************
c---- One-loop matrix element for W+b production, including b mass
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + bbar(p5)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ b(p5) 
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
      include 'ckm.f'
      include 'nflav.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),wprop,fac
      real(dp):: virtqg,virtgq,virtqbg,virtgqb
      real(dp):: twotDg,q(mxpart,4),dot
      complex(dp):: amp0(2,2),spp,spm,smp,smm,
     & virt_pp,virt_pm,virt_mp,virt_mm
      
      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      wprop=(s(3,4)-wmass**2)**2+(wmass*wwidth)**2

c--- calculate auxiliary momentum array - gq case
      twotDg=2._dp*dot(p,5,1)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(1,k)*mb**2/twotDg
      enddo
      
c---fill matrices of spinor products  
      call spinoru(6,q,za,zb)

      if     (nwz == +1) then
c---- basic process is g+c -> W+ + b 
        call tree(mb,1,2,3,4,6,amp0)
        spp=virt_pp(mb,1,2,3,4,5,q)/sqrt(wprop)
        spm=virt_pm(mb,1,2,3,4,5,q)/sqrt(wprop)
        smm=virt_mm(mb,1,2,3,4,5,q)/sqrt(wprop)
        smp=virt_mp(mb,1,2,3,4,5,q)/sqrt(wprop)
     
        virtgq=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &             +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

      elseif (nwz == -1) then
c---- basic process is g+c~ -> W- + b~ 
        call tree(mb,1,2,4,3,6,amp0)
        spp=virt_pp(mb,1,2,4,3,5,q)/sqrt(wprop)
        spm=virt_pm(mb,1,2,4,3,5,q)/sqrt(wprop)
        smm=virt_mm(mb,1,2,4,3,5,q)/sqrt(wprop)
        smp=virt_mp(mb,1,2,4,3,5,q)/sqrt(wprop)
     
        virtgqb=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &              +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

      else
        write(6,*) 'Problem with nwz in qqb_w_bjet_v.f: nwz=',nwz
        stop
      endif
      
c--- calculate auxiliary momentum array - qg case
      twotDg=2._dp*dot(p,5,2)
      do k=1,4
      do j=1,5
      q(j,k)=p(j,k)
      enddo
      q(6,k)=p(5,k)-p(2,k)*mb**2/twotDg
      enddo
      
c---fill matrices of spinor products  
      call spinoru(6,q,za,zb)

      if     (nwz == +1) then
c---- basic process is c+g -> W+ + b 
        call tree(mb,2,1,3,4,6,amp0)
        spp=virt_pp(mb,2,1,3,4,5,q)/sqrt(wprop)
        spm=virt_pm(mb,2,1,3,4,5,q)/sqrt(wprop)
        smm=virt_mm(mb,2,1,3,4,5,q)/sqrt(wprop)
        smp=virt_mp(mb,2,1,3,4,5,q)/sqrt(wprop)
     
        virtqg=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &             +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

      elseif (nwz == -1) then
c---- basic process is c~+g -> W- + b~ 
        call tree(mb,2,1,4,3,6,amp0)
        spp=virt_pp(mb,2,1,4,3,5,q)/sqrt(wprop)
        spm=virt_pm(mb,2,1,4,3,5,q)/sqrt(wprop)
        smm=virt_mm(mb,2,1,4,3,5,q)/sqrt(wprop)
        smp=virt_mp(mb,2,1,4,3,5,q)/sqrt(wprop)
     
        virtqbg=real(smm*conjg(amp0(1,1))+spp*conjg(amp0(2,2))
     &              +smp*conjg(amp0(1,2))+spm*conjg(amp0(2,1)))

      endif
      
      fac=ason2pi*gsq*gwsq**2*V*aveqg
      
      do j=-nflav,nflav
      do k=-nflav,nflav
      if (((j == 2) .or. (j == 4)) .and. (k == 0)) then
          msq(j,k)=fac*Vsq(j,-5)*virtqg
      elseif ((j == 0) .and. ((k == +2).or.(k == +4))) then
          msq(j,k)=fac*Vsq(-5,k)*virtgq
      elseif (((j == -2).or.(j == -4)) .and. (k == 0))then
          msq(j,k)=fac*Vsq(j,+5)*virtqbg
      elseif ((j == 0) .and. ((k == -2).or.(k == -4))) then
          msq(j,k)=fac*Vsq(+5,k)*virtgqb
      endif
      enddo
      enddo
       
      return
      end
 
