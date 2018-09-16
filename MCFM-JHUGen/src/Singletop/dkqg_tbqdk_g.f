      subroutine dkqg_tbqdk_g(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     [nwz=+1]                                                         *
*     q(-p1) +g(-p2)=nu(p3)+e+(p4)+b(p5)+bb(p6)+q'(p7)                 *
*     +g(p8) radiated from top in decay                                *
*                                                                      *
*     [nwz=-1]                                                         *
*     q(-p1) +g(-p2)=e-(p3)+nu~(p4)+bb(p5)+b(p6)+q'(p7)                *
*     +g(p8) radiated from antitop in decay                            *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*     Mass of bottom quark in decay is included.                       *
*                                                                      *
*     NOTE: this routine is a replacement for dkqg_tbqdk_g_old.f,      *
*           including the effect of the b-quark mass. In the massless  *
*           case it is approximately 3 times faster than that routine  *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'nwz.f'
      integer:: j,k,hb,hc,ht,ha,h2,hg
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,msq_qg,msq_gq,msq_qbg,msq_gqb
      complex(dp)::  prop
      complex(dp)::  mtop(2,2,2),manti(2,2,2),
     & mqg(2,2,2),mgq(2,2,2),mqbg(2,2,2),mgqb(2,2,2),
     & mtotqg(2,2,2,2),mtotgq(2,2,2,2),mtotqbg(2,2,2,2),mtotgqb(2,2,2,2)

C----set all elements to zero
      msq(:,:)=0._dp

      if (nwz == +1) then
        call singletoponshell(1,2,7,p,1,mqg)
        call singletoponshell(2,1,7,p,1,mgq)
        call singletoponshell(7,2,1,p,1,mqbg)
        call singletoponshell(7,1,2,p,1,mgqb)
        call tdecayg(p,3,4,5,8,mtop)
      else
        call singleatoponshell(1,2,7,p,-1,mqg)
        call singleatoponshell(2,1,7,p,-1,mgq)
        call singleatoponshell(7,2,1,p,-1,mqbg)
        call singleatoponshell(7,1,2,p,-1,mgqb)
        call adecayg(p,3,4,5,8,manti)
      endif

c--- q-g amplitudes
      do hb=1,2
      do h2=1,2
      do hg=1,2
      do hc=1,2
      mtotqg(hb,hg,h2,hc)=czip
      mtotgq(hb,hg,h2,hc)=czip
      mtotqbg(hb,hg,h2,hc)=czip
      mtotgqb(hb,hg,h2,hc)=czip

      if (nwz == +1) then

      do ht=1,2
      mtotqg(hb,hg,h2,hc)=mtotqg(hb,hg,h2,hc)
     & +mtop(hb,hg,ht)*mqg(ht,h2,hc)
      mtotgq(hb,hg,h2,hc)=mtotgq(hb,hg,h2,hc)
     & +mtop(hb,hg,ht)*mgq(ht,h2,hc)
      mtotqbg(hb,hg,h2,hc)=mtotqbg(hb,hg,h2,hc)
     & +mtop(hb,hg,ht)*mqbg(ht,h2,hc)
      mtotgqb(hb,hg,h2,hc)=mtotgqb(hb,hg,h2,hc)
     & +mtop(hb,hg,ht)*mgqb(ht,h2,hc)
      enddo

      else

      do ha=1,2
      mtotqg(hb,hg,h2,hc)=mtotqg(hb,hg,h2,hc)
     & +mqg(hb,h2,ha)*manti(ha,hg,hc)
      mtotgq(hb,hg,h2,hc)=mtotgq(hb,hg,h2,hc)
     & +mgq(hb,h2,ha)*manti(ha,hg,hc)
      mtotqbg(hb,hg,h2,hc)=mtotqbg(hb,hg,h2,hc)
     & +mqbg(hb,h2,ha)*manti(ha,hg,hc)
      mtotgqb(hb,hg,h2,hc)=mtotgqb(hb,hg,h2,hc)
     & +mgqb(hb,h2,ha)*manti(ha,hg,hc)
      enddo

      endif

      enddo
      enddo
      enddo
      enddo

      prop=cplx2(zip,mt*twidth)
      fac=V*xn*gwsq**4*gsq/abs(prop)**2*gsq*V/xn
c--- include factor for hadronic decays
c      if ((kcase==ktt_bbh) .or. (kcase==ktt_hdk)) fac=2._dp*xn*fac
      msq_qg=0._dp
      msq_gq=0._dp
      msq_qbg=0._dp
      msq_gqb=0._dp
      do hb=1,2
      do hg=1,2
      do h2=1,2
      do hc=1,2
      msq_qg=msq_qg+fac*aveqg*abs(mtotqg(hb,hg,h2,hc))**2
      msq_gq=msq_gq+fac*aveqg*abs(mtotgq(hb,hg,h2,hc))**2
      msq_qbg=msq_qbg+fac*aveqg*abs(mtotqbg(hb,hg,h2,hc))**2
      msq_gqb=msq_gqb+fac*aveqg*abs(mtotgqb(hb,hg,h2,hc))**2
      enddo
      enddo
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k == 0)) then
      msq(j,k)=Vsum(j)*msq_qg
      elseif ((j < 0) .and. (k == 0)) then
      msq(j,k)=Vsum(j)*msq_qbg
      elseif ((j == 0) .and. (k > 0)) then
      msq(j,k)=Vsum(k)*msq_gq
      elseif ((j == 0) .and. (k < 0)) then
      msq(j,k)=Vsum(k)*msq_gqb
      endif
      enddo
      enddo
      return
      end
