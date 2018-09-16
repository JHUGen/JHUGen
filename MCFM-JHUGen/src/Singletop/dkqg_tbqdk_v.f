      subroutine dkqg_tbqdk_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the contribution of virtual corrections to the process       *
*                                                                      *
*     [nwz=+1]                                                         *
*     q(-p1) +g(-p2)=nu(p3)+e+(p4)+b(p5)+bb(p6)+q'(p7)                 *
*                                                                      *
*     [nwz=-1]                                                         *
*     q(-p1) +g(-p2)=e-(p3)+v~(p4)+b~(p5)+b(p6)+q'(p7)                 *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
*     NOTE: this routine is a replacement for dkqg_tbqdk_v_old.f,      *
*           including the effect of the b-quark mass. In the massless  *
*           case it is approximately 10% faster than that routine      *
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
      integer:: j,k,hb,hc,ht,ha,h2
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,msq_qg,msq_gq,msq_qbg,msq_gqb
      complex(dp)::  prop
      complex(dp)::  mtop(2,2),mtopv(2,2),manti(2,2),mantiv(2,2),
     & mqg(2,2,2),mgq(2,2,2),mqbg(2,2,2),mgqb(2,2,2),
     & mtotqg(2,2,2),mtotgq(2,2,2),
     & mtotqbg(2,2,2),mtotgqb(2,2,2),
     & mtotqgv(2,2,2),mtotgqv(2,2,2),
     & mtotqbgv(2,2,2),mtotgqbv(2,2,2)

C----set all elements to zero
      msq(:,:)=0._dp

      if (nwz == +1) then
        call singletoponshell(1,2,7,p,0,mqg)
        call singletoponshell(2,1,7,p,0,mgq)
        call singletoponshell(7,2,1,p,0,mqbg)
        call singletoponshell(7,1,2,p,0,mgqb)
        call tdecay(p,3,4,5,mtop)
        call tdecay_v(p,3,4,5,mtopv)
      else
        call singleatoponshell(1,2,7,p,0,mqg)
        call singleatoponshell(2,1,7,p,0,mgq)
        call singleatoponshell(7,2,1,p,0,mqbg)
        call singleatoponshell(7,1,2,p,0,mgqb)
        call adecay(p,3,4,5,manti)
        call adecay_v(p,3,4,5,mantiv)
      endif

c--- q-g amplitudes
      do hb=1,2
      do h2=1,2
      do hc=1,2
      mtotqg(hb,h2,hc)=czip
      mtotgq(hb,h2,hc)=czip
      mtotqbg(hb,h2,hc)=czip
      mtotgqb(hb,h2,hc)=czip
      mtotqgv(hb,h2,hc)=czip
      mtotgqv(hb,h2,hc)=czip
      mtotqbgv(hb,h2,hc)=czip
      mtotgqbv(hb,h2,hc)=czip

      if (nwz == +1) then

      do ht=1,2
      mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)
     & +mtop(hb,ht)*mqg(ht,h2,hc)
      mtotgq(hb,h2,hc)=mtotgq(hb,h2,hc)
     & +mtop(hb,ht)*mgq(ht,h2,hc)
      mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)
     & +mtop(hb,ht)*mqbg(ht,h2,hc)
      mtotgqb(hb,h2,hc)=mtotgqb(hb,h2,hc)
     & +mtop(hb,ht)*mgqb(ht,h2,hc)

      mtotqgv(hb,h2,hc)=mtotqgv(hb,h2,hc)
     & +mtopv(hb,ht)*mqg(ht,h2,hc)
      mtotgqv(hb,h2,hc)=mtotgqv(hb,h2,hc)
     & +mtopv(hb,ht)*mgq(ht,h2,hc)
      mtotqbgv(hb,h2,hc)=mtotqbgv(hb,h2,hc)
     & +mtopv(hb,ht)*mqbg(ht,h2,hc)
      mtotgqbv(hb,h2,hc)=mtotgqbv(hb,h2,hc)
     & +mtopv(hb,ht)*mgqb(ht,h2,hc)
      enddo

      else

      do ha=1,2
      mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)
     & +mqg(hb,h2,ha)*manti(ha,hc)
      mtotgq(hb,h2,hc)=mtotgq(hb,h2,hc)
     & +mgq(hb,h2,ha)*manti(ha,hc)
      mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)
     & +mqbg(hb,h2,ha)*manti(ha,hc)
      mtotgqb(hb,h2,hc)=mtotgqb(hb,h2,hc)
     & +mgqb(hb,h2,ha)*manti(ha,hc)

      mtotqgv(hb,h2,hc)=mtotqgv(hb,h2,hc)
     & +mqg(hb,h2,ha)*mantiv(ha,hc)
      mtotgqv(hb,h2,hc)=mtotgqv(hb,h2,hc)
     & +mgq(hb,h2,ha)*mantiv(ha,hc)
      mtotqbgv(hb,h2,hc)=mtotqbgv(hb,h2,hc)
     & +mqbg(hb,h2,ha)*mantiv(ha,hc)
      mtotgqbv(hb,h2,hc)=mtotgqbv(hb,h2,hc)
     & +mgqb(hb,h2,ha)*mantiv(ha,hc)
      enddo

      endif

      enddo
      enddo
      enddo

      prop=cplx2(zip,mt*twidth)
      fac=V*xn*gwsq**4*gsq/abs(prop)**2*ason2pi*CF
c--- include factor for hadronic decays
c      if ((kcase==ktt_bbh) .or. (kcase==ktt_hdk)) fac=2._dp*xn*fac
      msq_qg=0._dp
      msq_gq=0._dp
      msq_qbg=0._dp
      msq_gqb=0._dp
      do hb=1,2
      do h2=1,2
      do hc=1,2
      msq_qg=msq_qg+fac*aveqg
     & *real(conjg(mtotqg(hb,h2,hc))*mtotqgv(hb,h2,hc))
      msq_gq=msq_gq+fac*aveqg
     & *real(conjg(mtotgq(hb,h2,hc))*mtotgqv(hb,h2,hc))
      msq_qbg=msq_qbg+fac*aveqg
     & *real(conjg(mtotqbg(hb,h2,hc))*mtotqbgv(hb,h2,hc))
      msq_gqb=msq_gqb+fac*aveqg
     & *real(conjg(mtotgqb(hb,h2,hc))*mtotgqbv(hb,h2,hc))
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
