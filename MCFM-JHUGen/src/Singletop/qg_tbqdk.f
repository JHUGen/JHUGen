      subroutine qg_tbqdk(p,msq)
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
*                                                                      *
*     [nwz=-1]                                                         *
*     q(-p1) +g(-p2)=e-(p3)+v~(p4)+b~(p5)+b(p6)+q'(p7)                 *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
*     NOTE: this routine is a replacement for qg_tbqdk_old.f,          *
*           including the effect of the b-quark mass. In the massive   *
*           case it is approximately 10% slower than that routine.     *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'stopscales.f'
      include 'ckm.f'
      include 'nwz.f'
      integer:: j,k,hb,hc,ht,ha,h2,hbmax,htmax,hcmin,hamin
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,msq_qg,msq_gq,msq_qbg,msq_gqb
      complex(dp):: prop
      complex(dp):: mtop(2,2),manti(2,2),
     & mqg(2,2,2),mgq(2,2,2),mqbg(2,2,2),mgqb(2,2,2),
     & mtotqg(2,2,2),mtotgq(2,2,2),mtotqbg(2,2,2),mtotgqb(2,2,2)
      parameter(hbmax=1,htmax=1)
      parameter(hcmin=2,hamin=2)

C----set all elements to zero
      msq(:,:)=0._dp

      if (nwz == +1) then
        call singletoponshell(1,2,7,p,0,mqg)
        call singletoponshell(2,1,7,p,0,mgq)
        call singletoponshell(7,2,1,p,0,mqbg)
        call singletoponshell(7,1,2,p,0,mgqb)
        call tdecay(p,3,4,5,mtop)
      else
        call singleatoponshell(1,2,7,p,0,mqg)
        call singleatoponshell(2,1,7,p,0,mgq)
        call singleatoponshell(7,2,1,p,0,mqbg)
        call singleatoponshell(7,1,2,p,0,mgqb)
        call adecay(p,3,4,5,manti)
      endif

      mtotqg(:,:,:)=czip
      mtotgq(:,:,:)=czip
      mtotqbg(:,:,:)=czip
      mtotgqb(:,:,:)=czip

      if (nwz == +1) then
      do hb=1,hbmax
      do h2=1,2
      do hc=1,2
      do ht=1,htmax
      mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)
     & +mtop(hb,ht)*mqg(ht,h2,hc)
      mtotgq(hb,h2,hc)=mtotgq(hb,h2,hc)
     & +mtop(hb,ht)*mgq(ht,h2,hc)
      mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)
     & +mtop(hb,ht)*mqbg(ht,h2,hc)
      mtotgqb(hb,h2,hc)=mtotgqb(hb,h2,hc)
     & +mtop(hb,ht)*mgqb(ht,h2,hc)
      enddo
      enddo
      enddo
      enddo

      else

      do hb=1,2
      do h2=1,2
      do hc=hcmin,2
      do ha=hamin,2
      mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)
     & +mqg(hb,h2,ha)*manti(ha,hc)
      mtotgq(hb,h2,hc)=mtotgq(hb,h2,hc)
     & +mgq(hb,h2,ha)*manti(ha,hc)
      mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)
     & +mqbg(hb,h2,ha)*manti(ha,hc)
      mtotgqb(hb,h2,hc)=mtotgqb(hb,h2,hc)
     & +mgqb(hb,h2,ha)*manti(ha,hc)
      enddo
      enddo
      enddo
      enddo

      endif


      prop=cplx2(zip,mt*twidth)
      fac=V*xn*gwsq**4*fourpi*as_H/abs(prop)**2
c--- include factor for hadronic decays
c      if ((kcase==ktt_bbh) .or. (kcase==ktt_hdk)) fac=2._dp*xn*fac
      msq_qg=0._dp
      msq_gq=0._dp
      msq_qbg=0._dp
      msq_gqb=0._dp
      do hb=1,2
      do h2=1,2
      do hc=1,2
      msq_qg=msq_qg+abs(mtotqg(hb,h2,hc))**2
      msq_gq=msq_gq+abs(mtotgq(hb,h2,hc))**2
      msq_qbg=msq_qbg+abs(mtotqbg(hb,h2,hc))**2
      msq_gqb=msq_gqb+abs(mtotgqb(hb,h2,hc))**2
      enddo
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k == 0)) then
      msq(j,k)=Vsum(j)*aveqg*fac*msq_qg
      elseif ((j < 0) .and. (k == 0)) then
      msq(j,k)=Vsum(j)*aveqg*fac*msq_qbg
      elseif ((j == 0) .and. (k > 0)) then
      msq(j,k)=Vsum(k)*aveqg*fac*msq_gq
      elseif ((j == 0) .and. (k < 0)) then
      msq(j,k)=Vsum(k)*aveqg*fac*msq_gqb
      endif
      enddo
      enddo
      return
      end
