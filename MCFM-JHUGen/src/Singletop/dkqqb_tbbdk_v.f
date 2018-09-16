      subroutine dkqqb_tbbdk_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the interference of virtual corrections                *
*     and integrated subtraction term                                  *
*     for the process                                                  *
*                                                                      *
*     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)                *
*     or                                                               *
*     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+bb(p5))+b(p6)                *
*                                                                      *
*     Top (antitop) is kept strictly on-shell                          *
*     although all spin correlations are retained.                     *
*                                                                      *
*     NOTE: this routine is a replacement for qqb_tbb_vdk.f, including *
*           the effect of the b-quark mass. In the massless case it is *
*           approximately the same speed as that routine               *
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
      integer:: j,k,hb,hc
c      integer:: jmax,jmin
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qqb,qbq
      complex(dp)::  prop
      complex(dp)::  mdecay(2,2),mdecayv(2,2),
     & mprodqa(2,2),mprodaq(2,2),
     & mtotqa(2,2),mtotaq(2,2),mtotqav(2,2),mtotaqv(2,2)
      prop=cplx2(zip,mt*twidth)
      fac=aveqq*xn**2*gwsq**4/abs(prop)**2*ason2pi*CF

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      if (nwz == +1) then
        call schantoponshell(1,2,p,0,mprodqa)
        call schantoponshell(2,1,p,0,mprodaq)
        call tdecay(p,3,4,5,mdecay)
        call tdecay_v(p,3,4,5,mdecayv)
      elseif (nwz == -1) then
        call schanatoponshell(1,2,p,0,mprodqa)
        call schanatoponshell(2,1,p,0,mprodaq)
        call adecay(p,3,4,5,mdecay)
        call adecay_v(p,3,4,5,mdecayv)
      endif


c--- Calculate complete amplitudes and square
      qqb=0._dp
      qbq=0._dp
      do hb=1,2
      do hc=1,2
      mtotqa(hb,hc)=czip
      mtotaq(hb,hc)=czip
      mtotqav(hb,hc)=czip
      mtotaqv(hb,hc)=czip
      if (nwz == +1) then
        do j=1,2
        mtotqa(hb,hc)=mtotqa(hb,hc)+mdecay(hb,j)*mprodqa(j,hc)
        mtotaq(hb,hc)=mtotaq(hb,hc)+mdecay(hb,j)*mprodaq(j,hc)
        mtotqav(hb,hc)=mtotqav(hb,hc)+mdecayv(hb,j)*mprodqa(j,hc)
        mtotaqv(hb,hc)=mtotaqv(hb,hc)+mdecayv(hb,j)*mprodaq(j,hc)
        enddo
      elseif (nwz == -1) then
        do j=1,2
        mtotqa(hb,hc)=mtotqa(hb,hc)+mprodqa(hb,j)*mdecay(j,hc)
        mtotaq(hb,hc)=mtotaq(hb,hc)+mprodaq(hb,j)*mdecay(j,hc)
        mtotqav(hb,hc)=mtotqav(hb,hc)+mprodqa(hb,j)*mdecayv(j,hc)
        mtotaqv(hb,hc)=mtotaqv(hb,hc)+mprodaq(hb,j)*mdecayv(j,hc)
        enddo
      endif
      qqb=qqb+real(conjg(mtotqa(hb,hc))*mtotqav(hb,hc))
      qbq=qbq+real(conjg(mtotaq(hb,hc))*mtotaqv(hb,hc))
      enddo
      enddo

C---fill qb-q and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=Vsq(j,k)*fac*qqb
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=Vsq(j,k)*fac*qbq
      endif
      enddo
      enddo

      return
      end
