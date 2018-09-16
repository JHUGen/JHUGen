      subroutine dkqqb_tbbdk_g(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the LO element squared and subtraction terms           *
*     for the process                                                  *
*                                                                      *
*     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)+g(p7)          *
*     or                                                               *
*     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+b~(p5))+b(p6)+g(p7)          *
*     with the radiation coming from the top (antitop) decay           *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
*     NOTE: this routine is a replacement for qqb_tbb_gdk.f, including *
*           the effect of the b-quark mass. In the massless case it is *
*           approximately 5% slower than that routine                  *
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
      integer:: j,k,hb,hg,hc
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qqb,qbq
      complex(dp)::  prop
      complex(dp)::  mdecay(2,2,2),mprodqa(2,2),mprodaq(2,2),
     & mtotqa(2,2,2),mtotaq(2,2,2)

      prop=cplx2(zip,mt*twidth)
      fac=aveqq*xn**2*gwsq**4/abs(prop)**2*gsq*V/xn

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      if (nwz == +1) then
      call schantoponshell(1,2,p,1,mprodqa)
      call schantoponshell(2,1,p,1,mprodaq)
      call tdecayg(p,3,4,5,7,mdecay)
      elseif (nwz == -1) then
      call schanatoponshell(1,2,p,1,mprodqa)
      call schanatoponshell(2,1,p,1,mprodaq)
      call adecayg(p,3,4,5,7,mdecay)
      endif

c--- Calculate complete amplitudes and square
      qqb=0._dp
      qbq=0._dp
      do hb=1,2
      do hg=1,2
      do hc=1,2
      mtotqa(hb,hg,hc)=czip
      mtotaq(hb,hg,hc)=czip

      if (nwz == +1) then
         do j=1,2
         mtotqa(hb,hg,hc)=mtotqa(hb,hg,hc)+mdecay(hb,hg,j)*mprodqa(j,hc)
         mtotaq(hb,hg,hc)=mtotaq(hb,hg,hc)+mdecay(hb,hg,j)*mprodaq(j,hc)
         enddo
      elseif (nwz == -1) then
         do j=1,2
         mtotqa(hb,hg,hc)=mtotqa(hb,hg,hc)+mprodqa(hb,j)*mdecay(j,hg,hc)
         mtotaq(hb,hg,hc)=mtotaq(hb,hg,hc)+mprodaq(hb,j)*mdecay(j,hg,hc)
         enddo
      endif
      qqb=qqb+abs(mtotqa(hb,hg,hc))**2
      qbq=qbq+abs(mtotaq(hb,hg,hc))**2
      enddo
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
