      subroutine qqb_tbbdk(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the LO element squared and subtraction terms           *
*     for the process                                                  *
*                                                                      *
*     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)                *
*     or                                                               *
*     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+bb(p5))+b(p6)                *
*                                                                      *
*     Top (antitop) is kept strictly on-shell                          *
*     although all spin correlations are retained.                     *
*                                                                      *
*     NOTE: this routine is a replacement for qqb_tbb.f that includes  *
*           the effect of the b-quark mass. In the massless case it is *
*           approximately 10% slower than that routine                 *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'nwz.f'
      integer j,k,hb,hc,jmax,jmin
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision fac,qqb,qbq
      double complex  prop
      double complex  mdecay(2,2),mprodqa(2,2),mprodaq(2,2),
     & mtotqa(2,2),mtotaq(2,2)
      parameter(jmax=1,jmin=2)

      prop=dcmplx(zip,mt*twidth)
      fac=aveqq*xn**2*gwsq**4/abs(prop)**2

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      if (nwz .eq. +1) then
        call schantoponshell(1,2,p,0,mprodqa)
        call schantoponshell(2,1,p,0,mprodaq)
        call tdecay(p,3,4,5,mdecay)
      elseif (nwz .eq. -1) then
        call schanatoponshell(1,2,p,0,mprodqa)
        call schanatoponshell(2,1,p,0,mprodaq)
        call adecay(p,3,4,5,mdecay)
      endif


c--- Calculate complete amplitudes and square
      qqb=0d0
      qbq=0d0
      do hb=1,2
      do hc=1,2
      mtotqa(hb,hc)=czip
      mtotaq(hb,hc)=czip
      if (nwz .eq. +1) then
          do j=1,jmax
          mtotqa(hb,hc)=mtotqa(hb,hc)+mdecay(hb,j)*mprodqa(j,hc)
          mtotaq(hb,hc)=mtotaq(hb,hc)+mdecay(hb,j)*mprodaq(j,hc)
          enddo
      elseif (nwz .eq. -1) then
          do j=jmin,2
          mtotqa(hb,hc)=mtotqa(hb,hc)+mprodqa(hb,j)*mdecay(j,hc)
          mtotaq(hb,hc)=mtotaq(hb,hc)+mprodaq(hb,j)*mdecay(j,hc)
          enddo
      endif
      qqb=qqb+abs(mtotqa(hb,hc))**2
      qbq=qbq+abs(mtotaq(hb,hc))**2
      enddo
      enddo

C---fill qb-q and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsq(j,k)*fac*qqb 
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsq(j,k)*fac*qbq 
      endif
      enddo
      enddo

      return
      end
