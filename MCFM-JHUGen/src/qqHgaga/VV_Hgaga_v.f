      subroutine VV_Hgaga_v(p,msqv)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: J. M. Campbell                                           *
*     June, 2002.                                                      *
*                                                                      *
*     Weak Boson Fusion : sums up WW and ZZ contributions              *
*     This routine calculates the virtual matrix element squared       *
*     for the process:                                                 *
*     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'scheme.f'
      include 'first.f'
      real(dp):: dot,p(mxpart,4),
     & msq0(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),xl15,xl26,facv
      integer:: j,k

      if (first) then
        first=.false.
      endif

c--- this result is in the MS-bar scheme
      scheme='dred'
      
* Virtual matrix elements are simply lowest order multiplied by factor
* Need to check this for pisq (cf. qqb_w_v.f) and I think there should
* be an additional +ason2pi*cf*(1+1) from DRED compared to MSBAR
      call VV_Hgaga(p,msq0)
      xl15=log(-2._dp*dot(p,1,5)/musq)
      xl26=log(-2._dp*dot(p,2,6)/musq)
      facv=ason2pi*cf*(
     & -4._dp*EPINV*EPINV2-(6._dp-2._dp*(xl15+xl26))*EPINV
     & +3._dp*(xl15+xl26)-(xl15**2+xl26**2)-14._dp)
      
      do j=-nf,nf
      do k=-nf,nf
        msqv(j,k)=facv*msq0(j,k)
      enddo
      enddo

      
      return
      end
