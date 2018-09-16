      subroutine VV_HWW_gs(p,msqc)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion : sums up WW and VV contributions
************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2002.                                                      *
*                                                                      *
*     Weak Boson Fusion : sums up WW and ZZ contributions              *
*     This routine calculates the dipole subtraction terms             *
*     for the process:                                                 *
*     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
*                           |                                          *
*                           |                                          *
*                           |                                          *
*                           ---> W+(nu(p3)+e+(p4))+W-(e-(p5)+nub(p6))  *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf),
     & msqc_ww(maxd,-nf:nf,-nf:nf),msqc_zz(maxd,-nf:nf,-nf:nf)
      integer:: j,k,nd

      call WW_HWW_gs(p,msqc_ww)
      call ZZ_HWW_gs(p,msqc_zz)

      do nd=1,ndmax
      do j=-nf,nf
      do k=-nf,nf
        msqc(nd,j,k)=msqc_ww(nd,j,k)+msqc_zz(nd,j,k)
      enddo
      enddo
      enddo

      return
      end

