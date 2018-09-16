      subroutine qqb_totttH_gvec(p,n,in,msqv)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     May, 2013.                                                       *
*----Matrix element for ttH production                                 *
*----averaged over initial colours and spins                           *
*    line in contracted with the vector n(mu)                          *
*     g(-p1)+g(-p2)--> t(p3)+tb(p4)+H(p5)                              *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'kpart.f'
      include 'couple.f'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'swapxz.f'
      include 'msqv_cs.f'
      include 'first.f'
C in is the label of the contracted line
      integer:: j,k,in,icol
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),res(0:2)
      real(dp):: n(4),fac,massfrun,mt_eff,ytsq
      save mt_eff
!$omp threadprivate(mt_eff)
      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

      call checkndotp(p,n,in)


      ytsq=0.5d0*gwsq*(mt_eff/wmass)**2
      fac=V/8d0*xn*gsq**2*ytsq

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      do icol=0,2
        msqv_cs(icol,j,k)=0d0
      enddo
      enddo
      enddo

      call ttggHampsqn(p,3,4,1,2,1,1,in,n,res)


      do icol=0,2
      msqv_cs(icol,0,0)=avegg*fac*res(icol)
      enddo
      msqv(0,0)=msqv_cs(0,0,0)+msqv_cs(1,0,0)+msqv_cs(2,0,0)
      
      return
      end
