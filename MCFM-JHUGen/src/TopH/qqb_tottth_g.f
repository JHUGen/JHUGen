      subroutine qqb_tottth_g(p,msq)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     May, 2013.                                                       *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     q(-p1) +qbar(-p2)=t(p3)+t(p4)+h(p5)+g(p6)                        *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'kpart.f'
      include 'masses.f'
      include 'swapxz.f'
      include 'msbarmasses.f'
      include 'first.f'
      
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: wtqqb,wtgg,wtqbq,wtqg,wtqbg,wtgq,wtgqb
      real(dp):: massfrun,mt_eff,ytsq,fac
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

      call ttgggHdriver(p,wtgg)
      ytsq=0.5d0*gwsq*(mt_eff/wmass)**2
      fac=V*gsq**3*ytsq/16d0

      call ttqqgHampsq(p,3,4,2,1,6,1,1,wtqbq)
      call ttqqgHampsq(p,3,4,1,2,6,1,1,wtqqb)
      call ttqqgHampsq(p,3,4,2,6,1,1,1,wtgq)
      call ttqqgHampsq(p,3,4,6,2,1,1,1,wtgqb)
      call ttqqgHampsq(p,3,4,1,6,2,1,1,wtqg)
      call ttqqgHampsq(p,3,4,6,1,2,1,1,wtqbg)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k == -j)) then
          msq(j,k)=aveqq*fac*wtqqb             !qqb
      elseif ((j < 0) .and. (k == -j)) then
          msq(j,k)=aveqq*fac*wtqbq             !qbq

      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=aveqg*fac*wtgq              !gq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=aveqg*fac*wtgqb             !gqb
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=aveqg*fac*wtqg              !qg
      elseif ((j < 0) .and. (k ==0)) then
          msq(j,k)=aveqg*fac*wtqbg             !qbg

      elseif ((j == 0)  .and. (k == 0)) then
          msq(j,k)=avegg*fac*wtgg        !gg
      endif
      enddo
      enddo
      return
      end
 
