      subroutine qqb_tottth_v(p,msq)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared 
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=t(p3)+t(p4)+h(p5)
C  
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
      include 'epinv.f'
      include 'swapxz.f'
      include 'msbarmasses.f'
      include 'msq_cs.f'
      include 'scheme.f'
      include 'first.f'
      
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ren,
     & wtqqb,wtqbq,wtgg0,massfrun,mt_eff,facqq,facgg,ytsq
      complex(dp):: cnab(-2:-1),cqed(-2:-1)
      save mt_eff
!$omp threadprivate(mt_eff) 

      scheme='dred'     
      
      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

      call ttqqHampsq(p,3,4,1,2,wtqqb)
      call ttqqHampsq(p,3,4,2,1,wtqbq)
      call singcoeffq(p,3,4,1,2,cnab,cqed)
      wtqqb=(real(cnab(-2))*epinv**2+real(cnab(-1))*epinv)*wtqqb
c      write(6,*) 'cnab(-2),cqed(-2)',cnab(-2),cqed(-2)
c      write(6,*) 'cnab(-1),cqed(-1)',cnab(-1),cqed(-1)
      call singcoeffq(p,3,4,2,1,cnab,cqed)
      wtqbq=(real(cnab(-2))*epinv**2+real(cnab(-1))*epinv)*wtqbq
c      write(6,*) 'cnab(-2),cqed(-2)',cnab(-2),cqed(-2)
c      write(6,*) 'cnab(-1),cqed(-1)',cnab(-1),cqed(-1)
      
      
      
      call ttggHsingdriver(p,wtgg0)

      ytsq=0.5d0*gwsq*(mt_eff/wmass)**2
      facqq=ason2pi*V/8d0*gsq**2*ytsq
      facgg=facqq*xn


c--- renormalization
c      ren=(
c     & +2d0*((nlf/xn-11d0/6d0*xn)*epinv+xn/6d0)
c     & +2d0*nhf/xn*(epinv+log(musq/mt**2))
c     & -Cf*(3d0*(epinv+log(musq/mt**2))+5d0))
c      resqqb1=resqqb1+ren*resqqb0
c      resqbq1=resqbq1+ren*resqqb0
     
c      ren=(
c     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
c     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/mt**2))+5d0))
c      resgg1=resgg1+ren*resgg0
      
C----set all elements to zero
      msq(:,:)=0d0

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j < 0) then
          msq(j,-j)=aveqq*facqq*wtqqb
      elseif (j == 0) then
          msq(j,j)=avegg*facgg*wtgg0
      elseif (j > 0) then
          msq(j,-j)=aveqq*facqq*wtqqb
      endif
      enddo
      
      return
      end
 
