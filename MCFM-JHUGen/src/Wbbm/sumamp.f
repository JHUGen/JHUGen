      subroutine sumamp(coeff,scints,amp,str)
      implicit none
      include 'types.f'
c--- routine to multiply scalar integrals (in scints) by the computed
c--- coefficients (in coeff) and return the result in amp
c--- the 6-character string "str" is only used as output when checking primitives
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'Wbbmlabels.f'
      integer:: iep,j,k
      character*6 str
      complex(dp):: amp(-2:0)
      logical:: numcheck
      common/numcheck/numcheck
!$omp threadprivate(/numcheck/)

c--- multiply scalar integrals by coefficients
c---  NB: only need sum over finite pieces, poles handled separately
c      do iep=-2,0
      do iep=0,0
      amp(iep)=czip
      do j=1,4
      do k=1,20
      amp(iep)=amp(iep)+coeff(j,k)*scints(j,k,iep)
c      write(6,*) iep,j,k,coeff(j,k),scints(j,k,iep)
      enddo
      enddo      
      enddo     
c--- add purely rational term
      amp(0)=amp(0)+coeff(0,irat)

      if (numcheck) then
        write(6,89) str//'(-2) =',amp(-2) 
        write(6,89) str//'(-1) =',amp(-1) 
        write(6,89) str//'( 0) =',amp( 0)
        write(6,*)
      endif

      return
      
   89 format(SP,a13,2e20.11,' i')   
      
      end
      
