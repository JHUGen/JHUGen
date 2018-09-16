      subroutine qqb_trigam(p,msq) 
      implicit none
      include 'types.f'
************************************************************************
*    Author: C. Williams                                               *
*    March, 2013.                                                      *
*    Lowest order matrix element squared, averaged over initial colors *
*    and spins                                                         *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + gam(p5)                 *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'ewcouple.f' 
      include 'ewcharge.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf) 
      integer:: j,h1,h2,h3,h4
      complex(dp):: qqb(2,2,2,2)
      real(dp):: qqbsum
      real(dp):: fac,statfac
      parameter(statfac=one/6._dp)

c--- initialize matrix elements
      msq(:,:)=0._dp 
      
      call spinoru(5,p,za,zb)
      
c--- fill qqb helicity amplitudes
      call amp_lo_3gam(1,2,3,4,5,za,zb,qqb)
c--- note that summed, squared qbq amplitudes are identical
c      call amp_lo_3gam(2,1,3,4,5,za,zb,qbq)

      qqbsum=0._dp
      do h1=1,2 
      do h2=1,2 
      do h3=1,2 
      do h4=1,2
        qqbsum=qqbsum+real(qqb(h1,h2,h3,h4)*conjg(qqb(h1,h2,h3,h4)))
      enddo
      enddo
      enddo     
      enddo
      
c--- overall factor except for photon charge which is applied below
      fac=esq**3*xn*8._dp*aveqq*statfac 
    
      do j=-nf,nf 
         if (j .ne. 0) then 
            msq(j,-j)=fac*Q(j)**6*qqbsum
         endif
      enddo
      
      return 
      end 
      
      
      subroutine amp_lo_3gam(p1,p2,p3,p4,p5,za,zb,amp) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5
      complex(dp):: amp(2,2,2,2),trigam

!======= amplitudes that are zero
      amp(2,1,1,1)=czip 
      amp(2,2,2,2)=czip
      amp(1,2,2,2)=czip
      amp(1,1,1,1)=czip

!======= 3,4,5 symmetry 
      amp(1,2,2,1)=trigam(p1,p2,p3,p4,p5,za,zb)
      amp(1,2,1,2)=trigam(p1,p2,p3,p5,p4,za,zb)
      amp(1,1,2,2)=trigam(p1,p2,p4,p5,p3,za,zb)
       
!======= line reversal 
      amp(2,2,2,1)=trigam(p2,p1,p3,p4,p5,za,zb)
      amp(2,2,1,2)=trigam(p2,p1,p3,p5,p4,za,zb)
      amp(2,1,2,2)=trigam(p2,p1,p4,p5,p3,za,zb)

!======= conjugation 
      amp(2,1,1,2)=-trigam(p1,p2,p3,p4,p5,zb,za)
      amp(2,1,2,1)=-trigam(p1,p2,p3,p5,p4,zb,za)
      amp(2,2,1,1)=-trigam(p1,p2,p4,p5,p3,zb,za)
!======= conjugation and line reversal 
      amp(1,1,1,2)=-trigam(p2,p1,p3,p4,p5,zb,za)
      amp(1,1,2,1)=-trigam(p2,p1,p3,p5,p4,zb,za)
      amp(1,2,1,1)=-trigam(p2,p1,p4,p5,p3,zb,za)

      return 
      end 


      function trigam(p1,p2,p3,p4,p5,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: trigam
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5 
!----- amplitude for q(-,-p1),qb(+,-p2),gam(3,+),gam(4,+),gam(5,-) 
!----- all momenta outgoing 
      
      trigam=za(p2,p1)*za(p1,p5)**2/
     &     (za(p1,p3)*za(p1,p4)*za(p2,p3)*za(p2,p4))
      
      return 
      end 
