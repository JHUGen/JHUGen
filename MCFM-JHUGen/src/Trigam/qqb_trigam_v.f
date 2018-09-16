      subroutine qqb_trigam_v(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Authors: C. Williams and J. M. Campbell                           *
*    March, 2013.                                                      *
*    Virtual matrix element squared, averaged over initial colors      *
*    and spins                                                         *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + gam(p5)                 *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'scheme.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qqb_lo(2,2,2,2),qqb_v(2,2,2,2)
      real(dp):: qqbsum
      integer:: j,h1,h2,h3,h4
      real(dp):: fac,statfac
      parameter(statfac=one/6._dp)
      
      scheme='dred'

c--- initialize matrix elements
      msq(:,:)=zip
      
      call spinoru(5,p,za,zb)

c--- fill qqb LO helicity amplitudes
      call amp_lo_3gam(1,2,3,4,5,za,zb,qqb_lo)
c--- fill qqb virtual helicity amplitudes
      call amp_virt_3gam(1,2,3,4,5,za,zb,qqb_v)
c--- note that the result for the qbq case will be identical

      qqbsum=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
        qqbsum=qqbsum
     &   +Dble(conjg(qqb_lo(h1,h2,h3,h4))*qqb_v(h1,h2,h3,h4))
      enddo
      enddo
      enddo
      enddo
      
      fac=esq**3*xn*8._dp*aveqq*statfac*Cf*ason2pi*2._dp

      do j=-nf,nf 
         if (j .ne. 0) then
            msq(j,-j)=fac*Q(j)**6*qqbsum
         endif
      enddo
      
      return
      end


      subroutine amp_virt_3gam(p1,p2,p3,p4,p5,za,zb,amp)
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5
      complex(dp):: amp(2,2,2,2),virt_trigam

!======= amplitudes that are zero
      amp(2,1,1,1)=czip 
      amp(2,2,2,2)=czip
      amp(1,2,2,2)=czip
      amp(1,1,1,1)=czip

c--- DEBUG ONLY
c      write(6,*)
c      amp(1,2,2,1)=virt_trigam_MHV(p1,p2,p3,p4,p5,za,zb)
c      amp(1,2,2,1)=virt_trigam(p1,p2,p3,p4,p5,za,zb)
c      write(6,*)
c      pause
c--- DEBUG ONLY


!======= 3,4,5 symmetry 
      amp(1,2,2,1)=virt_trigam(p1,p2,p3,p4,p5,za,zb)
      amp(1,2,1,2)=virt_trigam(p1,p2,p3,p5,p4,za,zb)
      amp(1,1,2,2)=virt_trigam(p1,p2,p4,p5,p3,za,zb)
      
!======= line reversal 
      amp(2,2,2,1)=virt_trigam(p2,p1,p3,p4,p5,za,zb)
      amp(2,2,1,2)=virt_trigam(p2,p1,p3,p5,p4,za,zb)
      amp(2,1,2,2)=virt_trigam(p2,p1,p4,p5,p3,za,zb)

!======= conjugation 
      amp(2,1,1,2)=-virt_trigam(p1,p2,p3,p4,p5,zb,za)
      amp(2,1,2,1)=-virt_trigam(p1,p2,p3,p5,p4,zb,za)
      amp(2,2,1,1)=-virt_trigam(p1,p2,p4,p5,p3,zb,za)
!======= conjugation and line reversal 
      amp(1,1,1,2)=-virt_trigam(p2,p1,p3,p4,p5,zb,za)
      amp(1,1,2,1)=-virt_trigam(p2,p1,p3,p5,p4,zb,za)
      amp(1,2,1,1)=-virt_trigam(p2,p1,p4,p5,p3,zb,za)

      return 
      end
