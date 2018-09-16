!--------------------------------------------------------------- 
!   This subroutine checks the number of external dipoles----
!---absorbing the correct number into the fragmenation functions 
!   it then returns the finite (msq_qcd*dip) ---------------------
!--------------------------------------------------------------- 
!--- Author C. Williams April 2012 
!-----------------------------------------------------------------
!---- modified March 2013 during frag re-write


      subroutine qqb_zaa_fragdips(p,p_phys,qcd_tree,msq_out) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: virt_dips,xl,dot,fsq 
      real(dp):: aewo2pi,fi_gaq
      external qcd_tree
            

      aewo2pi=esq/(fourpi*twopi)      
      
      fsq=frag_scale**2
     
      
     
      xl=log(-two*dot(p_phys,2,6)/fsq)
      virt_dips=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,6,2,2))
      

      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0._dp
            msq_out(j,k)=0._dp
         enddo
      enddo

  
      
      call qcd_tree(p,msq_qcd) 


      do j=-nf,nf
         do k=-nf,nf
            
!   factor of two cancelled by statistical factor because two photons
            if((j==0).and.(k.ne.0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dips
            elseif((j.ne.0).and.(k==0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dips                    
            endif
            
         enddo
      enddo
     

      return 
      end subroutine


        
