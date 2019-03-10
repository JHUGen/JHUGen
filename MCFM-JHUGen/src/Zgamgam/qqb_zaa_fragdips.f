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
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      double precision p(mxpart,4),p_phys(mxpart,4)
      double precision msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer j,k
      double precision virt_dips,xl,dot,fsq 
      double precision aewo2pi,fi_gaq
      external qcd_tree
            

      aewo2pi=esq/(fourpi*twopi)      
      
      fsq=frag_scale**2
     
      
     
      xl=dlog(-two*dot(p_phys,2,6)/fsq)
      virt_dips=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,6,2,2))
      

      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0d0
            msq_out(j,k)=0d0
         enddo
      enddo

  
      
      call qcd_tree(p,msq_qcd) 


      do j=-nf,nf
         do k=-nf,nf
            
!   factor of two cancelled by statistical factor because two photons
            if((j.eq.0).and.(k.ne.0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dips
            elseif((j.ne.0).and.(k.eq.0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dips                    
            endif
            
         enddo
      enddo
     

      return 
      end subroutine


        
