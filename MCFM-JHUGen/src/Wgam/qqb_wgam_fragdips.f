!--------------------------------------------------------------- 
!   This subroutine checks the number of external dipoles----
!---absorbing the correct number into the fragmenation functions 
!   it then returns the finite (msq_qcd*dip) ---------------------
!--------------------------------------------------------------- 
!--- Author C. Williams Dec 2010 
!-----------------------------------------------------------------

!-----As a first pass we drop the poles and return qcd*dips
!-----for proccess where u -> d (charge flip) 

      subroutine qqb_wgam_fragdips(p,p_phys,qcd_tree,msq_out) 
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

      xl=log(-two*dot(p_phys,2,5)/fsq)
      virt_dips=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,5,2,2))
      

      do j=-nf,nf
         do k=-nf,nf
            msq_out(j,k)=0._dp
         enddo
      enddo
      
      call qcd_tree(p,msq_qcd) 

      do j=-nf,nf
         do k=-nf,nf
            
            if((j==0).and.(k.ne.0)) then
               if(mod(abs(k),2)==1) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(2)**2*virt_dips
               else
                  msq_out(j,k)=msq_qcd(j,k)*Q(1)**2*virt_dips
               endif
            elseif((j.ne.0).and.(k==0)) then
               if(mod(abs(j),2)==1) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(2)**2*virt_dips
               else
                  msq_out(j,k)=msq_qcd(j,k)*Q(1)**2*virt_dips
               endif
            endif
            
         enddo
      enddo
     

      return 
      end subroutine


        
