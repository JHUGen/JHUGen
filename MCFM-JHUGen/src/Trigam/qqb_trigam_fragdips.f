!--------------------------------------------------------------- 
!   This subroutine checks the number of external dipoles----
!---absorbing the correct number into the fragmenation functions 
!   it then returns the finite (msq_qcd*dip) ---------------------
!--------------------------------------------------------------- 
!--- Author C. Williams Mar 2013 
!-----------------------------------------------------------------

c--- Passed in
c---   p:      array of momenta to evaluate matrix elements
c---   p_phys: array of momenta to evaluate integrated dipoles

      subroutine qqb_trigam_fragdips(p,p_phys,qcd_tree,msq_out) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'lastphot.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: virt_dip,xl,dot,fsq 
      real(dp):: aewo2pi,fi_gaq
      external qcd_tree

      aewo2pi=esq/(fourpi*twopi)            
      fsq=frag_scale**2

c--- assemble integrated subtraction terms     
      xl=log(-two*dot(p_phys,1,lastphot)/fsq)
      virt_dip=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,lastphot,1,2))

c--- initialize array      
      msq_out(:,:)=0._dp
      
c--- fill underlying QCD matrix elements     
      call qcd_tree(p,msq_qcd) 

c--- fill output array
      do j=-nf,nf
      do k=-nf,nf
            
!   factor of three cancelled by statistical factor because three photons
         if    ((j==0).and.(k.ne.0)) then
            msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dip
         elseif((j.ne.0).and.(k==0)) then
            msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dip
         endif
         
      enddo
      enddo
     
      return 
      end


        
