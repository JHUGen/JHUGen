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
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'lastphot.f'
      double precision p(mxpart,4),p_phys(mxpart,4)
      double precision msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer j,k
      double precision virt_dip,xl,dot,fsq 
      double precision aewo2pi,fi_gaq
      external qcd_tree

      aewo2pi=esq/(fourpi*twopi)            
      fsq=frag_scale**2

c--- assemble integrated subtraction terms     
      xl=dlog(-two*dot(p_phys,1,lastphot)/fsq)
      virt_dip=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,lastphot,1,2))

c--- initialize array      
      msq_out(:,:)=0d0
      
c--- fill underlying QCD matrix elements     
      call qcd_tree(p,msq_qcd) 

c--- fill output array
      do j=-nf,nf
      do k=-nf,nf
            
!   factor of three cancelled by statistical factor because three photons
         if    ((j.eq.0).and.(k.ne.0)) then
            msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dip
         elseif((j.ne.0).and.(k.eq.0)) then
            msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dip
         endif
         
      enddo
      enddo
     
      return 
      end


        
