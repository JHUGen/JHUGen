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

      subroutine qqb_gmgmjt_frag_combo(p,p_phys,msq_out) 
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
      include 'msqbits.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      real(dp):: msq(-nf:nf,-nf:nf),msq_qcd_swap(-nf:nf,-nf:nf) 
      real(dp):: msqbits_qcd(12),msqbits_qcd_swap(12)
      integer:: j,k,i
      real(dp):: virt_dip,xl,dot,fsq 
      real(dp):: D(0:5),d_sum_q
      real(dp):: aewo2pi,fi_gaq,facgg
      real(dp):: facDgg
      external qcd_tree
      common/D/D
!$omp threadprivate(/D/)

      aewo2pi=esq/(fourpi*twopi)            
      fsq=frag_scale**2
      msq_out(:,:)=0._dp

c--- assemble integrated subtraction terms     
      xl=log(-two*dot(p_phys,1,lastphot)/fsq)
      virt_dip=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,lastphot,1,2))


!======== fragmentation functions 
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0._dp
         if     (fragset == 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
         elseif (fragset == 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset == 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0) 
         elseif (fragset == 'GdRG_NLO') then 
            call GGdR_frag(z_frag,i,D(i),1) 
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop        
         endif
      enddo

!==== Debug turn off frag funcs for comp to Johns routine 
!      D(:)=0._dp 

      



      call qqb_dirgam_g(p,msq_qcd)
      msqbits_qcd(:)=msqbits(:)
!====== swapped routine 
      call qqb_dirgam_g_swap(p,msq_qcd_swap) 
      msqbits_qcd_swap(:)=msqbits(:)

    
c--- this is the factor to apply the right couplings
c--- in the gluon-gluon contribution
      facgg=
     &  (2._dp*Q(2)**4+real(nf-2,dp)*Q(1)**4)
     & /(2._dp*Q(2)**2+real(nf-2,dp)*Q(1)**2)
      
      facDgg=((D(1)+D(3)+D(5))*Q(1)**2+(D(2)+D(4))*Q(2)**2)
     & /(2._dp*Q(2)**2+real(nf-2,dp)*Q(1)**2)

    
c--- fill output array
      do j=-nf,nf
      do k=-nf,nf
            
     
         if    ((j==0).and.(k.ne.0)) then
            msq_out(j,k)=msq_qcd(j,k)*(Q(k)**2*virt_dip+D(abs(k)))            
         elseif((j.ne.0).and.(k==0)) then
            msq_out(j,k)=msq_qcd(j,k)*(Q(j)**2*virt_dip+D(abs(j)))
         elseif((j==0).and.(k==0)) then
            msq_out(j,k)=msq_qcd(j,k)*2._dp*(facgg*virt_dip+facDgg)
         elseif((j>0).and.(k>0).and.(j.ne.-k)) then
            msq_out(j,k)=msq_qcd(j,k)*(Q(j)**2*virt_dip+D(j))
     &           +msq_qcd_swap(j,k)*(Q(k)**2*virt_dip+D(k))
         elseif((j<0).and.(k<0).and.(j.ne.-k)) then
            msq_out(j,k)=msq_qcd(j,k)*(Q(j)**2*virt_dip+D(-j))
     &           +msq_qcd_swap(j,k)*(Q(k)**2*virt_dip+D(-k))
         elseif((j>0).and.(k<0).and.(j.ne.-k)) then 
            msq_out(j,k)=msq_qcd(j,k)*(Q(j)**2*virt_dip+D(j))
     &           +msq_qcd_swap(j,k)*(Q(k)**2*virt_dip+D(-k))
         elseif((j<0).and.(k>0).and.(j.ne.-k)) then 
            msq_out(j,k)=msq_qcd_swap(j,k)*(Q(j)**2*virt_dip+D(-j))
     &           +msq_qcd(j,k)*(Q(k)**2*virt_dip+D(-k))
         elseif((j==-k).and.(j.ne.0).and.(k.ne.0)) then 
!=========== u type 
            if(abs(j)==2) then
               msq_out(j,k)=
     &               +(Q(2)**2*virt_dip+D(2))*msqbits_qcd(uub_uub)
     &                +(Q(2)**2*virt_dip+D(4))*msqbits_qcd(uub_ccb)
     &      +(3._dp*Q(1)**2*virt_dip+D(1)+D(3)+D(5))*msqbits_qcd(uub_ddb)
     &               +(Q(2)**2*virt_dip+D(2))*msqbits_qcd_swap(uub_uub)
     &                +(Q(2)**2*virt_dip+D(4))*msqbits_qcd_swap(uub_ccb)
     & +(3._dp*Q(1)**2*virt_dip+D(1)+D(3)+D(5))*msqbits_qcd_swap(uub_ddb)
!=========c type 
            elseif(abs(j)==4) then
               msq_out(j,k)=
     &               +(Q(2)**2*virt_dip+D(4))*msqbits_qcd(uub_uub)
     &                +(Q(2)**2*virt_dip+D(2))*msqbits_qcd(uub_ccb)
     &      +(3._dp*Q(1)**2*virt_dip+D(1)+D(3)+D(5))*msqbits_qcd(uub_ddb)
     &               +(Q(2)**2*virt_dip+D(4))*msqbits_qcd_swap(uub_uub)
     &                +(Q(2)**2*virt_dip+D(2))*msqbits_qcd_swap(uub_ccb)
     & +(3._dp*Q(1)**2*virt_dip+D(1)+D(3)+D(5))*msqbits_qcd_swap(uub_ddb)
  
!==========d type             
            elseif(abs(j)==1) then
              msq_out(j,k)=
     &               +(Q(1)**2*virt_dip+D(1))*msqbits_qcd(ddb_ddb)
     &           +(2._dp*Q(1)**2*virt_dip+D(3)+D(5))*msqbits_qcd(ddb_ssb)
     &         +(2._dp*Q(2)**2*virt_dip+D(2)+D(4))*msqbits_qcd(ddb_uub)
     &               +(Q(1)**2*virt_dip+D(1))*msqbits_qcd_swap(ddb_ddb)
     &      +(2._dp*Q(1)**2*virt_dip+D(3)+D(5))*msqbits_qcd_swap(ddb_ssb)
     &      +(2._dp*Q(2)**2*virt_dip+D(2)+D(4))*msqbits_qcd_swap(ddb_uub)
    

!====s type
            elseif(abs(j)==3) then
              msq_out(j,k)=
     &               +(Q(1)**2*virt_dip+D(3))*msqbits_qcd(ddb_ddb)
     &            +(2._dp*Q(1)**2*virt_dip+D(1)+D(5))*msqbits_qcd(ddb_ssb)
     &         +(2._dp*Q(2)**2*virt_dip+D(2)+D(4))*msqbits_qcd(ddb_uub)
     &               +(Q(1)**2*virt_dip+D(3))*msqbits_qcd_swap(ddb_ddb)
     &      +(2._dp*Q(1)**2*virt_dip+D(1)+D(5))*msqbits_qcd_swap(ddb_ssb)
     &      +(2._dp*Q(2)**2*virt_dip+D(2)+D(4))*msqbits_qcd_swap(ddb_uub)
    
!============b type 
                   elseif(abs(j)==5) then
              msq_out(j,k)=
     &               +(Q(1)**2*virt_dip+D(5))*msqbits_qcd(ddb_ddb)
     &            +(2._dp*Q(1)**2*virt_dip+D(3)+D(1))*msqbits_qcd(ddb_ssb)
     &         +(2._dp*Q(2)**2*virt_dip+D(2)+D(4))*msqbits_qcd(ddb_uub)
     &             +(Q(1)**2*virt_dip+D(5))*msqbits_qcd_swap(ddb_ddb)
     &       +(2._dp*Q(1)**2*virt_dip+D(3)+D(1))*msqbits_qcd_swap(ddb_ssb)
     &     +(2._dp*Q(2)**2*virt_dip+D(2)+D(4))*msqbits_qcd_swap(ddb_uub)
              
            endif

!======== higher order gluon corrections 
         elseif((j==0).and.(k.ne.0)) then 
            msq_out(j,k)=msq_out(j,k)+D(0)*msq_qcd_swap(j,k)
         elseif((j.ne.0).and.(k==0)) then
            msq_out(j,k)=msq_out(j,k)+D(0)*msq_qcd_swap(j,k)
         elseif((j==0).and.(k==0)) then 
            msq_out(j,k)=msq_out(j,k)
     &           +D(0)*(msq_qcd(j,k)+msq_qcd_swap(j,k)) 
         endif
         
 11      continue
      enddo
      enddo
     
      
      return 
      end


        
