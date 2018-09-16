!--------------------------------------------------------------- 
! Subroutine for pp->Z+gamma + jet integrated frag dipoles-------------------
!--------------------------------------------------------------- 
!--- Author C. Williams Feb 2011 
!-----------------------------------------------------------------

!==== modified Dec 13 to remove call to outdated rescale_pjet and return_pjet 
!==== instead p_phys is now passed as in other fragmentaiton routines 

      subroutine qqb_zaj_fragdips(p,p_phys,qcd_tree,msq_out) 
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
      real(dp):: msq_qcd_s(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: virt_dips(3),xl(3),dot,fsq 
      real(dp):: aewo2pi,fi_gaq,ff_gaq
      external qcd_tree
      external qcd_tree_s
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp)::  mdum1(0:2,fn:nf,fn:nf) ! mqq
      real(dp)::  mdum2(0:2,fn:nf,fn:nf) ! mqq
      integer:: f
      aewo2pi=esq/(fourpi*twopi)      
      
      fsq=frag_scale**2
      
      xl(:)=0._dp 
      virt_dips(:)=0._dp
     
!---- Integrated dipoles are functions of p_gamma = z * pjet so need to rescale pjet

!==== below routine is outdated
!      call rescale_pjet(p) 
     
         xl(1)=log(-two*dot(p_phys,1,5)/fsq)
         virt_dips(1)=+aewo2pi*(fi_gaq(z_frag,p_phys,xl(1),5,1,2))
!
    
!---- Matrix elements conserve momenta thro pjet = sum of rest so return orignal pjet
!==== below routine is outdated
!     call return_pjet(p)
     
      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0._dp
            msq_out(j,k)=0._dp
         enddo
      enddo      

      
      
!------- msq_qcd should be qqb_z2jetx       
      call qcd_tree(p,msq_qcd,mdum2,msqx_cs,mdum1)
    
  
      do j=-nf,nf
         do k=-nf,nf
            
            
            
            if((k==0).and.(j.ne.0)) then 
               msq_out(j,k)=Q(abs(j))**2*msq_qcd(j,k)*virt_dips(1) 
            elseif((j==0).and.(k.ne.0)) then 
               msq_out(j,k)=Q(abs(k))**2*msq_qcd(j,k)*virt_dips(1)
            elseif((j==0).and.(k==0)) then 
               msq_out(j,k)=2._dp*Q(2)**2*virt_dips(1)*
     & (msqx_cs(0,j,k,2,-2)+msqx_cs(1,j,k,2,-2)+msqx_cs(2,j,k,2,-2))
     &             +3._dp*Q(1)**2*virt_dips(1)
     & *(msqx_cs(0,j,k,1,-1)+msqx_cs(1,j,k,1,-1)+msqx_cs(2,j,k,1,-1))
               msq_out(j,k)=2._dp*msq_out(j,k) 
              
              elseif (((j > 0).and.(k > 0)) .or. 
     &              ((j < 0).and.(k < 0))) then
                 msq_out(j,k)=Q(abs(j))**2*virt_dips(1)
     & *(msqx_cs(0,j,k,j,k)+msqx_cs(1,j,k,j,k)+msqx_cs(2,j,k,j,k))
                 msq_out(j,k)=msq_out(j,k)+Q(abs(k))**2*virt_dips(1)
     & *(msqx_cs(0,j,k,k,j)+msqx_cs(1,j,k,k,j)+msqx_cs(2,j,k,k,j))
              elseif((j>0).and.(k<0)) then 
                 if(j==-k) then 
                    do f=1,5
        msq_out(j,k)=msq_out(j,k)+2._dp*Q(f)**2*virt_dips(1)*
     &   (msqx_cs(0,j,k,f,-f)+msqx_cs(1,j,k,f,-f)+msqx_cs(2,j,k,f,-f))
                    enddo
                 else
           msq_out(j,k)=Q(abs(j))**2*virt_dips(1)
     &   *(msqx_cs(0,j,k,j,k)+msqx_cs(1,j,k,j,k)+msqx_cs(2,j,k,j,k))
           msq_out(j,k)=msq_out(j,k)+Q(abs(k))**2*virt_dips(1)
     &   *(msqx_cs(0,j,k,k,j)+msqx_cs(1,j,k,k,j)+msqx_cs(2,j,k,k,j))
                 endif
              elseif((j<0).and.(k>0)) then 
                 if(j==-k) then 
                    do f=1,5
        msq_out(j,k)=msq_out(j,k)+2._dp*Q(f)**2*virt_dips(1)*
     &   (msqx_cs(0,j,k,f,-f)+msqx_cs(1,j,k,f,-f)+msqx_cs(2,j,k,f,-f))
                    enddo
                 else
           msq_out(j,k)=Q(abs(j))**2*virt_dips(1)
     &   *(msqx_cs(0,j,k,j,k)+msqx_cs(1,j,k,j,k)+msqx_cs(2,j,k,j,k))
           msq_out(j,k)=msq_out(j,k)+Q(abs(k))**2*virt_dips(1)
     &   *(msqx_cs(0,j,k,k,j)+msqx_cs(1,j,k,k,j)+msqx_cs(2,j,k,k,j))
                 endif
              endif

        enddo
      enddo
     
     
      return 
      end subroutine


        
