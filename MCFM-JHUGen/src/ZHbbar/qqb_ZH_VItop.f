!===== C. Williams Sept 2015 
!===== routine which calcualtes the process 
!====   q(i1)^++qb(i2)^-=>ell(i3)^-+ell(i4)^+ 
!     + Higgs where superscripts denote helicity. 
!===== this is the piece which goes like the Higgs-top coupling in the EFT (i.e. O(alpha_s^2) in QCD)

!===== Note Higgs decay not handeled in this routine 
!===== routine returns msq with no color averaging  
      function qqb_ZH_VItop(i1,i2,i3,i4,p,j) 
      implicit none 
      include 'types.f' 
      include 'constants.f' 
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f' 
      include 'qcdcouple.f'
      include 'nf.f'
      include 'zcouple.f'
      include 'cutoff.f'
      real(dp)::qqb_ZH_VItop
      integer:: i1,i2,i3,i4,j
      real(dp) :: p(mxpart,4),cutoff_orig
      complex(dp) :: qqb_WH_HtopEFT_loop
      complex(dp) :: qqb_WH_treeamp
      complex(dp) :: fac_eft,fac_tree
      complex(dp) :: prop_12,prop_34,Ch
      complex(dp) :: loops(2,2),tree(2,2) 
      integer h1,h2
      real(dp):: cl(2),cq(2)
      include "cplx.h"

      qqb_ZH_VItop=zip     
      cl(1)=l1
      cl(2)=r1
      cq(1)=l(j) 
      cq(2)=r(j)
      
      Ch=as/(3._dp*pi*sqrt(vevsq))

!==== factors 
      fac_eft=-ason4pi*(2._dp*cf)*Ch*esq*2._dp
  
      fac_tree=gw/(one-xw)*wmass*esq*2._dp

    
!==== just need spinors for 1-4  
      call spinoru(4,p,za,zb) 

!=====this routine runs in NNLO mode, and can have very small
!=====values of s, since its non-singular lets ensure stability
!==== by doing a smalls cut.
      cutoff_orig=cutoff
      cutoff=1.E-3_dp
      call smalls(s,4,*999)
      cutoff=cutoff_orig

!======= props 
      prop_34=s(i3,i4)/cplx2(s(i3,i4)-zmass**2,zmass*zwidth)
      prop_12=s(i1,i2)/cplx2(s(i1,i2)-zmass**2,zmass*zwidth)
 
!======== call amplitude for the loop part 
!      EFT_loop=qqb_WH_HtopEFT_loop(i1,i2,i3,i4,za,zb) 


      tree(1,1)=qqb_WH_treeamp(i2,i1,i3,i4,za,zb) 
      tree(2,1)=qqb_WH_treeamp(i1,i2,i3,i4,za,zb) 
      tree(1,2)=qqb_WH_treeamp(i2,i1,i4,i3,za,zb) 
      tree(2,2)=qqb_WH_treeamp(i1,i2,i4,i3,za,zb) 

      loops(1,1)=qqb_WH_HtopEFT_loop(i2,i1,i3,i4,za,zb) 
      loops(2,1)=qqb_WH_HtopEFT_loop(i1,i2,i3,i4,za,zb) 
      loops(1,2)=qqb_WH_HtopEFT_loop(i2,i1,i4,i3,za,zb) 
      loops(2,2)=qqb_WH_HtopEFT_loop(i1,i2,i4,i3,za,zb) 

!======= call amplitude for the LO intf. 
!      tree=qqb_WH_treeamp(i1,i2,i3,i4,za,zb) 

!===== dress with factors and propagators 
      tree(:,:)=tree(:,:)*fac_tree*prop_12*prop_34
      
      loops(:,:)=loops(:,:)*fac_EFT*prop_34
!      Loops(:,:)=tree(:,:)
! ====== perform interference and include overall color charge

      qqb_ZH_VItop=zip
      do h1=1,2
         do h2=1,2
            qqb_ZH_VItop=qqb_ZH_VItop+
     &       xn*(cq(h1)*cl(h2))**2*
     &       real(conjg(tree(h1,h2))*loops(h1,h2)
     &           +conjg(loops(h1,h2))*tree(h1,h2),dp)

         enddo
      enddo
!     qqb_WH_HtopEFT=xn*(conjg(tree)*EFT_loop+conjg(EFT_loop)*tree)

      return 

 999  continue
      cutoff=cutoff_orig
      return 

      end 

   
