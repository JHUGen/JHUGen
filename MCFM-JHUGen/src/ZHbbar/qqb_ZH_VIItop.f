      function qqb_ZH_VIItop(i1,i2,i3,i4,p,j) 
!===== C. Williams Sept 2015 
!===== routine which calculates the VII style contributions to 
!==== qqb => ZH (see e.g. Brein, Harlander Wisemann and Zirke) 
      implicit none 
      include 'types.f' 
      real(dp) :: qqb_ZH_VIItop 
      include 'constants.f' 
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f' 
      include 'qcdcouple.f'
      include 'nproc.f'
      include 'nf.f'
      include 'zcouple.f' 
      include 'cutoff.f'
      integer:: i1,i2,i3,i4,j
      integer h1,h2
      real(dp) :: p(mxpart,4),cl(2),cq(2),cutoff_orig
      complex(dp) :: qqb_ZH_asymt,asympt(2,2)
      complex(dp) :: qqb_WH_treeamp,tree(2,2) 
      complex(dp) :: fac_asym,fac_tree
      complex(dp) :: prop_12,prop_34,Ch
      real(dp) :: ccou(2,2),s134,s234,corr
      
      include "cplx.h"
      
      qqb_ZH_VIItop=zip
      

!==== debug, routine is written for leptons, needs generalization... 
      cl(1)=l1
      cl(2)=r1
!      if((nproc.ne.101)
!     &     .or.(nproc.ne.104).or.(nproc.ne.106).or.(nproc.ne.109)) then 
!         write(6,*) 'Warning, only leptonic decays implemented in '
!         write(6,*) 'file qqb_ZH_VIItop.f' 
!         write(6,*) 'nproc = 'nproc,' not yet supported' 
!         stop 
!      endif

      cq(1)=-l(j)
      cq(2)=r(j)


!==== factors 
      fac_asym=-(ason4pi)**2*Gf*wmass*(2._dp*cf)*sqrt(esq)*(l(4)-r(4))*8._dp
      fac_tree=esq*wmass*gw/(one-xw)*2._dp

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
 
!======== cal
!==== fill tree array and assympt 
      
      tree(1,1)=qqb_WH_treeamp(i2,i1,i3,i4,za,zb)
      tree(2,1)=qqb_WH_treeamp(i1,i2,i3,i4,za,zb)
      tree(1,2)=qqb_WH_treeamp(i2,i1,i4,i3,za,zb)
      tree(2,2)=qqb_WH_treeamp(i1,i2,i4,i3,za,zb)

      asympt(1,1)=qqb_ZH_asymt(i2,i1,i3,i4,za,zb) 
      asympt(2,1)=qqb_ZH_asymt(i1,i2,i3,i4,za,zb) 
      asympt(1,2)=qqb_ZH_asymt(i2,i1,i4,i3,za,zb) 
      asympt(2,2)=qqb_ZH_asymt(i1,i2,i4,i3,za,zb) 


      asympt(:,:)=asympt(:,:)*prop_34*fac_asym
      tree(:,:)=tree(:,:)*prop_12*prop_34*fac_tree

!===== chiral couplings 
      qqb_ZH_VIItop=zip
      do h1=1,2
         do h2=1,2 
            ccou(h1,h2)=cq(h1)*cl(h2)**2
            qqb_ZH_VIItop=qqb_ZH_VIItop+ 
     &    xn*ccou(h1,h2)*real((conjg(tree(h1,h2))*asympt(h1,h2)
     &           +conjg(asympt(h1,h2))*tree(h1,h2)),dp)

         enddo
      enddo

      return 

 999  continue
      cutoff=cutoff_orig
      return 

      end

      
      function qqb_ZH_asymt(i1,i2,i3,i4,za,zb) 
!==== amplitude for intf. with EFT loop 
!     q(i1)^++qb(i2)^-=>ell(i3)^-+ell(i4)^+
      implicit none 
      include 'types.f' 
      complex(kind=dp)  qqb_ZH_asymt
      include 'constants.f'  
      include 'mxpart.f' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer i1,i2,i3,i4 
      

      qqb_ZH_asymt=zb(i1,i4)*za(i3,i2)/(s(i3,i4))

      return 
      end 

