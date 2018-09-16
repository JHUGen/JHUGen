      subroutine qqb_dm_monophot_g_Axamps(p,i1,i2,i3,i4,i5,i6,amp) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4)
      complex(dp):: amp(2,2,2,2,2)
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: Ax_monophot_helppC
      complex(dp):: Ax_monophot_helpmC

!=========== q(i1)+g(i2)+qb(i3)+gamma(i4)+x(i5)+x~(i6)
!---- order is quark helicity,gluon helicity photon helciity dm 
      if(xmass>1d-8) then 
!---------generate massless phase space 
         call gen_masslessvecs(p,q,i5,i6)
!--------- generate spinors 
         call spinoru(6,q,za,zb)
      else
!--------massless dm can use usual spinoru
         call spinoru(6,p,za,zb)       
      endif
      
      amp(:,:,:,:,:)=czip
     

!===== helicity conserving same sign gluon/photon 
      amp(2,2,2,2,1)=Ax_monophot_helppC(i1,i2,i3,i4,i5,i6,za,zb) 
      amp(2,2,2,1,2)=-Ax_monophot_helppC(i1,i2,i3,i4,i6,i5,za,zb) 
!===== line flip 
      amp(1,2,2,2,1)=Ax_monophot_helppC(i3,i2,i1,i4,i5,i6,za,zb)
      amp(1,2,2,1,2)=-Ax_monophot_helppC(i3,i2,i1,i4,i6,i5,za,zb)

!===== helicity violaiting same sign gluon/photon  = 0 

!====== helicity conserving opposite sign gluon/photon 
      amp(1,2,1,2,1)=Ax_monophot_helpmC(i1,i2,i3,i4,i5,i6,za,zb) 
      amp(1,2,1,1,2)=-Ax_monophot_helpmC(i1,i2,i3,i4,i6,i5,za,zb) 
!====== line flip 
      amp(2,2,1,2,1)=Ax_monophot_helpmC(i3,i2,i1,i4,i5,i6,za,zb)
      amp(2,2,1,1,2)=-Ax_monophot_helpmC(i3,i2,i1,i4,i6,i5,za,zb)

!======= Conjugate
      amp(1,1,1,1,2)=conjg(amp(2,2,2,2,1))
      amp(1,1,1,2,1)=conjg(amp(2,2,2,1,2)) 
      amp(2,1,1,1,2)=conjg(amp(1,2,2,2,1)) 
      amp(2,1,1,2,1)=conjg(amp(1,2,2,1,2))

      amp(1,1,2,2,1)=conjg(amp(2,2,1,1,2))
      amp(1,1,2,1,2)=conjg(amp(2,2,1,2,1))
      amp(2,1,2,2,1)=conjg(amp(1,2,1,1,2))
      amp(2,1,2,1,2)=conjg(amp(1,2,1,2,1))

      return 
      end 


      function Ax_monophot_helppC(i1,i2,i3,i4,i5,i6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Ax_monophot_helppC
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'dm_params.f' 
      integer:: i1,i2,i3,i4,i5,i6 
!===== helicity amplitude for q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (+,-) Current
      real(dp):: bp,beta,s56

      s56=real(za(i6,i5)*zb(i5,i6))
      beta=sqrt(1d0-4d0*xmass**2/s56) 
      bp=0.5d0*(one+beta)
 

      Ax_monophot_helppC=
     &((cone - cone*2d0*bp)*za(i1,i3)*za(i3,i6)**2*zb(i6,i5))/
     &   (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4))

      return 
      end 
      function Ax_monophot_helpmC(i1,i2,i3,i4,i5,i6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Ax_monophot_helpmC
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'dm_params.f' 
      integer:: i1,i2,i3,i4,i5,i6 
!===== helicity amplitude for q(i1)^+g(i2)^+qb(i3)^-+gamma(i4)^+ * Axial (+,+) Current
      real(dp):: bp,beta,s56
      integer:: i,j,k 
      real(dp):: st(mxpart,mxpart,mxpart) 

      s56=real(za(i6,i5)*zb(i5,i6))
      beta=sqrt(1d0-4d0*xmass**2/s56) 
      bp=0.5d0*(one+beta)

      st(:,:,:)=0d0 
      
      do i=1,6
         do j=1,6 
            do k=1,6
               st(i,j,k)=s(i,j)+s(j,k)+s(i,k) 
            enddo
         enddo
      enddo
      
!====== A terms 

      Ax_monophot_helpmC=  
     &-(((-cone + cone*2d0*bp)*(-(za(i1,i4)*za(i1,i6)*zb(i3,i1)) - 
     &        za(i1,i4)*za(i2,i6)*zb(i3,i2) + 
     &   za(i1,i2)*za(i4,i6)*zb(i3,i2) + za(i1,i4)*za(i4,i6)*zb(i4,i3))
     &       *(-(za(i1,i2)*zb(i5,i2)) - za(i1,i3)*zb(i5,i3)))/
     &    (za(i1,i2)*za(i2,i3)*zb(i4,i3)*
     &      (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2) + 
     &        za(i3,i4)*zb(i4,i3))))
      Ax_monophot_helpmC=Ax_monophot_helpmC+
     &  ((-cone + ctwo*bp)*za(i1,i6)*zb(i3,i2)*
     &    (za(i2,i4)*zb(i5,i2) + za(i3,i4)*zb(i5,i3)))/
     &  (st(i2,i3,i4)*za(i2,i3)*zb(i4,i3))

!==== Bterms 

      Ax_monophot_helpmC=Ax_monophot_helpmC+
     &((-cone + 2d0*cone*bp)*za(i1,i4)*(za(i1,i6)*zb(i2,i1) 
     &     - za(i4,i6)*zb(i4,i2))*
     &    zb(i5,i3))/(st(i1,i2,i4)*za(i1,i2)*zb(i4,i1))


      Ax_monophot_helpmC=Ax_monophot_helpmC
     &+((-cone + 2d0*cone*bp)*(-(za(i1,i2)*zb(i5,i2)) 
     &     - za(i1,i3)*zb(i5,i3))*
     &    (s(i2,i3)*za(i4,i6) + za(i2,i4)*za(i5,i6)*zb(i5,i2) + 
     &      za(i3,i4)*za(i5,i6)*zb(i5,i3)))/
     &  (za(i1,i2)*za(i2,i3)*zb(i4,i1)*
     &    (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2) 
     & + za(i3,i4)*zb(i4,i3)))

      return 
      end 
