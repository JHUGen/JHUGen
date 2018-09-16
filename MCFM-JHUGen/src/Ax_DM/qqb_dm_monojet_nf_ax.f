      subroutine qqb_dm_monojet_nf_ax(p,i1,i2,i3,i4,i5,qgqb)
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'first.f'
    
      real(dp):: p(mxpart,4),q(mxpart,4)
      integer:: i1,i2,i3,i4,i5
      real(dp):: qgqb(2)
      complex(dp):: amp_tree(2,2,2,2),amp(2,2,2,2) 
      integer:: h1,h2,h3,h4
      complex(dp):: L1
      real(dp):: bp,beta,s34,s123
      complex(dp):: sc123
!========= amplitudes for the axial part of q(i1)+g(i2)+qb(i3)+x(i4)+~x(i5) (triange loop) 
      integer:: j 
      real(dp):: nffac
      qgqb(:)=zip 
      amp(:,:,:,:)=czip
      if(xmass>1d-8) then 
!--------- generate massless phase space 
      call gen_masslessvecs(p,q,i4,i5)
!--------- generate spinors 
      call spinoru(5,q,za,zb)
      else
!-------- massless dm can use usual spinoru
         call spinoru(5,p,za,zb)          
      endif

!bp = 1/2(1+beta) 
! beta = sqrt(1-4xmass**2/s34) 
      s34=Dble(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34) 
      bp=0.5d0*(one+beta)

      sc123=cplx2(s(i1,i2)+s(i2,i3)+s(i3,i1),zip)
      s123=s(i1,i2)+s(i2,i3)+s(i3,i1)
!====== amplitudes : helicity conserving 
      amp(1,1,1,2)=(L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i4)*
     -     zb(i5,i3))/sc123 - 
     -  (2*bp*L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i4)*
     -     zb(i5,i3))/sc123  
      amp(1,2,1,2)= (L1(-s(i1,i3),-s123)*za(i1,i4)*zb(i3,i2)*
     -     zb(i5,i2))/sc123 - 
     -  (2*bp*L1(-s(i1,i3),-s123)*za(i1,i4)*zb(i3,i2)*
     -     zb(i5,i2))/sc123
      amp(2,1,1,2)= -((L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i4)*
     -       zb(i5,i1))/sc123)
     -   + (2*bp*L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i4)*
     -     zb(i5,i1))/sc123
      amp(2,2,1,2)= -((L1(-s(i1,i3),-s123)*za(i3,i4)*zb(i2,i1)*
     -       zb(i5,i2))/sc123)
     -   + (2*bp*L1(-s(i1,i3),-s123)*za(i3,i4)*zb(i2,i1)*
     -     zb(i5,i2))/sc123
      amp(1,1,2,1)= -((L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i5)*
     -       zb(i4,i3))/sc123)
     -   + (2*bp*L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i5)*
     -     zb(i4,i3))/sc123
      amp(2,1,2,1)=(L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i5)*
     -     zb(i4,i1))/sc123 - 
     -  (2*bp*L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i5)*
     -     zb(i4,i1))/sc123
      amp(1,2,2,1)= -((L1(-s(i1,i3),-s123)*za(i1,i5)*zb(i3,i2)*
     -       zb(i4,i2))/sc123)
     -   + (2*bp*L1(-s(i1,i3),-s123)*za(i1,i5)*zb(i3,i2)*
     -     zb(i4,i2))/sc123
      amp(2,2,2,1)=  (L1(-s(i1,i3),-s123)*za(i3,i5)*zb(i2,i1)*
     -     zb(i4,i2))/sc123 - 
     -  (2*bp*L1(-s(i1,i3),-s123)*za(i3,i5)*zb(i2,i1)*
     -     zb(i4,i2))/sc123

!======amplitudes : helicity violating 
      amp(1,1,1,1)=-((xmass*L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i4)*
     -       zb(i4,i3))/
     -     (sc123*zb(i5,i4)))
     -   - (xmass*L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i5)*
     -     zb(i5,i3))/
     -   (sc123*zb(i5,i4))
      amp(1,2,1,1)= -((xmass*L1(-s(i1,i3),-s123)*za(i1,i4)*zb(i3,i2)*
     -       zb(i4,i2))/
     -     (sc123*zb(i5,i4)))
     -   - (xmass*L1(-s(i1,i3),-s123)*za(i1,i5)*zb(i3,i2)*
     -     zb(i5,i2))/
     -   (sc123*zb(i5,i4))
      amp(2,1,1,1)=(xmass*L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i4)*
     -     zb(i4,i1))/
     -   (sc123*zb(i5,i4)) + 
     -  (xmass*L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i5)*
     -     zb(i5,i1))/
     -   (sc123*zb(i5,i4))
      amp(2,2,1,1)=(xmass*L1(-s(i1,i3),-s123)*za(i3,i4)*zb(i2,i1)*
     -     zb(i4,i2))/
     -   (sc123*zb(i5,i4)) + 
     -  (xmass*L1(-s(i1,i3),-s123)*za(i3,i5)*zb(i2,i1)*
     -     zb(i5,i2))/
     -   (sc123*zb(i5,i4))
      amp(1,1,2,2)=-((xmass*L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i4)*
     -       zb(i4,i3))/
     -     (sc123*za(i4,i5)))
     -   - (xmass*L1(-s(i1,i3),-s123)*za(i1,i2)*za(i2,i5)*
     -     zb(i5,i3))/
     -   (sc123*za(i4,i5))
      amp(2,1,2,2)= (xmass*L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i4)*
     -     zb(i4,i1))/
     -   (sc123*za(i4,i5)) + 
     -  (xmass*L1(-s(i1,i3),-s123)*za(i2,i3)*za(i2,i5)*
     -     zb(i5,i1))/
     -   (sc123*za(i4,i5))
      amp(2,2,2,2)=(xmass*L1(-s(i1,i3),-s123)*za(i3,i4)*zb(i2,i1)*
     -     zb(i4,i2))/
     -   (sc123*za(i4,i5)) + 
     -  (xmass*L1(-s(i1,i3),-s123)*za(i3,i5)*zb(i2,i1)*
     -     zb(i5,i2))/
     -   (sc123*za(i4,i5))

      call qqb_dm_monojet_Axamps(p,i1,i2,i3,i4,i5,amp_tree) 

!---- bulid amplitudes (sum over nf) 
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  qgqb(h1)=qgqb(h1)-
     & ason2pi/xn*Dble(conjg(amp_tree(h1,h2,h3,h4))*amp(h1,h2,h3,h4))
               enddo
            enddo
         enddo
      enddo

      nffac=zip
!---- now sum over active flavors 
      do j=1,nf 
        if(first) then 
           first=.false. 
           call check_dmAxC
        endif 
        nffac=nffac+0.5d0*(dmL(j)-dmR(j))
      enddo

      
      qgqb(1)=nffac*qgqb(1)
      qgqb(2)=nffac*qgqb(2)
      
      return 
      end 
