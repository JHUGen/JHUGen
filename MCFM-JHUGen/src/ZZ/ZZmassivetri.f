      subroutine ZZmassivetri(k1,k2,k3,k4,k5,k6,za,zb,mt,totrat,drat,
     & tri)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'ggZZintegrals.f'
      include 'ZZclabels.f'
      include 'docheck.f'
C      include 'first.f'
C---Indices of coefficients are CLL(integral,h1,h2,h3,h5)
      complex(dp):: CLL(6,2,2,2,2)
c      complex(dp):: CLR(6,2,2,2,2),CRL(6,2,2,2,2),CRR(6,2,2,2,2)
C---Indices of xpp etc are xpp(h3,h5) 
      complex(dp):: xapp(2,2),xamp(2,2),xapm(2,2),xamm(2,2),
     & Xrat(2,2,2,2),crat(2,2,2,2,6),drat(2,2,2,2,3),
     & totrat(2,2,2,2),tri12_34rat(2,2,2,2),tri(2,2,2,2,-2:0)
      real(dp):: mt
      integer:: k1,k2,k3,k4,k5,k6,h1,h2,h3,h5,j

c--- QCDLoop already initialized from call to massivebox
c      if (first) then
c      first=.false. 
c      write(6,*) 'mt',mt
c      write(6,*) 'musq',musq
c      call qlinit
c      pause
c      endif

c--- Initialize all triangle coefficients to zero
      CLL(:,:,:,:,:)=czip
c      CLR(:,:,:,:,:)=czip
c      CRL(:,:,:,:,:)=czip
c      CRR(:,:,:,:,:)=czip

C----12_34
c      call ZZtri12_34LR(k1,k2,k3,k4,k5,k6,za,zb,mt,xapp,xamp)
c      CLR(c12_34,2,2,:,:)=xapp(:,:)
c      CLR(c12_34,1,2,:,:)=xamp(:,:)
C----1_2
c      call ZZtri1_2LR(k1,k2,k3,k4,k5,k6,za,zb,mt,xapp,xamp)
c      CLR(c1_2,2,2,:,:)=xapp(:,:)
c      CLR(c1_2,1,2,:,:)=xamp(:,:)
C----1_34
c      call ZZtri1_34LR(k1,k2,k3,k4,k5,k6,za,zb,mt,xapp,xamp)
c      write(6,*) '1_34:xapp(1,1)',xapp(1,1)
c      write(6,*) '1_34:xamp(1,1)',xamp(1,1)
c      CLR(c1_34,2,2,:,:)=xapp(:,:)
c      CLR(c1_34,1,2,:,:)=xamp(:,:)
C----1_56
c      call ZZtri1_34LR(k1,k2,k5,k6,k3,k4,za,zb,mt,xapp,xamp)
c      write(6,*) '1_56:xapp(1,1)',xapp(1,1)
c      write(6,*) '1_56:xamp(1,1)',xamp(1,1)
c      CLR(c1_56,2,2,:,:)=xapp(:,:)
c      CLR(c1_56,1,2,:,:)=xamp(:,:)

C----2_34
c      call ZZtri1_34LR(k2,k1,k4,k3,k6,k5,zb,za,mt,xamm,xamp)
c      write(6,*) '2_34:xamm(1,1)',xamm(1,1)
c      write(6,*) '2_34:xamp(1,1)',xamp(1,1)
c      CLR(c1_56,1,1,:,:)=xamm(:,:)
c      CLR(c1_56,1,2,:,:)=xamp(:,:)

C----2_56
c      call ZZtri1_34LR(k2,k1,k6,k5,k4,k3,zb,za,mt,xamm,xamp)
c      write(6,*) '2_56:xamm(1,1)',xamm(1,1)
c      write(6,*) '2_56:xamp(1,1)',xamp(1,1)
c      CLR(c1_56,2,2,:,:)=xamm(:,:)
c      CLR(c1_56,1,2,:,:)=xamp(:,:)



c--- triangle (1,2)     
      call ZZtri1_2LL(k1,k2,k3,k4,k5,k6,za,zb,mt,xapp,xamp,xapm,xamm,
     & Xrat)
      CLL(c1_2,1,1,:,:)=xamm(:,:)
      CLL(c1_2,1,2,:,:)=xamp(:,:)
      CLL(c1_2,2,1,:,:)=xapm(:,:)
      CLL(c1_2,2,2,:,:)=xapp(:,:)
      crat(:,:,:,:,c1_2)=Xrat(:,:,:,:)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
c      write(6,*) 'Xrat(1,2,1,1)',Xrat(1,2,1,1)

c--- triangle (1,34)      
      call ZZtri1_34LL(k1,k2,k3,k4,k5,k6,za,zb,mt,xapp,xamp,xapm,xamm,
     & Xrat)
c      write(6,*) 'Xapp(1,1)',Xapp(1,1)
c      write(6,*) 'Xamp(1,1)',Xamp(1,1)
      CLL(c1_34,1,1,:,:)=xamm(:,:)
      CLL(c1_34,1,2,:,:)=xamp(:,:)
      CLL(c1_34,2,1,:,:)=xapm(:,:)
      CLL(c1_34,2,2,:,:)=xapp(:,:)
      crat(:,:,:,:,c1_34)=Xrat(:,:,:,:)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
c      write(6,*) 'Xrat(1,2,1,1)',Xrat(1,2,1,1)
      
c--- triangle (1,56)      
      call ZZtri1_34LL(k1,k2,k5,k6,k3,k4,za,zb,mt,xapp,xamp,xapm,xamm,
     & Xrat)
c      write(6,*) 'Xapp(1,1)',Xapp(1,1)
c      write(6,*) 'Xamp(1,1)',Xamp(1,1)
c--- since 34 and 56 have been interchanged in call, swap h3 and h5 elements
      do h3=1,2
      do h5=1,2
      CLL(c1_56,1,1,h3,h5)=xamm(h5,h3)
      CLL(c1_56,1,2,h3,h5)=xamp(h5,h3)
      CLL(c1_56,2,1,h3,h5)=xapm(h5,h3)
      CLL(c1_56,2,2,h3,h5)=xapp(h5,h3)
      crat(:,:,h3,h5,c1_56)=Xrat(:,:,h5,h3)
      enddo
      enddo
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
c      write(6,*) 'Xrat(1,2,1,1)',Xrat(1,2,1,1)
      
c--- triangle (2,34) - note interchange of xapm and xamp since 1 and 2 swapped
      call ZZtri1_34LL(k2,k1,k3,k4,k5,k6,za,zb,mt,xapp,xapm,xamp,xamm,
     & Xrat)
c      write(6,*) 'Xapp(1,1)',Xapp(1,1)
c      write(6,*) 'Xamp(1,1)',Xamp(1,1)
      CLL(c2_34,1,1,:,:)=xamm(:,:)
      CLL(c2_34,1,2,:,:)=xamp(:,:)
      CLL(c2_34,2,1,:,:)=xapm(:,:)
      CLL(c2_34,2,2,:,:)=xapp(:,:)
      do h1=1,2
      do h2=1,2
      crat(h1,h2,:,:,c2_34)=Xrat(h2,h1,:,:) ! due to interchange of 1 and 2
      enddo
      enddo
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
c      write(6,*) 'Xrat(1,2,1,1)',Xrat(1,2,1,1)

c--- triangle (2,56) - note interchange of xapm and xamp since 1 and 2 swapped    
      call ZZtri1_34LL(k2,k1,k5,k6,k3,k4,za,zb,mt,xapp,xapm,xamp,xamm,
     & Xrat)
c      write(6,*) 'Xapp(1,1)',Xapp(1,1)
c      write(6,*) 'Xamp(1,1)',Xamp(1,1)
c--- since 34 and 56 have been interchanged in call, swap h3 and h5 elements
      do h3=1,2
      do h5=1,2
      CLL(c2_56,1,1,h3,h5)=xamm(h5,h3)
      CLL(c2_56,1,2,h3,h5)=xamp(h5,h3)
      CLL(c2_56,2,1,h3,h5)=xapm(h5,h3)
      CLL(c2_56,2,2,h3,h5)=xapp(h5,h3)
      do h1=1,2
      do h2=1,2
      crat(h1,h2,h3,h5,c2_56)=Xrat(h2,h1,h5,h3) ! due to interchange of 1 and 2 as well
      enddo
      enddo
      enddo
      enddo
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
c      write(6,*) 'Xrat(1,2,1,1)',Xrat(1,2,1,1)


c--- compute rational contribution of tri12_34 from known total rational
c--- piece and contributions of boxes and other triangles
      tri12_34rat(:,:,:,:)=totrat(:,:,:,:)
     & -(drat(:,:,:,:,1)+drat(:,:,:,:,2)+drat(:,:,:,:,3))*(-1._dp/6._dp)
     & -(crat(:,:,:,:,c1_34)+crat(:,:,:,:,c1_56)+crat(:,:,:,:,c2_34)
     &  +crat(:,:,:,:,c2_56)+crat(:,:,:,:,c1_2))*(-0.5_dp)
c--- NB: coefficient of mt^2 is -(-2)=2 times this number

c      write(6,*) 'Reconstr,1,1',two*tri12_34rat(1,2,1,1)
c      write(6,*) 'Reconstr,1,2',two*tri12_34rat(1,2,1,2)
c      write(6,*) 'Reconstr,2,1',two*tri12_34rat(1,2,2,1)
c      write(6,*) 'Reconstr,2,2',two*tri12_34rat(1,2,2,2)

c--- result in this file is only correct for mt=0      
      call ZZtri12_34LL(k1,k2,k3,k4,k5,k6,za,zb,mt,xapp,xamp,xapm,xamm)
c      pause
c--- transfer and add mt^2 coefficient calculated above      
      CLL(c12_34,2,2,:,:)=xapp(:,:)+mt**2*two*tri12_34rat(2,2,:,:)
      CLL(c12_34,1,2,:,:)=xamp(:,:)+mt**2*two*tri12_34rat(1,2,:,:)
      CLL(c12_34,2,1,:,:)=xapm(:,:)+mt**2*two*tri12_34rat(2,1,:,:)
      CLL(c12_34,1,1,:,:)=xamm(:,:)+mt**2*two*tri12_34rat(1,1,:,:)
c      write(6,*) '12_34LL:xapp(1,1)',xapp(1,1)
c      write(6,*) '12_34LL:xamp(1,1)',xamp(1,1)
c      write(6,*) '12_34LL:xamp(1,2)',xamp(1,2)
c      write(6,*) '12_34LL:xamp(2,1)',xamp(2,1)
c      write(6,*) '12_34LL:xamp(2,2)',xamp(2,2)
 
c      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
c      s156=s(k1,k5)+s(k1,k6)+s(k5,k6)
c      s34=s(k3,k4)
c      s56=s(k5,k6)
c      write(6,*) 'CLL(c1_2,1,1,1,1)/s(k1,k2)',
c     & CLL(c1_2,1,1,1,1)/s(k1,k2),C0(c1_2)*s(k1,k2)
c      write(6,*) 'CLL(c1_34,1,1,1,1)/(s134-s34)',
c     & CLL(c1_34,1,1,1,1)/(s134-s34),C0(c1_34)*(s134-s34)
c      write(6,*) 'CLL(c2_34,1,1,1,1)/(s234-s34)',
c     & CLL(c2_34,1,1,1,1)/(s156-s34),C0(c2_34)*(s156-s34)
c      write(6,*) 'CLL(c1_56,1,1,1,1)/(s156-s56)',
c     & CLL(c1_56,1,1,1,1)/(s156-s56),C0(c1_56)*(s156-s56)
c      write(6,*) 'CLL(c2_56,1,1,1,1)/(s256-s56)',
c     & CLL(c2_56,1,1,1,1)/(s134-s56),C0(c2_56)*(s134-s56)

      if (docheck) then
        do h3=1,2
        do h5=1,2
        call ggZZcapture('12x34pp',h3,h5,1,2,3,4,5,6,
     &   CLL(c12_34,2,2,h3,h5)-two*mt**2*tri12_34rat(2,2,h3,h5),
     &   +two*tri12_34rat(2,2,h3,h5),czip)
        enddo
        enddo
      endif

 
      if (docheck) then
        call ggZZparsecheck('LL',-1,-1,CLL(c1_34,1,1,:,:),
     &   '(s134, s34,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,+1,CLL(c1_34,1,2,:,:),
     &   '(s134, s34,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,CLL(c1_34,2,1,:,:),
     &   '(s134, s34,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,CLL(c1_34,2,2,:,:),
     &   '(s134, s34,   0;mtsq,mtsq,mtsq)')

        call ggZZparsecheck('LL',-1,-1,CLL(c1_56,1,1,:,:),
     &   '(s156, s56,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,+1,CLL(c1_56,1,2,:,:),
     &   '(s156, s56,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,CLL(c1_56,2,1,:,:),
     &   '(s156, s56,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,CLL(c1_56,2,2,:,:),
     &   '(s156, s56,   0;mtsq,mtsq,mtsq)')

        call ggZZparsecheck('LL',-1,-1,CLL(c2_34,1,1,:,:),
     &   '( s34,s156,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,+1,CLL(c2_34,1,2,:,:),
     &   '( s34,s156,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,CLL(c2_34,2,1,:,:),
     &   '( s34,s156,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,CLL(c2_34,2,2,:,:),
     &   '( s34,s156,   0;mtsq,mtsq,mtsq)')

        call ggZZparsecheck('LL',-1,-1,CLL(c2_56,1,1,:,:),
     &   '( s56,s134,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,+1,CLL(c2_56,1,2,:,:),
     &   '( s56,s134,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,CLL(c2_56,2,1,:,:),
     &   '( s56,s134,   0;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,CLL(c2_56,2,2,:,:),
     &   '( s56,s134,   0;mtsq,mtsq,mtsq)')
     
        call ggZZparsecheck('LL',-1,-1,CLL(c12_34,1,1,:,:),
     &   '( s56, s34, s12;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',-1,+1,CLL(c12_34,1,2,:,:),
     &   '( s56, s34, s12;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,-1,CLL(c12_34,2,1,:,:),
     &   '( s56, s34, s12;mtsq,mtsq,mtsq)')
        call ggZZparsecheck('LL',+1,+1,CLL(c12_34,2,2,:,:),
     &   '( s56, s34, s12;mtsq,mtsq,mtsq)')
      endif
      
      tri(:,:,:,:,:)=czip

c--- only compute finite part here
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      do j=1,6
      tri(h1,h2,h3,h5,0)=tri(h1,h2,h3,h5,0)+CLL(j,h1,h2,h3,h5)*C0(j)
      enddo
      enddo
      enddo
      enddo
      enddo
      
c      write(6,*) 'tri(1,1,1,1,0)',tri(1,1,1,1,0)      

      return
      end

