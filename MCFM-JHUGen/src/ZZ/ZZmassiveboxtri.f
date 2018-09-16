      subroutine ZZmassiveboxtri(j1,j2,j3,j4,j5,j6,za,zb,mt,totrat,
     & box,tri)
      implicit none
      include 'types.f'
c--- coefficients of boxes and triangles using the integral basis
c--- that includes 6-dimensional boxes
c--- (-,+) and (+,-) helicities of gluons only
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'docheck.f'
      include 'ggZZintegrals.f'
      include 'ZZdlabels.f'
      include 'ZZclabels.f'
      include 'first.f'
      integer:: j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,j,itot,irat
      real(dp):: mt,mtsq
      complex(dp):: Xpp(2,2),Xmp(2,2),Xpm(2,2),Xmm(2,2),
     & Zmp(2,2,2),Zpm(2,2,2),totrat(2,2,2,2),
     & box(2,2,2,2,-2:0),tri(2,2,2,2,-2:0),d(2,2,2,2,5),c(2,2,2,2,6),
     & Xrat(2,2,2,2),drat(2,2,2,2,3),crat(2,2,2,2,6)
      parameter(itot=1,irat=2)
      
      mtsq=mt**2
      
      if (first) then
c      write(6,*) 'mtsq',mtsq
c      write(6,*) 'musq',musq
c        first=.false. 
c        call qlinit
c      pause
      endif


c--- Box 1: note that coefficient of this box is unchanged
      call ZZbox1LL(j1,j2,j3,j4,j5,j6,za,zb,mt,Xpp,Xmp,Xpm,Xmm,Xrat)
      d(1,2,:,:,d1_34_2)=Xmp(:,:)
      d(2,1,:,:,d1_34_2)=Xpm(:,:)


c--- d=6 Box 2
      call ZZD062x1x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      d(1,2,:,:,d6_2_1_34)=Xmp(:,:)
      d(2,1,:,:,d6_2_1_34)=Xpm(:,:)
      
c--- d=6 Box 3
      call ZZD062x1x34LLmp(j2,j1,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      d(1,2,:,:,d6_1_2_34)=Xpm(:,:)
      d(2,1,:,:,d6_1_2_34)=Xmp(:,:)


c---- From here on, returned arrays include complete result and
c---- also the contribution of that integral to the rational part
c---- (note change from Xmp,Xpm to Zmp,Zpm)

c---- Box 2
      call ZZD02x1x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Zmp,Zpm)
      d(1,2,:,:,d2_1_34)=Zmp(:,:,itot)
      d(2,1,:,:,d2_1_34)=Zpm(:,:,itot)
      drat(1,2,:,:,d2_1_34)=Zmp(:,:,irat)
      drat(2,1,:,:,d2_1_34)=Zpm(:,:,irat)


c---- Box 3
      call ZZD02x1x34LLmp(j2,j1,j3,j4,j5,j6,za,zb,mt,Zmp,Zpm)
      d(1,2,:,:,d1_2_34)=Zpm(:,:,itot)
      d(2,1,:,:,d1_2_34)=Zmp(:,:,itot)
      drat(1,2,:,:,d1_2_34)=Zpm(:,:,irat)
      drat(2,1,:,:,d1_2_34)=Zmp(:,:,irat)
   
       
c---- Triangle 1
      call ZZC01x2LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Zmp,Zpm)
      c(1,2,:,:,c1_2)=Zmp(:,:,itot)
      c(2,1,:,:,c1_2)=Zpm(:,:,itot)
      crat(1,2,:,:,c1_2)=Zmp(:,:,irat)
      crat(2,1,:,:,c1_2)=Zpm(:,:,irat)
      
c---- Triangle 3
      call ZZC01x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Zmp,Zpm)
      c(1,2,:,:,c1_34)=Zmp(:,:,itot)
      c(2,1,:,:,c1_34)=Zpm(:,:,itot)
      crat(1,2,:,:,c1_34)=Zmp(:,:,irat)
      crat(2,1,:,:,c1_34)=Zpm(:,:,irat)

c---- Triangle 4
      call ZZC01x34LLmp(j2,j1,j3,j4,j5,j6,za,zb,mt,Zmp,Zpm)
      c(1,2,:,:,c2_34)=Zpm(:,:,itot)
      c(2,1,:,:,c2_34)=Zmp(:,:,itot)
      crat(1,2,:,:,c2_34)=Zpm(:,:,irat)
      crat(2,1,:,:,c2_34)=Zmp(:,:,irat)

c---- Triangle 5
      call ZZC01x34LLmp(j1,j2,j5,j6,j3,j4,za,zb,mt,Zmp,Zpm)
      do h3=1,2
      do h4=1,2
      c(1,2,h3,h4,c1_56)=Zmp(h4,h3,itot)
      c(2,1,h3,h4,c1_56)=Zpm(h4,h3,itot)
      crat(1,2,h3,h4,c1_56)=Zmp(h4,h3,irat)
      crat(2,1,h3,h4,c1_56)=Zpm(h4,h3,irat)
      enddo
      enddo

c---- Triangle 6
      call ZZC01x34LLmp(j2,j1,j5,j6,j3,j4,za,zb,mt,Zmp,Zpm)
      do h3=1,2
      do h4=1,2
      c(1,2,h3,h4,c2_56)=Zpm(h4,h3,itot)
      c(2,1,h3,h4,c2_56)=Zmp(h4,h3,itot)
      crat(1,2,h3,h4,c2_56)=Zpm(h4,h3,irat)
      crat(2,1,h3,h4,c2_56)=Zmp(h4,h3,irat)
      enddo
      enddo

c---- Triangle 2      

c------ Fill the coefficients for the massless case
      call ZZC012x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Zmp,Zpm)
      c(1,2,:,:,c12_34)=Zmp(:,:,itot)
      c(2,1,:,:,c12_34)=Zpm(:,:,itot)
      crat(1,2,:,:,c12_34)=Zmp(:,:,irat)
      crat(2,1,:,:,c12_34)=Zpm(:,:,irat)

c------ Compute mt^2 dependence from rational contributions and add
      do h1=1,2
      h2=3-h1
      do h3=1,2
      do h4=1,2
c        write(6,*) h1,h2,h3,h4,crat(h1,h2,h3,h4,c12_34)
        crat(h1,h2,h3,h4,c12_34)=totrat(h1,h2,h3,h4)
     &    -3._dp*drat(h1,h2,h3,h4,d2_1_34)-crat(h1,h2,h3,h4,c1_2)
     &    -crat(h1,h2,h3,h4,c1_34)-crat(h1,h2,h3,h4,c2_34)
     &    -crat(h1,h2,h3,h4,c1_56)-crat(h1,h2,h3,h4,c2_56)
c        write(6,*) h1,h2,h3,h4,crat(h1,h2,h3,h4,c12_34)
c        write(6,*)
        c(h1,h2,h3,h4,c12_34)=c(h1,h2,h3,h4,c12_34)
     &   +2._dp*mtsq*crat(h1,h2,h3,h4,c12_34)
      enddo
      enddo
      enddo

      if (docheck) then
        do h3=1,2
        do h4=1,2
        call ggZZcapture('12x34',h3,h4,1,2,3,4,5,6,
     &   c(1,2,h3,h4,c12_34)-2._dp*mtsq*crat(1,2,h3,h4,c12_34),
     &   +2._dp*crat(1,2,h3,h4,c12_34),czip)
        enddo
        enddo
      endif

c--- only compute finite part here
      do h1=1,2
      h2=3-h1
      box(h1,h2,:,:,:)=czip
      tri(h1,h2,:,:,:)=czip
      do h3=1,2
      do h4=1,2
c---  boxes
      do j=1,5
      box(h1,h2,h3,h4,0)=box(h1,h2,h3,h4,0)+d(h1,h2,h3,h4,j)*D0(j)
      enddo
c---  triangles
      do j=1,6
      tri(h1,h2,h3,h4,0)=tri(h1,h2,h3,h4,0)+c(h1,h2,h3,h4,j)*C0(j)
      enddo
      enddo
      enddo
      enddo

      return
      end
      
