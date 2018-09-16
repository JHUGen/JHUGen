      subroutine ZZmassivebox(j1,j2,j3,j4,j5,j6,za,zb,mt,box,drat)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'ggZZintegrals.f'
      include 'ZZdlabels.f'
      integer:: j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,j
      real(dp):: mt
      complex(dp):: Xpp(2,2),Xmp(2,2),Xpm(2,2),Xmm(2,2),
     & box(2,2,2,2,-2:0),d(2,2,2,2,3),
     & Xrat(2,2,2,2),drat(2,2,2,2,3)
      

      call ZZbox2LL(j1,j2,j3,j4,j5,j6,za,zb,mt,Xpp,Xmp,Xpm,Xmm,Xrat)
      d(2,2,:,:,d2_1_34)=Xpp(:,:)
      d(1,2,:,:,d2_1_34)=Xmp(:,:)
      d(2,1,:,:,d2_1_34)=Xpm(:,:)
      d(1,1,:,:,d2_1_34)=Xmm(:,:)
      drat(:,:,:,:,d2_1_34)=Xrat(:,:,:,:)
c      write(6,*) 'Xmp(1,1)',Xmp(1,1)
c      write(6,*) 'Xmp(1,2)',Xmp(1,2)
c      write(6,*) 'Xmp(2,1)',Xmp(2,1)
c      write(6,*) 'Xmp(2,2)',Xmp(2,2)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)
      
      call ZZbox2LL(j1,j2,j5,j6,j3,j4,za,zb,mt,Xpp,Xmp,Xpm,Xmm,Xrat)
      d(2,2,:,:,d1_2_34)=Xpp(:,:)
      d(1,2,:,:,d1_2_34)=Xmp(:,:)
      d(2,1,:,:,d1_2_34)=Xpm(:,:)
      d(1,1,:,:,d1_2_34)=Xmm(:,:)
      drat(:,:,:,:,d1_2_34)=Xrat(:,:,:,:)
c      write(6,*) 'Xmp(1,1)',Xmp(1,1)
c      write(6,*) 'Xmp(1,2)',Xmp(1,2)
c      write(6,*) 'Xmp(2,1)',Xmp(2,1)
c      write(6,*) 'Xmp(2,2)',Xmp(2,2)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)

      call ZZbox1LL(j1,j2,j3,j4,j5,j6,za,zb,mt,Xpp,Xmp,Xpm,Xmm,Xrat)
      d(2,2,:,:,d1_34_2)=Xpp(:,:)
      d(1,2,:,:,d1_34_2)=Xmp(:,:)
      d(2,1,:,:,d1_34_2)=Xpm(:,:)
      d(1,1,:,:,d1_34_2)=Xmm(:,:)
      drat(:,:,:,:,d1_34_2)=Xrat(:,:,:,:)
c      write(6,*) 'Xpp(1,1)',Xpp(1,1)
c      write(6,*) 'Xpp(1,2)',Xpp(1,2)
c      write(6,*) 'Xpp(2,1)',Xpp(2,1)
c      write(6,*) 'Xpp(2,2)',Xpp(2,2)
c      write(6,*) 'Xrat(1,1,1,1)',Xrat(1,1,1,1)

c--- These integrals are now initialized in ZZintegraleval      
c      do e=-2,0
c      Dint(1,e)
c     & =qlI4(s34,zip,zip,s56,s134,s12,mtsq,mtsq,mtsq,mtsq,musq,e)
c      enddo
c      do e=-2,0
c      Dint(2,e)
c     & =qlI4(s56,zip,zip,s34,s156,s12,mtsq,mtsq,mtsq,mtsq,musq,e)
c      enddo
c      do e=-2,0
c      Dint(3,e)
c     & =qlI4(s56,zip,s34,zip,s156,s134,mtsq,mtsq,mtsq,mtsq,musq,e)
c      enddo

c      write(6,*) 'check',D0(d2_1_34)/Dint(1,0)
c      write(6,*) 'check',D0(d1_2_34)/Dint(2,0)
c      write(6,*) 'check',D0(d1_34_2)/Dint(3,0)
c      pause

c--- DEBUG
c      write(6,*) 'd(1,1,1,1,d1_2_34)/s12/s234,D0(d1_2_34)*s12*s234',
c     & d(1,1,1,1,d1_2_34)/s12/s156,D0(d1_2_34)*s12*s156
c      write(6,*) 'd(1,1,1,1,d2_1_34)/s12/s134,D0(d2_1_34)*s12*s134',
c     & d(1,1,1,1,d2_1_34)/s12/s134,D0(d2_1_34)*s12*s134
c--- DEBUG

      box(:,:,:,:,:)=czip

c--- only compute finite part here
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,3
      box(h1,h2,h3,h4,0)=box(h1,h2,h3,h4,0)+d(h1,h2,h3,h4,j)*D0(j)
      enddo
      enddo
      enddo
      enddo
      enddo
      
c      write(6,*) 'box(1,1,1,1,0)',box(1,1,1,1,0)

      return
      end
      
