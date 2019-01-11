      subroutine WWZZ(n1,n2,n3,n4,n5,n6,n7,n8,za,zb,WWZZamp,srWW)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'runstring.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
!      include 'first.f'
      include 'WWbits.f'
      include 'spinzerohiggs_anomcoupl.f'
      double precision s34,s56,s17,s28,s3456,t3,t4,
     & xl1,xr1,xq1,xl2,xr2,xq2
      double complex WWZZamp(2,2),ggWW(2,2),propX3456,
     & prop34,prop56,propw17,propw28,propWBF,prop3456,propz3456,
     & zab2,zba2,srWW(2,2),rxw,
     & sqzmass,Amp_S_DK,Amp_S_PR,facHiggs,srL,srR,
     & srgmWW34(2,2),srgmWW56(2,2),srZWW34(2,2),srZWW56(2,2)
!      & ,srggWW34(2,2),srggWW56(2,2)
      integer h34,h56,i1,i2,i3,i4,i5,i6,i7,i8,
     & n1,n2,n3,n4,n5,n6,n7,n8
      double complex ZZ3456(2,2)
      double complex ZA3456(2,2)
      double complex AZ3456(2,2)
      double complex AA3456(2,2)
      double complex anomhzzamp,anomhzaamp,anomhaaamp,anomhwwamp
      double complex srL_anom,srR_anom
!$omp threadprivate(ZZ3456)
      t4(i1,i2,i3,i4)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)
     & +s(i2,i3)+s(i2,i4)+s(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)

      srL(i1,i2,i3,i4,i5,i6,i7,i8)=
     & (2d0*zab2(i7,i2,i8,i1)*za(i3,i8)*zb(i6,i4)*zba2(i2,i4,i6,i5)
     & -2d0*zab2(i8,i1,i7,i2)*za(i3,i7)*zb(i6,i4)*zba2(i1,i4,i6,i5)
     & +((zab2(i3,i1,i7,i6)-zab2(i3,i2,i8,i6))*za(i6,i5)
     & + (zab2(i3,i1,i7,i4)-zab2(i3,i2,i8,i4))*za(i4,i5))
     & *za(i7,i8)*zb(i2,i1)*zb(i6,i4))/t3(i4,i5,i6)
     &+(2d0*zab2(i8,i1,i7,i2)*za(i3,i5)*zb(i1,i4)*zba2(i6,i3,i5,i7)
     & -2d0*zab2(i7,i2,i8,i1)*za(i3,i5)*zb(i2,i4)*zba2(i6,i3,i5,i8)
     & +(zb(i6,i5)*(zab2(i5,i2,i8,i4)-zab2(i5,i1,i7,i4))
     & + zb(i6,i3)*(zab2(i3,i2,i8,i4)-zab2(i3,i1,i7,i4)))
     & *za(i7,i8)*za(i3,i5)*zb(i2,i1))/t3(i3,i5,i6)

      srR(i1,i2,i3,i4,i5,i6,i7,i8)=
     & (2d0*zab2(i7,i2,i8,i1)*zab2(i8,i4,i5,i6)*za(i5,i4)*zb(i3,i2)
     & -2d0*zab2(i7,i4,i5,i6)*zab2(i8,i1,i7,i2)*za(i5,i4)*zb(i3,i1)
     & +(zb(i5,i6)*(zba2(i3,i1,i7,i5)-zba2(i3,i2,i8,i5))
     & + zb(i4,i6)*(zba2(i3,i1,i7,i4)-zba2(i3,i2,i8,i4)))
     & *za(i7,i8)*za(i5,i4)*zb(i2,i1))/t3(i4,i5,i6)
     &+(2d0*zab2(i8,i1,i7,i2)*zab2(i5,i3,i6,i1)*za(i7,i4)*zb(i3,i6)
     & -2d0*zab2(i7,i2,i8,i1)*zab2(i5,i3,i6,i2)*za(i8,i4)*zb(i3,i6)
     & +(za(i5,i6)*(zba2(i6,i2,i8,i4)-zba2(i6,i1,i7,i4))
     & +za(i5,i3)*(zba2(i3,i2,i8,i4)-zba2(i3,i1,i7,i4)))
     & *za(i7,i8)*zb(i2,i1)*zb(i3,i6))/t3(i3,i5,i6)
C---end statement functions

c--- special fix for Madgraph check
!      if (index(runstring,'mad') .gt. 0) then
!       sqzmass=dcmplx(zmass**2,0d0)
!      else
         sqzmass=czmass2
!      endif
!      facHiggs = 2d0*sqzmass/cxw**2
      ! Multiply by 1=(4d0*cwmass2/vevsq*cxw/esq)
      facHiggs = 2d0*sqzmass*cxw**(-2)*(4d0*cwmass2/vevsq*cxw/esq)

C---setting up couplings dependent on whether we are doing 34-line or 56-line
      if ((n3+n4 == 7) .or. (n3+n4 == 9)) then
      xl1=l1
      xr1=r1
      xq1=q1
      xl2=l2
      xr2=r2
      xq2=q2
      elseif (n3+n4 == 11) then
      xl1=l2
      xr1=r2
      xq1=q2
      xl2=l1
      xr2=r1
      xq2=q1
      else
      write(6,*) 'Unexpected case jtwo3456.f'
      stop
      endif

      srWW(:,:)=czip
      ZZ3456(1,1)=xl1*xl2
      ZZ3456(1,2)=xl1*xr2
      ZZ3456(2,1)=xr1*xl2
      ZZ3456(2,2)=xr1*xr2
      ZA3456(1,1)=xl1*xq2
      ZA3456(1,2)=xl1*xq2
      ZA3456(2,1)=xr1*xq2
      ZA3456(2,2)=xr1*xq2
      AZ3456(1,1)=xq1*xl2
      AZ3456(1,2)=xq1*xr2
      AZ3456(2,1)=xq1*xl2
      AZ3456(2,2)=xq1*xr2
      AA3456(:,:)=xq1*xq2

      rxw=sqrt((cone-cxw)/cxw)
      s34=s(n3,n4)
      s56=s(n5,n6)
      s17=s(n1,n7)
      s28=s(n2,n8)
      s3456=t4(n3,n4,n5,n6)
      prop34=dcmplx(s34-zmass**2,zmass*zwidth)
      prop56=dcmplx(s56-zmass**2,zmass*zwidth)
      propw17=dcmplx(s17-wmass**2,wmass*wwidth)
      propw28=dcmplx(s28-wmass**2,wmass*wwidth)
      propz3456=dcmplx(s3456-zmass**2,zmass*zwidth)
      prop3456=dcmplx(s3456-hmass**2,hmass*hwidth)
      propWBF=propw17*propw28*prop34*prop56

      propX3456=dcmplx(s3456-h2mass**2,h2mass*h2width)
C----setup couplings and propagators
c      ggWW(1,1)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*l1**2)/prop34/prop56
c      ggWW(1,2)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*l1*r1)/prop34/prop56
c      ggWW(2,1)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*r1*l1)/prop34/prop56
c      ggWW(2,2)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*r1**2)/prop34/prop56



c--- Make sure WWZA vertices included
!       ggWW(1,1)=(dcmplx(xq1/s34)+rxw*dcmplx(xl1)/prop34)
!      &         *(dcmplx(xq2/s56)+rxw*dcmplx(xl2)/prop56)
!       ggWW(1,2)=(dcmplx(xq1/s34)+rxw*dcmplx(xl1)/prop34)
!      &         *(dcmplx(xq2/s56)+rxw*dcmplx(xr2)/prop56)
!       ggWW(2,1)=(dcmplx(xq1/s34)+rxw*dcmplx(xr1)/prop34)
!      &         *(dcmplx(xq2/s56)+rxw*dcmplx(xl2)/prop56)
!       ggWW(2,2)=(dcmplx(xq1/s34)+rxw*dcmplx(xr1)/prop34)
!      &         *(dcmplx(xq2/s56)+rxw*dcmplx(xr2)/prop56)
! 
!       print *, "before",ggWW(:,:)
      
      ggWW(1,1)=dcmplx(xq1/s34)*dcmplx(xq2/s56)
     &         +rxw*dcmplx(xl1)/prop34*dcmplx(xq2/s56)* dZAWpWm
     &         +dcmplx(xq1/s34)*rxw*dcmplx(xl2)/prop56* dZAWpWm
     &         +rxw*dcmplx(xl1)/prop34*rxw*dcmplx(xl2)/prop56* dZZWpWm
     
     
      ggWW(1,2)=dcmplx(xq1/s34)*dcmplx(xq2/s56)
     &        +rxw*dcmplx(xl1)/prop34*dcmplx(xq2/s56)* dZAWpWm
     &        +dcmplx(xq1/s34)*rxw*dcmplx(xr2)/prop56* dZAWpWm
     &        +rxw*dcmplx(xl1)/prop34*rxw*dcmplx(xr2)/prop56* dZZWpWm
     
     
      ggWW(2,1)=dcmplx(xq1/s34)*dcmplx(xq2/s56)
     &         +rxw*dcmplx(xr1)/prop34*dcmplx(xq2/s56)* dZAWpWm
     &         +dcmplx(xq1/s34)*rxw*dcmplx(xl2)/prop56* dZAWpWm
     &         +rxw*dcmplx(xr1)/prop34*rxw*dcmplx(xl2)/prop56* dZZWpWm
     
     
      ggWW(2,2)=dcmplx(xq1/s34)*dcmplx(xq2/s56)
     &         +rxw*dcmplx(xr1)/prop34*dcmplx(xq2/s56)* dZAWpWm
     &         +dcmplx(xq1/s34)*rxw*dcmplx(xr2)/prop56* dZAWpWm
     &         +rxw*dcmplx(xr1)/prop34*rxw*dcmplx(xr2)/prop56* dZZWpWm
     
!       print *, "after ",ggWW(:,:)
!       pause

c--- This is WW->Z/A->f fb(->Z/A->f'fb')+WW->Z/A->fb f(->Z/A->f'fb')
! !     old
!       srggWW34(1,1)=(dcmplx(xq2/s3456)+rxw*dcmplx(xl2)/propz3456)
!      &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xl1)/prop34)
!       srggWW34(1,2)=(dcmplx(xq2/s3456)+rxw*dcmplx(xl2)/propz3456)
!      &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xr1)/prop34)
!       srggWW34(2,1)=(dcmplx(xq2/s3456)+rxw*dcmplx(xr2)/propz3456)
!      &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xl1)/prop34)
!       srggWW34(2,2)=(dcmplx(xq2/s3456)+rxw*dcmplx(xr2)/propz3456)
!      &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xr1)/prop34)
! 
!       srggWW56(1,1)=(dcmplx(xq1/s3456)+rxw*dcmplx(xl1)/propz3456)
!      &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xl2)/prop56)
!       srggWW56(1,2)=(dcmplx(xq1/s3456)+rxw*dcmplx(xl1)/propz3456)
!      &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xr2)/prop56)
!       srggWW56(2,1)=(dcmplx(xq1/s3456)+rxw*dcmplx(xr1)/propz3456)
!      &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xl2)/prop56)
!       srggWW56(2,2)=(dcmplx(xq1/s3456)+rxw*dcmplx(xr1)/propz3456)
!      &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xr2)/prop56)

     
!     new
      srgmWW34(1,1)=(dcmplx(xq2/s3456))
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xl1)/prop34)
      srgmWW34(1,2)=(dcmplx(xq2/s3456))
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xr1)/prop34)
      srgmWW34(2,1)=(dcmplx(xq2/s3456))
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xl1)/prop34)
      srgmWW34(2,2)=(dcmplx(xq2/s3456))
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xr1)/prop34)

      srgmWW56(1,1)=(dcmplx(xq1/s3456))
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xl2)/prop56)
      srgmWW56(1,2)=(dcmplx(xq1/s3456))
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xr2)/prop56)
      srgmWW56(2,1)=(dcmplx(xq1/s3456))
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xl2)/prop56)
      srgmWW56(2,2)=(dcmplx(xq1/s3456))
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xr2)/prop56)

     

      srZWW34(1,1)=(rxw*dcmplx(xl2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xl1)/prop34)
      srZWW34(1,2)=(rxw*dcmplx(xl2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xr1)/prop34)
      srZWW34(2,1)=(rxw*dcmplx(xr2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xl1)/prop34)
      srZWW34(2,2)=(rxw*dcmplx(xr2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xr1)/prop34)

      srZWW56(1,1)=(rxw*dcmplx(xl1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xl2)/prop56)
      srZWW56(1,2)=(rxw*dcmplx(xl1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xr2)/prop56)
      srZWW56(2,1)=(rxw*dcmplx(xr1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xl2)/prop56)
      srZWW56(2,2)=(rxw*dcmplx(xr1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xr2)/prop56)
     
     
      i1=n1
      i2=n2
      i7=n7
      i8=n8
      do h56=1,2
        if (h56.eq.1) then
          i5=n5
          i6=n6
        elseif (h56.eq.2) then
          i5=n6
          i6=n5
        endif
        do h34=1,2
         if (h34.eq.1) then
            i3=n3
            i4=n4
         elseif (h34.eq.2) then
            i3=n4
            i4=n3
         endif

         WWZZamp(h34,h56)=
     &    -2d0/cxw*(
     &     (2d0*za(i3,i5)*zb(i6,i4)*za(i7,i8)*zb(i2,i1)
     &        -za(i3,i7)*zb(i1,i4)*za(i5,i8)*zb(i2,i6)
     &        -za(i3,i8)*zb(i2,i4)*za(i5,i7)*zb(i1,i6)))
     &    *ggWW(h34,h56)/(propw17*propw28)*Bbit

C----Higgs contribution
C----First resonance
C-- MARKUS: this is the old (original) MCFM code
!          WWZZamp(h34,h56)=WWZZamp(h34,h56)
!      &    -2d0*sqzmass/cxw**2*ZZ3456(h34,h56)
!      &    *za(i7,i8)*zb(i2,i1)*za(i3,i5)*zb(i6,i4)
!      &    /(propWBF*prop3456)*Hbit
!

         Amp_S_PR=czip
         Amp_S_DK=czip
         if( hmass.ge.zip .and. channeltoggle_stu.ne.1 ) then
           if( AnomalCouplPR.eq.1 ) then
      Amp_S_PR=-anomhwwamp(i7,i1,i8,i2,1,s3456,s(i7,i1),s(i8,i2),za,zb)
     &         /(propw17*propw28)
           else
      Amp_S_PR=za(i7,i8)*zb(i2,i1)/(propw17*propw28)
           endif
           if( AnomalCouplDK.eq.1 ) then
      Amp_S_DK=
     & -anomhzzamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*prop56)*ZZ3456(h34,h56)
     & +anomhzaamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*s(i5,i6))*ZA3456(h34,h56)
     & +anomhzaamp(i5,i6,i3,i4,1,s3456,s(i5,i6),s(i3,i4),za,zb)
     & /(s(i3,i4)*prop56)*AZ3456(h34,h56)
     & -anomhaaamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(s(i3,i4)*s(i5,i6))*AA3456(h34,h56)
           else
      Amp_S_DK=za(i3,i5)*zb(i6,i4)/(prop34*prop56)*ZZ3456(h34,h56)
           endif
      WWZZamp(h34,h56)=WWZZamp(h34,h56)
     & -facHiggs*Amp_S_DK*Amp_S_PR/prop3456*Hbit
         endif

C----Second resonance
         Amp_S_PR=czip
         Amp_S_DK=czip
         if( h2mass.ge.zip .and. channeltoggle_stu.ne.1 ) then
           if( AnomalCouplPR.eq.1 ) then
      Amp_S_PR=-anomhwwamp(i7,i1,i8,i2,2,s3456,s(i7,i1),s(i8,i2),za,zb)
     &         /(propw17*propw28)
           else
      Amp_S_PR=za(i7,i8)*zb(i2,i1)/(propw17*propw28)
           endif
           if( AnomalCouplDK.eq.1 ) then
      Amp_S_DK=
     & -anomhzzamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*prop56)*ZZ3456(h34,h56)
     & +anomhzaamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*s(i5,i6))*ZA3456(h34,h56)
     & +anomhzaamp(i5,i6,i3,i4,2,s3456,s(i5,i6),s(i3,i4),za,zb)
     & /(s(i3,i4)*prop56)*AZ3456(h34,h56)
     & -anomhaaamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(s(i3,i4)*s(i5,i6))*AA3456(h34,h56)
           else
      Amp_S_DK=za(i3,i5)*zb(i6,i4)/(prop34*prop56)*ZZ3456(h34,h56)
           endif
      WWZZamp(h34,h56)=WWZZamp(h34,h56)
     & -facHiggs*Amp_S_DK*Amp_S_PR/propX3456*Hbit
         endif
        enddo
      enddo


C----Background contribution


!       do h56=1,2
!       if (h56==1) then 
!         i5=n5
!         i6=n6
!       elseif (h56==2) then
!         i5=n6
!         i6=n5
!       endif
!       i3=n3
!       i4=n4
!       srWW(1,h56)=2d0/(cxw*propw17*propw28)
!      & *srggWW56(1,h56)*srL(i1,i2,i3,i4,i5,i6,i7,i8)*BBit
!       srWW(2,h56)=2d0/(cxw*propw17*propw28)
!      & *srggWW56(2,h56)*srR(i1,i2,i3,i4,i5,i6,i7,i8)*BBit
!       enddo
! 
!       do h34=1,2
!       if (h34==1) then 
!         i3=n3
!         i4=n4
!       elseif (h34==2) then
!         i3=n4
!         i4=n3
!       endif
!       i5=n5
!       i6=n6
!       srWW(h34,1)=srWW(h34,1)+2d0/(cxw*propw17*propw28)
!      & *srggWW34(1,h34)*srL(i1,i2,i5,i6,i3,i4,i7,i8)*BBit
!       srWW(h34,2)=srWW(h34,2)+2d0/(cxw*propw17*propw28)
!      & *srggWW34(2,h34)*srR(i1,i2,i5,i6,i3,i4,i7,i8)*BBit
!       enddo
!       
!       print *, "old srWW",srWW(:,:)



      do h56=1,2
      if (h56.eq.1) then
        i5=n5
        i6=n6
      elseif (h56.eq.2) then
        i5=n6
        i6=n5
      endif
      i3=n3
      i4=n4
      srWW(1,h56)=2d0/(cxw*propw17*propw28)*(
     & +srgmWW56(1,h56)
     & *srL_anom(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,dP_A,dM_A,dV_A,dFour_A)
     & +srZWW56(1,h56)
     & *srL_anom(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,dP_Z,dM_Z,dV_Z,dFour_Z)
     & )*BBit
     
      srWW(2,h56)=2d0/(cxw*propw17*propw28)*(
     & +srgmWW56(2,h56)
     & *srR_anom(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,dP_A,dM_A,dV_A,dFour_A)
     & +srZWW56(2,h56)
     & *srR_anom(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,dP_Z,dM_Z,dV_Z,dFour_Z)
     & )*BBit
        

      enddo

      
      
      do h34=1,2
      if (h34.eq.1) then
        i3=n3
        i4=n4
      elseif (h34.eq.2) then
        i3=n4
        i4=n3
      endif
      i5=n5
      i6=n6
      srWW(h34,1)=srWW(h34,1)+2d0/(cxw*propw17*propw28)*(
     & +srgmWW34(1,h34)
     & *srL_anom(i1,i2,i5,i6,i3,i4,i7,i8,za,zb,dP_A,dM_A,dV_A,dFour_A)
     & +srZWW34(1,h34)
     & *srL_anom(i1,i2,i5,i6,i3,i4,i7,i8,za,zb,dP_Z,dM_Z,dV_Z,dFour_Z)
     & )*BBit
     
      srWW(h34,2)=srWW(h34,2)+2d0/(cxw*propw17*propw28)*(
     & +srgmWW34(2,h34)
     & *srR_anom(i1,i2,i5,i6,i3,i4,i7,i8,za,zb,dP_A,dM_A,dV_A,dFour_A)
     & +srZWW34(2,h34)
     & *srR_anom(i1,i2,i5,i6,i3,i4,i7,i8,za,zb,dP_Z,dM_Z,dV_Z,dFour_Z)
     & )*BBit
     
      enddo

!       print *, "new srWW",srWW(:,:)
!       pause      
      
      return
      end

      
      
      
      
      !srL/R_anom replaces srL/R and includes four anomalous couplings (SM: dP=dM,dV=1, dFour=0)
      function srL_anom(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,dP,dM,dV,dFour)
      implicit none 
      include 'constants.f'     
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex srL_anom,dFour,dP,dM,dV
      double complex, parameter :: cI=(0d0,1d0)
      integer i1,i2,i3,i4,i5,i6,i7,i8
      double precision t356,t456

        t356=s(i3,i5)+s(i3,i6)+s(i5,i6)
        t456=s(i4,i5)+s(i4,i6)+s(i5,i6)

        srL_anom = (cI*dFour*(-(za(i3,i8)*
     -          (za(i4,i5)*zb(i1,i4) - za(i5,i6)*zb(i1,i6))*
     -          (za(i4,i7)*zb(i2,i4) + za(i5,i7)*zb(i2,i5) + 
     -            za(i6,i7)*zb(i2,i6))) + 
     -       za(i3,i7)*((za(i4,i8)*zb(i1,i4) + 
     -             za(i5,i8)*zb(i1,i5) + za(i6,i8)*zb(i1,i6))*
     -           (za(i4,i5)*zb(i2,i4) - za(i5,i6)*zb(i2,i6))
     -           + za(i3,i8)*
     -           (za(i4,i5)*
     -              (-(zb(i1,i4)*zb(i2,i3)) + 
     -                zb(i1,i3)*zb(i2,i4)) + 
     -             za(i5,i6)*
     -              (zb(i1,i6)*zb(i2,i3) - 
     -                zb(i1,i3)*zb(i2,i6)))))*zb(i4,i6))/t456
     -   + (cI*dFour*za(i3,i5)*
     -     ((za(i3,i7)*(za(i3,i8)*zb(i1,i4)*zb(i2,i3) - 
     -             (za(i3,i8)*zb(i1,i3) + 
     -                za(i4,i8)*zb(i1,i4) + 
     -                za(i5,i8)*zb(i1,i5) + 
     -                za(i6,i8)*zb(i1,i6))*zb(i2,i4)) + 
     -          za(i3,i8)*zb(i1,i4)*
     -           (za(i4,i7)*zb(i2,i4) + za(i5,i7)*zb(i2,i5) + 
     -             za(i6,i7)*zb(i2,i6)))*zb(i3,i6) + 
     -       (za(i3,i7)*za(i5,i8)*zb(i1,i4)*zb(i2,i3) - 
     -          za(i3,i8)*za(i5,i7)*zb(i1,i3)*zb(i2,i4) - 
     -          za(i4,i8)*za(i5,i7)*zb(i1,i4)*zb(i2,i4) + 
     -          za(i4,i7)*za(i5,i8)*zb(i1,i4)*zb(i2,i4) - 
     -          za(i5,i7)*za(i5,i8)*zb(i1,i5)*zb(i2,i4) - 
     -          za(i5,i7)*za(i6,i8)*zb(i1,i6)*zb(i2,i4) + 
     -          za(i5,i7)*za(i5,i8)*zb(i1,i4)*zb(i2,i5) + 
     -          za(i5,i8)*za(i6,i7)*zb(i1,i4)*zb(i2,i6))*
     -        zb(i5,i6)))/t356 + 
     -  (za(i3,i5)*(zb(i3,i6)*
     -        (dP*za(i1,i3)*za(i7,i8)*zb(i1,i2)*zb(i1,i4) + 
     -          (dM*za(i2,i7)*za(i3,i8)*zb(i1,i2) - 
     -             dV*za(i3,i8)*
     -              (za(i4,i7)*zb(i1,i4) + 
     -                za(i5,i7)*zb(i1,i5) + 
     -                za(i6,i7)*zb(i1,i6)) - 
     -             dM*za(i7,i8)*
     -              (za(i2,i3)*zb(i1,i2) + 
     -                za(i3,i8)*zb(i1,i8)))*zb(i2,i4) + 
     -          za(i3,i7)*
     -           (dV*za(i3,i8)*
     -              (zb(i1,i4)*zb(i2,i3) - 
     -                zb(i1,i3)*zb(i2,i4)) + 
     -             zb(i1,i4)*
     -              (dV*za(i4,i8)*zb(i2,i4) + 
     -                dV*za(i5,i8)*zb(i2,i5) + 
     -                dV*za(i6,i8)*zb(i2,i6) - 
     -                dP*za(i7,i8)*zb(i2,i7)) + 
     -             dP*za(i7,i8)*zb(i1,i2)*zb(i4,i7)) - 
     -          dM*za(i3,i8)*za(i7,i8)*zb(i1,i2)*zb(i4,i8)) + 
     -       (dP*za(i1,i5)*za(i7,i8)*zb(i1,i2)*zb(i1,i4) + 
     -          dV*za(i3,i8)*za(i5,i7)*zb(i1,i4)*zb(i2,i3) + 
     -          dM*(za(i2,i7)*za(i5,i8) - 
     -             za(i2,i5)*za(i7,i8))*zb(i1,i2)*zb(i2,i4) - 
     -          dV*za(i3,i7)*za(i5,i8)*zb(i1,i3)*zb(i2,i4) + 
     -          dV*za(i4,i8)*za(i5,i7)*zb(i1,i4)*zb(i2,i4) - 
     -          dV*za(i4,i7)*za(i5,i8)*zb(i1,i4)*zb(i2,i4) - 
     -          dV*za(i5,i7)*za(i5,i8)*zb(i1,i5)*zb(i2,i4) - 
     -          dV*za(i5,i8)*za(i6,i7)*zb(i1,i6)*zb(i2,i4) - 
     -          dM*za(i5,i8)*za(i7,i8)*zb(i1,i8)*zb(i2,i4) + 
     -          dV*za(i5,i7)*za(i5,i8)*zb(i1,i4)*zb(i2,i5) + 
     -          dV*za(i5,i7)*za(i6,i8)*zb(i1,i4)*zb(i2,i6) - 
     -          dP*za(i5,i7)*za(i7,i8)*zb(i1,i4)*zb(i2,i7) + 
     -          dP*za(i5,i7)*za(i7,i8)*zb(i1,i2)*zb(i4,i7) - 
     -          dM*za(i5,i8)*za(i7,i8)*zb(i1,i2)*zb(i4,i8))*
     -        zb(i5,i6) + 
     -       dP*za(i1,i8)*zb(i1,i2)*zb(i1,i4)*
     -        (za(i3,i7)*zb(i3,i6) + za(i5,i7)*zb(i5,i6))))/
     -   t356 + (zb(i4,i6)*
     -     (dP*za(i1,i8)*za(i3,i7)*zb(i1,i2)*
     -        (-(za(i4,i5)*zb(i1,i4)) + za(i5,i6)*zb(i1,i6))
     -        + dP*za(i1,i3)*za(i7,i8)*zb(i1,i2)*
     -        (-(za(i4,i5)*zb(i1,i4)) + za(i5,i6)*zb(i1,i6))
     -        - dV*za(i3,i7)*za(i3,i8)*za(i4,i5)*zb(i1,i4)*
     -        zb(i2,i3) + 
     -       dV*za(i3,i7)*za(i3,i8)*za(i5,i6)*zb(i1,i6)*
     -        zb(i2,i3) - 
     -       dM*za(i2,i7)*za(i3,i8)*za(i4,i5)*zb(i1,i2)*
     -        zb(i2,i4) + 
     -       dM*za(i2,i3)*za(i4,i5)*za(i7,i8)*zb(i1,i2)*
     -        zb(i2,i4) + 
     -       dV*za(i3,i7)*za(i3,i8)*za(i4,i5)*zb(i1,i3)*
     -        zb(i2,i4) + 
     -       dV*za(i3,i8)*za(i4,i5)*za(i4,i7)*zb(i1,i4)*
     -        zb(i2,i4) - 
     -       dV*za(i3,i7)*za(i4,i5)*za(i4,i8)*zb(i1,i4)*
     -        zb(i2,i4) + 
     -       dV*za(i3,i8)*za(i4,i5)*za(i5,i7)*zb(i1,i5)*
     -        zb(i2,i4) + 
     -       dV*za(i3,i7)*za(i4,i8)*za(i5,i6)*zb(i1,i6)*
     -        zb(i2,i4) + 
     -       dV*za(i3,i8)*za(i4,i5)*za(i6,i7)*zb(i1,i6)*
     -        zb(i2,i4) + 
     -       dM*za(i3,i8)*za(i4,i5)*za(i7,i8)*zb(i1,i8)*
     -        zb(i2,i4) - 
     -       dV*za(i3,i7)*za(i4,i5)*za(i5,i8)*zb(i1,i4)*
     -        zb(i2,i5) + 
     -       dV*za(i3,i7)*za(i5,i6)*za(i5,i8)*zb(i1,i6)*
     -        zb(i2,i5) + 
     -       dM*za(i2,i7)*za(i3,i8)*za(i5,i6)*zb(i1,i2)*
     -        zb(i2,i6) - 
     -       dM*za(i2,i3)*za(i5,i6)*za(i7,i8)*zb(i1,i2)*
     -        zb(i2,i6) - 
     -       dV*za(i3,i7)*za(i3,i8)*za(i5,i6)*zb(i1,i3)*
     -        zb(i2,i6) - 
     -       dV*za(i3,i8)*za(i4,i7)*za(i5,i6)*zb(i1,i4)*
     -        zb(i2,i6) - 
     -       dV*za(i3,i7)*za(i4,i5)*za(i6,i8)*zb(i1,i4)*
     -        zb(i2,i6) - 
     -       dV*za(i3,i8)*za(i5,i6)*za(i5,i7)*zb(i1,i5)*
     -        zb(i2,i6) - 
     -       dV*za(i3,i8)*za(i5,i6)*za(i6,i7)*zb(i1,i6)*
     -        zb(i2,i6) + 
     -       dV*za(i3,i7)*za(i5,i6)*za(i6,i8)*zb(i1,i6)*
     -        zb(i2,i6) - 
     -       dM*za(i3,i8)*za(i5,i6)*za(i7,i8)*zb(i1,i8)*
     -        zb(i2,i6) + 
     -       dP*za(i3,i7)*za(i4,i5)*za(i7,i8)*zb(i1,i4)*
     -        zb(i2,i7) - 
     -       dP*za(i3,i7)*za(i5,i6)*za(i7,i8)*zb(i1,i6)*
     -        zb(i2,i7) - 
     -       dP*za(i3,i7)*za(i4,i5)*za(i7,i8)*zb(i1,i2)*
     -        zb(i4,i7) + 
     -       dM*za(i3,i8)*za(i4,i5)*za(i7,i8)*zb(i1,i2)*
     -        zb(i4,i8) + 
     -       dP*za(i3,i7)*za(i5,i6)*za(i7,i8)*zb(i1,i2)*
     -        zb(i6,i7) - 
     -       dM*za(i3,i8)*za(i5,i6)*za(i7,i8)*zb(i1,i2)*
     -        zb(i6,i8)))/t456  
            
      return
      end

      function srR_anom(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,dP,dM,dV,dFour)
      implicit none 
      include 'constants.f'     
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex srR_anom,dFour,dP,dM,dV
      double complex, parameter :: cI=(0d0,1d0)      
      integer i1,i2,i3,i4,i5,i6,i7,i8
      double precision t356,t456
      

        t356=s(i3,i5)+s(i3,i6)+s(i5,i6)
        t456=s(i4,i5)+s(i4,i6)+s(i5,i6)     
      
      srR_anom = (cI*dFour*za(i4,i5)*
     -     (-((za(i7,i8)*zb(i1,i3)*
     -             (za(i3,i4)*zb(i2,i3) - 
     -               za(i4,i5)*zb(i2,i5) - za(i4,i6)*zb(i2,i6)
     -               ) + 
     -            za(i4,i7)*zb(i1,i2)*
     -             (za(i4,i8)*zb(i3,i4) + 
     -               za(i5,i8)*zb(i3,i5) + za(i6,i8)*zb(i3,i6)
     -               ))*zb(i4,i6)) - 
     -       (za(i7,i8)*zb(i1,i3)*
     -           (za(i3,i5)*zb(i2,i3) + za(i4,i5)*zb(i2,i4) - 
     -             za(i5,i6)*zb(i2,i6)) + 
     -          za(i5,i7)*zb(i1,i2)*
     -           (za(i4,i8)*zb(i3,i4) + za(i5,i8)*zb(i3,i5) + 
     -             za(i6,i8)*zb(i3,i6)))*zb(i5,i6)))/t456 + 
     -  (za(i4,i5)*((-(dP*za(i1,i4)*za(i7,i8)*zb(i1,i2)*
     -             zb(i1,i3)) - 
     -          dM*za(i2,i7)*za(i4,i8)*zb(i1,i2)*zb(i2,i3) + 
     -          dM*za(i2,i4)*za(i7,i8)*zb(i1,i2)*zb(i2,i3) - 
     -          dV*za(i3,i8)*za(i4,i7)*zb(i1,i3)*zb(i2,i3) + 
     -          dV*za(i3,i7)*za(i4,i8)*zb(i1,i3)*zb(i2,i3) + 
     -          dV*za(i4,i7)*za(i4,i8)*zb(i1,i4)*zb(i2,i3) + 
     -          dV*za(i4,i8)*za(i5,i7)*zb(i1,i5)*zb(i2,i3) + 
     -          dV*za(i4,i8)*za(i6,i7)*zb(i1,i6)*zb(i2,i3) + 
     -          dM*za(i4,i8)*za(i7,i8)*zb(i1,i8)*zb(i2,i3) - 
     -          dV*za(i4,i7)*za(i4,i8)*zb(i1,i3)*zb(i2,i4) - 
     -          dV*za(i4,i7)*za(i5,i8)*zb(i1,i3)*zb(i2,i5) - 
     -          dV*za(i4,i7)*za(i6,i8)*zb(i1,i3)*zb(i2,i6) + 
     -          dP*za(i4,i7)*za(i7,i8)*zb(i1,i3)*zb(i2,i7) - 
     -          dP*za(i4,i7)*za(i7,i8)*zb(i1,i2)*zb(i3,i7) + 
     -          dM*za(i4,i8)*za(i7,i8)*zb(i1,i2)*zb(i3,i8))*
     -        zb(i4,i6) + 
     -       (-(dP*za(i1,i5)*za(i7,i8)*zb(i1,i2)*zb(i1,i3)) - 
     -          dM*za(i2,i7)*za(i5,i8)*zb(i1,i2)*zb(i2,i3) + 
     -          dM*za(i2,i5)*za(i7,i8)*zb(i1,i2)*zb(i2,i3) - 
     -          dV*za(i3,i8)*za(i5,i7)*zb(i1,i3)*zb(i2,i3) + 
     -          dV*za(i3,i7)*za(i5,i8)*zb(i1,i3)*zb(i2,i3) + 
     -          dV*za(i4,i7)*za(i5,i8)*zb(i1,i4)*zb(i2,i3) + 
     -          dV*za(i5,i7)*za(i5,i8)*zb(i1,i5)*zb(i2,i3) + 
     -          dV*za(i5,i8)*za(i6,i7)*zb(i1,i6)*zb(i2,i3) + 
     -          dM*za(i5,i8)*za(i7,i8)*zb(i1,i8)*zb(i2,i3) - 
     -          dV*za(i4,i8)*za(i5,i7)*zb(i1,i3)*zb(i2,i4) - 
     -          dV*za(i5,i7)*za(i5,i8)*zb(i1,i3)*zb(i2,i5) - 
     -          dV*za(i5,i7)*za(i6,i8)*zb(i1,i3)*zb(i2,i6) + 
     -          dP*za(i5,i7)*za(i7,i8)*zb(i1,i3)*zb(i2,i7) - 
     -          dP*za(i5,i7)*za(i7,i8)*zb(i1,i2)*zb(i3,i7) + 
     -          dM*za(i5,i8)*za(i7,i8)*zb(i1,i2)*zb(i3,i8))*
     -        zb(i5,i6) - 
     -       dP*za(i1,i8)*zb(i1,i2)*zb(i1,i3)*
     -        (za(i4,i7)*zb(i4,i6) + za(i5,i7)*zb(i5,i6))))/
     -   t456 + (cI*dFour*zb(i3,i6)*
     -     (za(i3,i4)*za(i7,i8)*
     -        (za(i3,i5)*zb(i1,i3) - za(i5,i6)*zb(i1,i6))*
     -        zb(i2,i3) + 
     -       za(i3,i5)*(-(za(i7,i8)*zb(i1,i3)*
     -             (za(i4,i5)*zb(i2,i5) + za(i4,i6)*zb(i2,i6))
     -             ) + za(i4,i7)*zb(i1,i2)*
     -           (za(i4,i8)*zb(i3,i4) + za(i5,i8)*zb(i3,i5) + 
     -             za(i6,i8)*zb(i3,i6))) + 
     -       za(i5,i6)*(za(i7,i8)*zb(i1,i6)*
     -           (za(i4,i5)*zb(i2,i5) + za(i4,i6)*zb(i2,i6))
     -           + za(i4,i7)*zb(i1,i2)*
     -           (za(i3,i8)*zb(i3,i6) + za(i4,i8)*zb(i4,i6) + 
     -             za(i5,i8)*zb(i5,i6)))))/t356 + 
     -  (zb(i3,i6)*(dP*za(i1,i8)*za(i4,i7)*zb(i1,i2)*
     -        (za(i3,i5)*zb(i1,i3) - za(i5,i6)*zb(i1,i6)) + 
     -       dP*za(i1,i4)*za(i7,i8)*zb(i1,i2)*
     -        (za(i3,i5)*zb(i1,i3) - za(i5,i6)*zb(i1,i6)) + 
     -       dM*za(i2,i7)*za(i3,i5)*za(i4,i8)*zb(i1,i2)*
     -        zb(i2,i3) - 
     -       dM*za(i2,i4)*za(i3,i5)*za(i7,i8)*zb(i1,i2)*
     -        zb(i2,i3) + 
     -       dV*za(i3,i5)*za(i3,i8)*za(i4,i7)*zb(i1,i3)*
     -        zb(i2,i3) - 
     -       dV*za(i3,i5)*za(i3,i7)*za(i4,i8)*zb(i1,i3)*
     -        zb(i2,i3) - 
     -       dV*za(i3,i5)*za(i4,i7)*za(i4,i8)*zb(i1,i4)*
     -        zb(i2,i3) - 
     -       dV*za(i3,i5)*za(i4,i8)*za(i5,i7)*zb(i1,i5)*
     -        zb(i2,i3) - 
     -       dV*za(i3,i8)*za(i4,i7)*za(i5,i6)*zb(i1,i6)*
     -        zb(i2,i3) - 
     -       dV*za(i3,i5)*za(i4,i8)*za(i6,i7)*zb(i1,i6)*
     -        zb(i2,i3) - 
     -       dM*za(i3,i5)*za(i4,i8)*za(i7,i8)*zb(i1,i8)*
     -        zb(i2,i3) + 
     -       dV*za(i3,i5)*za(i4,i7)*za(i4,i8)*zb(i1,i3)*
     -        zb(i2,i4) - 
     -       dV*za(i4,i7)*za(i4,i8)*za(i5,i6)*zb(i1,i6)*
     -        zb(i2,i4) + 
     -       dV*za(i3,i5)*za(i4,i7)*za(i5,i8)*zb(i1,i3)*
     -        zb(i2,i5) - 
     -       dV*za(i4,i7)*za(i5,i6)*za(i5,i8)*zb(i1,i6)*
     -        zb(i2,i5) - 
     -       dM*za(i2,i7)*za(i4,i8)*za(i5,i6)*zb(i1,i2)*
     -        zb(i2,i6) + 
     -       dM*za(i2,i4)*za(i5,i6)*za(i7,i8)*zb(i1,i2)*
     -        zb(i2,i6) + 
     -       dV*za(i3,i7)*za(i4,i8)*za(i5,i6)*zb(i1,i3)*
     -        zb(i2,i6) + 
     -       dV*za(i3,i5)*za(i4,i7)*za(i6,i8)*zb(i1,i3)*
     -        zb(i2,i6) + 
     -       dV*za(i4,i7)*za(i4,i8)*za(i5,i6)*zb(i1,i4)*
     -        zb(i2,i6) + 
     -       dV*za(i4,i8)*za(i5,i6)*za(i5,i7)*zb(i1,i5)*
     -        zb(i2,i6) + 
     -       dV*za(i4,i8)*za(i5,i6)*za(i6,i7)*zb(i1,i6)*
     -        zb(i2,i6) - 
     -       dV*za(i4,i7)*za(i5,i6)*za(i6,i8)*zb(i1,i6)*
     -        zb(i2,i6) + 
     -       dM*za(i4,i8)*za(i5,i6)*za(i7,i8)*zb(i1,i8)*
     -        zb(i2,i6) - 
     -       dP*za(i3,i5)*za(i4,i7)*za(i7,i8)*zb(i1,i3)*
     -        zb(i2,i7) + 
     -       dP*za(i4,i7)*za(i5,i6)*za(i7,i8)*zb(i1,i6)*
     -        zb(i2,i7) + 
     -       dP*za(i3,i5)*za(i4,i7)*za(i7,i8)*zb(i1,i2)*
     -        zb(i3,i7) - 
     -       dM*za(i3,i5)*za(i4,i8)*za(i7,i8)*zb(i1,i2)*
     -        zb(i3,i8) - 
     -       dP*za(i4,i7)*za(i5,i6)*za(i7,i8)*zb(i1,i2)*
     -        zb(i6,i7) + 
     -       dM*za(i4,i8)*za(i5,i6)*za(i7,i8)*zb(i1,i2)*
     -        zb(i6,i8)))/t356      
      
      
      return
      end

      
            
      
      
