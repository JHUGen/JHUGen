      subroutine WWZZ(n1,n2,n3,n4,n5,n6,n7,n8,za,zb,WWZZamp,srWW)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'zcouple.f'
      include 'runstring.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'first.f'
      include 'WWbits.f'
      double precision s34,s56,s17,s28,s3456,t3,t4,
     & xl1,xr1,xq1,xl2,xr2,xq2
      double complex WWZZamp(2,2),ggWW(2,2),
     & prop34,prop56,propw17,propw28,propWBF,prop3456,propz3456,
     & zab2,zba2,srWW(2,2),srL,srR,srggWW34(2,2),srggWW56(2,2),rxw,
     & sqzmass
      integer h34,h56,i1,i2,i3,i4,i5,i6,i7,i8,
     & n1,n2,n3,n4,n5,n6,n7,n8
      double complex, save:: ZZ3456(2,2)
      double complex anomzzamp,anomwwamp
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
      if (index(runstring,'mad') .gt. 0) then
        sqzmass=dcmplx(zmass**2,0d0)
      else
        sqzmass=czmass2
      endif

C---setting up couplings dependent on whether we are doing 34-line or 56-line
      if (n3+n4 == 7) then
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
      if (first) then
        first=.false.
        ZZ3456(1,1)=xl1*xl2
        ZZ3456(1,2)=xl1*xr2
        ZZ3456(2,1)=xr1*xl2
        ZZ3456(2,2)=xr1*xr2
      endif
      
      rxw=sqrt((cone-cxw)/cxw)
      s34=s(n3,n4)
      s56=s(n5,n6)
      s17=s(n1,n7)
      s28=s(n2,n8)
      s3456=t4(n3,n4,n5,n6)
      prop34=dcmplx(s34)-czmass2
      prop56=dcmplx(s56)-czmass2
      propw17=dcmplx(s17)-cwmass2
      propw28=dcmplx(s28)-cwmass2
      propz3456=dcmplx(s3456)-czmass2
      prop3456=dcmplx(s3456-hmass**2,hmass*hwidth)
      propWBF=propw17*propw28*prop34*prop56

!       print *,""
!       print *, "ww prop17",1d0/cdabs(propw17)**2
!       print *, "ww prop28",1d0/cdabs(propw28)**2
!       print *, "ww prop34",1d0/cdabs(prop34)**2
!       print *, "ww prop56",1d0/cdabs(prop56)**2
!       print *, "ww prop3456",1d0/cdabs(propz3456)**2
      
C----setup couplings and propagators
c      ggWW(1,1)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*l1**2)/prop34/prop56
c      ggWW(1,2)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*l1*r1)/prop34/prop56
c      ggWW(2,1)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*r1*l1)/prop34/prop56
c      ggWW(2,2)=dcmplx(q1**2/(s34*s56))+dcmplx(rxw*r1**2)/prop34/prop56

c--- Make sure WWZA vertices included
      ggWW(1,1)=(dcmplx(xq1/s34)+rxw*dcmplx(xl1)/prop34)
     &         *(dcmplx(xq2/s56)+rxw*dcmplx(xl2)/prop56)
      ggWW(1,2)=(dcmplx(xq1/s34)+rxw*dcmplx(xl1)/prop34)
     &         *(dcmplx(xq2/s56)+rxw*dcmplx(xr2)/prop56)
      ggWW(2,1)=(dcmplx(xq1/s34)+rxw*dcmplx(xr1)/prop34)
     &         *(dcmplx(xq2/s56)+rxw*dcmplx(xl2)/prop56)
      ggWW(2,2)=(dcmplx(xq1/s34)+rxw*dcmplx(xr1)/prop34)
     &         *(dcmplx(xq2/s56)+rxw*dcmplx(xr2)/prop56)

      srggWW34(1,1)=(dcmplx(xq2/s3456)+rxw*dcmplx(xl2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xl1)/prop34)
      srggWW34(1,2)=(dcmplx(xq2/s3456)+rxw*dcmplx(xl2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xl2*xr1)/prop34)
      srggWW34(2,1)=(dcmplx(xq2/s3456)+rxw*dcmplx(xr2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xl1)/prop34)
      srggWW34(2,2)=(dcmplx(xq2/s3456)+rxw*dcmplx(xr2)/propz3456)
     &         *(dcmplx(xq1*xq2/s34)+dcmplx(xr2*xr1)/prop34)

      srggWW56(1,1)=(dcmplx(xq1/s3456)+rxw*dcmplx(xl1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xl2)/prop56)
      srggWW56(1,2)=(dcmplx(xq1/s3456)+rxw*dcmplx(xl1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xl1*xr2)/prop56)
      srggWW56(2,1)=(dcmplx(xq1/s3456)+rxw*dcmplx(xr1)/propz3456)
     &         *(dcmplx(xq1*xq2/s56)+dcmplx(xr1*xl2)/prop56)
      srggWW56(2,2)=(dcmplx(xq1/s3456)+rxw*dcmplx(xr1)/propz3456)
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






! !      MARKUS: removed decay     
!          WWZZamp(h34,h56)=WWZZamp(h34,h56)
!      & -2d0*sqzmass/cxw**2                   !   <----  where is the factor 2.0 coming from::  ZZH amplitudes also have a global factor 2 in ZZ1728 (for no good reason)
!      & *za(i7,i8)*zb(i2,i1)
!      & /(propw17*propw28)
!      & /zmass*sqrt(cxw)*sqrt(1d0-cxw) !-- only  MCFM PR




C----Higgs contribution
        WWZZamp(h34,h56)=WWZZamp(h34,h56)
     &    -2d0*sqzmass/cxw**2*ZZ3456(h34,h56)
     &    *za(i7,i8)*zb(i2,i1)*za(i3,i5)*zb(i6,i4) !--F
!      &     *anomzzamp(i3,i4,i5,i6,s3456,s(i3,i4),s(i5,i6),za,zb) !--F
!      &     *anomwwamp(i7,i1,i8,i2,s3456,s(i7,i1),s(i8,i2),za,zb) !--F
     &    /(propWBF*prop3456)*Hbit


        enddo
      enddo

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
      srWW(1,h56)=2d0/(cxw*propw17*propw28)
     & *srggWW56(1,h56)*srL(i1,i2,i3,i4,i5,i6,i7,i8)*BBit
      srWW(2,h56)=2d0/(cxw*propw17*propw28)
     & *srggWW56(2,h56)*srR(i1,i2,i3,i4,i5,i6,i7,i8)*BBit
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
      srWW(h34,1)=srWW(h34,1)+2d0/(cxw*propw17*propw28)
     & *srggWW34(1,h34)*srL(i1,i2,i5,i6,i3,i4,i7,i8)*BBit
      srWW(h34,2)=srWW(h34,2)+2d0/(cxw*propw17*propw28)
     & *srggWW34(2,h34)*srR(i1,i2,i5,i6,i3,i4,i7,i8)*BBit
      enddo
      
      return
      end
