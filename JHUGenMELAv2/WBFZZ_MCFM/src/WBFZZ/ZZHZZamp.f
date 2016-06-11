      subroutine ZZHZZamp(n1,n2,n3,n4,n5,n6,n7,n8,za,zb,ZZHamp)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
c      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'first.f'
      double precision t4,s3456,s1734,s1756,htheta,
     & ZZ3456(2,2),ZZ1728(2,2,2,2),ZZ1734(2,2,2),ZZ2856(2,2,2)
      double complex ZZHamp(2,2,2,2,2,2),propWBF,
     & prop34,prop56,prop17,prop28,prop3456,prop1734,prop1756,fac
      save ZZ3456,ZZ1734,ZZ2856,ZZ1728,fac
C---order of indices jdu1,jdu2,h17,h28,h34,h56)
      integer h17,h28,h34,h56,i1,i2,i3,i4,i5,i6,i7,i8,
     & n1,n2,n3,n4,n5,n6,n7,n8,jdu1,jdu2
      double complex anomzzamp
C---begin statement functions
      t4(i1,i2,i3,i4)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)
     & +s(i2,i3)+s(i2,i4)+s(i3,i4)
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
c      htheta(s3456)=half+sign(half,s3456)
      htheta(s3456)=one ! propdebug
C---end statement functions
!$omp threadprivate(ZZ3456,ZZ1734,ZZ2856,ZZ1728,fac)
      
      if (first) then
      first=.false.
      ZZ3456(1,1)=2d0*l1*l2
      ZZ3456(1,2)=2d0*l1*r2
      ZZ3456(2,1)=2d0*r1*l2
      ZZ3456(2,2)=2d0*r1*r2
      do jdu1=1,2
      ZZ1734(jdu1,1,1)=2d0*L(jdu1)*l1
      ZZ1734(jdu1,1,2)=2d0*L(jdu1)*r1
      ZZ1734(jdu1,2,1)=2d0*R(jdu1)*l1
      ZZ1734(jdu1,2,2)=2d0*R(jdu1)*r1
      ZZ2856(jdu1,1,1)=2d0*L(jdu1)*l2
      ZZ2856(jdu1,1,2)=2d0*L(jdu1)*r2
      ZZ2856(jdu1,2,1)=2d0*R(jdu1)*l2
      ZZ2856(jdu1,2,2)=2d0*R(jdu1)*r2
      do jdu2=1,2
      ZZ1728(jdu1,jdu2,1,1)=2d0*L(jdu1)*L(jdu2)
      ZZ1728(jdu1,jdu2,1,2)=2d0*L(jdu1)*R(jdu2)
      ZZ1728(jdu1,jdu2,2,1)=2d0*R(jdu1)*L(jdu2)
      ZZ1728(jdu1,jdu2,2,2)=2d0*R(jdu1)*R(jdu2)
      enddo
      enddo
      endif
      fac=-zmass**2/(cxw*(cone-cxw)) ! this is the HZZ coupling squared

C---setup propagators
      s3456=t4(n3,n4,n5,n6)
      s1734=t4(n1,n7,n3,n4)
      s1756=t4(n1,n7,n5,n6)

      prop17=dcmplx(s(n1,n7))-czmass2
      prop28=dcmplx(s(n2,n8))-czmass2
      prop34=dcmplx(s(n3,n4))-czmass2
      prop56=dcmplx(s(n5,n6))-czmass2
      prop3456=dcmplx(s3456-hmass**2,htheta(s3456)*hmass*hwidth)
      prop1734=dcmplx(s1734-hmass**2,htheta(s1734)*hmass*hwidth)
      prop1756=dcmplx(s1756-hmass**2,htheta(s1756)*hmass*hwidth)
      propWBF=prop17*prop28*prop34*prop56

!       print *,""
!       print *, "zz prop17",1d0/cdabs(prop17)**2
!       print *, "zz prop28",1d0/cdabs(prop28)**2
!       print *, "zz prop34",1d0/cdabs(prop34)**2
!       print *, "zz prop56",1d0/cdabs(prop56)**2
!       print *, "zz prop3456",1d0/cdabs(prop3456)**2


      do h17=1,2
         if (h17.eq.1) then 
            i1=n1
            i7=n7
         elseif (h17.eq.2) then
            i1=n7
            i7=n1
         endif
      do h28=1,2
         if (h28.eq.1) then 
            i2=n2
            i8=n8
         elseif (h28.eq.2) then
            i2=n8
            i8=n2
         endif
      do h34=1,2
         if (h34.eq.1) then 
            i3=n3
            i4=n4
         elseif (h34.eq.2) then
            i3=n4
            i4=n3
         endif
      do h56=1,2
         if (h56.eq.1) then 
            i5=n5
            i6=n6
         elseif (h56.eq.2) then
            i5=n6
            i6=n5
         endif
         
      do jdu1=1,2
      
      do jdu2=1,2


! ! C---s-channel   MARKUS: removed decay
!       ZZHamp(jdu1,jdu2,h17,h28,h34,h56)=
!      & +fac*ZZ1728(jdu1,jdu2,h17,h28)
! !      & *za(i7,i8)*zb(i2,i1)*za(i3,i5)*zb(i6,i4) !-- removed from MCFM
!      & *za(i7,i8)*zb(i2,i1)
!      & /(prop17*prop28)/( zmass/sqrt(cxw*(cone-cxw)) )  !-- only  MCFM PR
! !      & *anomzzamp(i3,i4,i5,i6,s3456,s(i3,i4),s(i5,i6),za,zb) !--added by Fabrizio (anomal. DK)
! !      & *anomzzamp(i7,i1,i8,i2,s3456,s(i7,i1),s(i8,i2),za,zb) !--added by Fabrizio (anomal. PR)
! !      & /(propWBF*prop3456)



C---s-channel
      ZZHamp(jdu1,jdu2,h17,h28,h34,h56)=
     & +fac*ZZ3456(h34,h56)*ZZ1728(jdu1,jdu2,h17,h28)
     & *za(i7,i8)*zb(i2,i1)*za(i3,i5)*zb(i6,i4) !--F
!      & *anomzzamp(i3,i4,i5,i6,s3456,s(i3,i4),s(i5,i6),za,zb) !--F
!      & *anomzzamp(i7,i1,i8,i2,s3456,s(i7,i1),s(i8,i2),za,zb) !--F
     & /(propWBF*prop3456)
C--t-channel
     & +fac*ZZ1734(jdu1,h17,h34)*ZZ2856(jdu2,h28,h56)
     & *za(i7,i3)*zb(i4,i1)*za(i8,i5)*zb(i6,i2) !--F
!      & *anomzzamp(i3,i4,i7,i1,s1734,s(i3,i4),s(i7,i1),za,zb) !--F
!      & *anomzzamp(i5,i6,i8,i2,s1734,s(i5,i6),s(i8,i2),za,zb) !--F
     & /(propWBF*prop1734)
C---u-channel
     & +fac*ZZ2856(jdu1,h17,h56)*ZZ1734(jdu2,h28,h34)
     & *za(i7,i5)*zb(i6,i1)*za(i8,i3)*zb(i4,i2) !--F
!      & *anomzzamp(i5,i6,i7,i1,s1756,s(i5,i6),s(i7,i1),za,zb) !--F
!      & *anomzzamp(i3,i4,i8,i2,s1756,s(i3,i4),s(i8,i2),za,zb) !--F
     & /(propWBF*prop1756)



      enddo
      enddo

      enddo
      enddo
      enddo
      enddo

      return
      end
