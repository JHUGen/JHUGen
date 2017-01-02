      subroutine ZZHZZamp(n1,n2,n3,n4,n5,n6,n7,n8,za,zb,ZZHamp)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zacouplejk.f'
!      include 'first.f'
      include 'spinzerohiggs_anomcoupl.f'
      double precision t4,s3456,s1734,s1756,htheta,
     & ZZ3456(2,2),ZZ1728(2,2,2,2),ZZ1734(2,2,2),ZZ2856(2,2,2),
     & ZA3456(2,2),ZA1728(2,2,2,2),ZA1734(2,2,2),ZA2856(2,2,2),
     & AZ3456(2,2),AZ1728(2,2,2,2),AZ1734(2,2,2),AZ2856(2,2,2),
     & AA3456(2,2),AA1728(2,2,2,2),AA1734(2,2,2),AA2856(2,2,2),
     & ZZ1756(2,2,2),ZZ2834(2,2,2),
     & ZA1756(2,2,2),ZA2834(2,2,2),
     & AZ1756(2,2,2),AZ2834(2,2,2),
     & AA1756(2,2,2),AA2834(2,2,2)
      double complex ZZHamp(2,2,2,2,2,2),
     & prop34,prop56,prop17,prop28,prop3456,prop1734,prop1756,fac,
     & propX3456,propX1734,propX1756
C---order of indices jdu1,jdu2,h17,h28,h34,h56)
      integer h17,h28,h34,h56,i1,i2,i3,i4,i5,i6,i7,i8,
     & n1,n2,n3,n4,n5,n6,n7,n8,jdu1,jdu2
      double complex Amp_S_PR,Amp_S_DK
      double complex Amp_T_PR,Amp_T_DK
      double complex Amp_U_PR,Amp_U_DK
      double complex anomhzzamp,anomhzaamp,anomhaaamp
C---begin statement functions
      t4(i1,i2,i3,i4)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)
     & +s(i2,i3)+s(i2,i4)+s(i3,i4)
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
c      htheta(s3456)=half+sign(half,s3456)
      htheta(s3456)=one ! propdebug
C---end statement functions
!$omp threadprivate(ZZ3456,ZZ1734,ZZ2856,ZZ1728,fac)


      ZZ3456(1,1)=2d0*l1*l2
      ZZ3456(1,2)=2d0*l1*r2
      ZZ3456(2,1)=2d0*r1*l2
      ZZ3456(2,2)=2d0*r1*r2
      ZA3456(1,:)=2d0*l1*q2
      ZA3456(2,:)=2d0*r1*q2
      AZ3456(:,1)=2d0*q1*l2
      AZ3456(:,2)=2d0*q1*r2
      AA3456(:,:)=2d0*q1*q2
      do jdu1=1,2
      ZZ1734(jdu1,1,1)=2d0*L_jk(n1,n7,jdu1)*l1
      ZZ1734(jdu1,1,2)=2d0*L_jk(n1,n7,jdu1)*r1
      ZZ1734(jdu1,2,1)=2d0*R_jk(n1,n7,jdu1)*l1
      ZZ1734(jdu1,2,2)=2d0*R_jk(n1,n7,jdu1)*r1
      ZA1734(jdu1,1,:)=2d0*L_jk(n1,n7,jdu1)*q1
      ZA1734(jdu1,2,:)=2d0*R_jk(n1,n7,jdu1)*q1
      AZ1734(jdu1,:,1)=2d0*Q_jk(n1,n7,jdu1)*l1
      AZ1734(jdu1,:,2)=2d0*Q_jk(n1,n7,jdu1)*r1
      AA1734(jdu1,:,:)=2d0*Q_jk(n1,n7,jdu1)*q1

      ZZ1756(jdu1,1,1)=2d0*L_jk(n1,n7,jdu1)*l2
      ZZ1756(jdu1,1,2)=2d0*L_jk(n1,n7,jdu1)*r2
      ZZ1756(jdu1,2,1)=2d0*R_jk(n1,n7,jdu1)*l2
      ZZ1756(jdu1,2,2)=2d0*R_jk(n1,n7,jdu1)*r2
      ZA1756(jdu1,1,:)=2d0*L_jk(n1,n7,jdu1)*q2
      ZA1756(jdu1,2,:)=2d0*R_jk(n1,n7,jdu1)*q2
      AZ1756(jdu1,:,1)=2d0*Q_jk(n1,n7,jdu1)*l2
      AZ1756(jdu1,:,2)=2d0*Q_jk(n1,n7,jdu1)*r2
      AA1756(jdu1,:,:)=2d0*Q_jk(n1,n7,jdu1)*q2

      ZZ2834(jdu1,1,1)=2d0*L_jk(n2,n8,jdu1)*l1
      ZZ2834(jdu1,1,2)=2d0*L_jk(n2,n8,jdu1)*r1
      ZZ2834(jdu1,2,1)=2d0*R_jk(n2,n8,jdu1)*l1
      ZZ2834(jdu1,2,2)=2d0*R_jk(n2,n8,jdu1)*r1
      ZA2834(jdu1,1,:)=2d0*L_jk(n2,n8,jdu1)*q1
      ZA2834(jdu1,2,:)=2d0*R_jk(n2,n8,jdu1)*q1
      AZ2834(jdu1,:,1)=2d0*Q_jk(n2,n8,jdu1)*l1
      AZ2834(jdu1,:,2)=2d0*Q_jk(n2,n8,jdu1)*r1
      AA2834(jdu1,:,:)=2d0*Q_jk(n2,n8,jdu1)*q1

      ZZ2856(jdu1,1,1)=2d0*L_jk(n2,n8,jdu1)*l2
      ZZ2856(jdu1,1,2)=2d0*L_jk(n2,n8,jdu1)*r2
      ZZ2856(jdu1,2,1)=2d0*R_jk(n2,n8,jdu1)*l2
      ZZ2856(jdu1,2,2)=2d0*R_jk(n2,n8,jdu1)*r2
      ZA2856(jdu1,1,:)=2d0*L_jk(n2,n8,jdu1)*q2
      ZA2856(jdu1,2,:)=2d0*R_jk(n2,n8,jdu1)*q2
      AZ2856(jdu1,:,1)=2d0*Q_jk(n2,n8,jdu1)*l2
      AZ2856(jdu1,:,2)=2d0*Q_jk(n2,n8,jdu1)*r2
      AA2856(jdu1,:,:)=2d0*Q_jk(n2,n8,jdu1)*q2
      do jdu2=1,2
      ZZ1728(jdu1,jdu2,1,1)=2d0*L_jk(n1,n7,jdu1)*L_jk(n2,n8,jdu2)
      ZZ1728(jdu1,jdu2,1,2)=2d0*L_jk(n1,n7,jdu1)*R_jk(n2,n8,jdu2)
      ZZ1728(jdu1,jdu2,2,1)=2d0*R_jk(n1,n7,jdu1)*L_jk(n2,n8,jdu2)
      ZZ1728(jdu1,jdu2,2,2)=2d0*R_jk(n1,n7,jdu1)*R_jk(n2,n8,jdu2)
      ZA1728(jdu1,jdu2,1,:)=2d0*L_jk(n1,n7,jdu1)*Q_jk(n2,n8,jdu2)
      ZA1728(jdu1,jdu2,2,:)=2d0*R_jk(n1,n7,jdu1)*Q_jk(n2,n8,jdu2)
      AZ1728(jdu1,jdu2,:,1)=2d0*Q_jk(n1,n7,jdu1)*L_jk(n2,n8,jdu2)
      AZ1728(jdu1,jdu2,:,2)=2d0*Q_jk(n1,n7,jdu1)*R_jk(n2,n8,jdu2)
      AA1728(jdu1,jdu2,:,:)=2d0*Q_jk(n1,n7,jdu1)*Q_jk(n2,n8,jdu2)
      enddo
      enddo

      !fac=-zmass**2/(cxw*(cone-cxw))
      fac=-czmass2/(cxw*(cone-cxw))
     & *(4d0*czmass2/vevsq*cxw*(1-cxw)/esq)

C---setup propagators
      s3456=t4(n3,n4,n5,n6)
      s1734=t4(n1,n7,n3,n4)
      s1756=t4(n1,n7,n5,n6)

      prop17=dcmplx(s(n1,n7))-dcmplx(zmass**2,-zmass*zwidth)
      prop28=dcmplx(s(n2,n8))-dcmplx(zmass**2,-zmass*zwidth)
      prop34=dcmplx(s(n3,n4))-dcmplx(zmass**2,-zmass*zwidth)
      prop56=dcmplx(s(n5,n6))-dcmplx(zmass**2,-zmass*zwidth)
      prop3456=dcmplx(s3456-hmass**2,htheta(s3456)*hmass*hwidth)
      prop1734=dcmplx(s1734-hmass**2,htheta(s1734)*hmass*hwidth)
      prop1756=dcmplx(s1756-hmass**2,htheta(s1756)*hmass*hwidth)

      propX3456=dcmplx(s3456-h2mass**2,htheta(s3456)*h2mass*h2width)
      propX1734=dcmplx(s1734-h2mass**2,htheta(s1734)*h2mass*h2width)
      propX1756=dcmplx(s1756-h2mass**2,htheta(s1756)*h2mass*h2width)
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

C-- MARKUS: this is the old (original) MCFM code
! C---s-channel
!       ZZHamp(jdu1,jdu2,h17,h28,h34,h56)=
!      & +fac*ZZ3456(h34,h56)*ZZ1728(jdu1,jdu2,h17,h28)
!      & *za(i7,i8)*zb(i2,i1)*za(i3,i5)*zb(i6,i4)
!      & /(propZZZZ*prop3456)
! C---t-channel
!      & +fac*ZZ1734(jdu1,h17,h34)*ZZ2856(jdu2,h28,h56)
!      & *za(i7,i3)*zb(i4,i1)*za(i8,i5)*zb(i6,i2)
!      & /(propZZZZ*prop1734)
! C---u-channel
!      & +fac*ZZ2856(jdu1,h17,h56)*ZZ1734(jdu2,h28,h34)
!      & *za(i7,i5)*zb(i6,i1)*za(i8,i3)*zb(i4,i2)
!      & /(propZZZZ*prop1756)

      ZZHamp(jdu1,jdu2,h17,h28,h34,h56)=czip

      Amp_S_PR=czip
      Amp_S_DK=czip
      Amp_T_PR=czip
      Amp_T_DK=czip
      Amp_U_PR=czip
      Amp_U_DK=czip
      if( hmass.ge.zip ) then

        if(channeltoggle_stu.ne.1) then
          if( AnomalCouplPR.eq.1 ) then
      Amp_S_PR=Amp_S_PR
     & -anomhzzamp(i7,i1,i8,i2,1,s3456,s(i7,i1),s(i8,i2),za,zb)
     & /(prop17*prop28)*ZZ1728(jdu1,jdu2,h17,h28)
     & +anomhzaamp(i7,i1,i8,i2,1,s3456,s(i7,i1),s(i8,i2),za,zb)
     & /(prop17*s(n2,n8))*ZA1728(jdu1,jdu2,h17,h28)
     & +anomhzaamp(i8,i2,i7,i1,1,s3456,s(i8,i2),s(i7,i1),za,zb)
     & /(s(n1,n7)*prop28)*AZ1728(jdu1,jdu2,h17,h28)
     & -anomhaaamp(i7,i1,i8,i2,1,s3456,s(i7,i1),s(i8,i2),za,zb)
     & /(s(n1,n7)*s(n2,n8))*AA1728(jdu1,jdu2,h17,h28)
          else
      Amp_S_PR=Amp_S_PR
     & +za(i7,i8)*zb(i2,i1)/(prop17*prop28)*ZZ1728(jdu1,jdu2,h17,h28)
          endif
          if( AnomalCouplDK.eq.1 ) then
      Amp_S_DK=Amp_S_DK
     & -anomhzzamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*prop56)*ZZ3456(h34,h56)
     & +anomhzaamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*s(n5,n6))*ZA3456(h34,h56)
     & +anomhzaamp(i5,i6,i3,i4,1,s3456,s(i5,i6),s(i3,i4),za,zb)
     & /(s(n3,n4)*prop56)*AZ3456(h34,h56)
     & -anomhaaamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(s(n3,n4)*s(n5,n6))*AA3456(h34,h56)
          else
      Amp_S_DK=Amp_S_DK
     & +za(i3,i5)*zb(i6,i4)/(prop34*prop56)*ZZ3456(h34,h56)
          endif
        endif

        if(channeltoggle_stu.ne.0) then
          if( AnomalCouplPR.eq.1 ) then
      Amp_T_PR=Amp_T_PR
     & -anomhzzamp(i3,i4,i7,i1,1,s1734,s(i3,i4),s(i7,i1),za,zb)
     & /(prop17*prop34)*ZZ1734(jdu1,h17,h34)
     & +anomhzaamp(i3,i4,i7,i1,1,s1734,s(i3,i4),s(i7,i1),za,zb)
     & /(s(n1,n7)*prop34)*AZ1734(jdu1,h17,h34)
     & +anomhzaamp(i7,i1,i3,i4,1,s1734,s(i7,i1),s(i3,i4),za,zb)
     & /(prop17*s(n3,n4))*ZA1734(jdu1,h17,h34)
     & -anomhaaamp(i3,i4,i7,i1,1,s1734,s(i3,i4),s(i7,i1),za,zb)
     & /(s(n1,n7)*s(n3,n4))*AA1734(jdu1,h17,h34)
      Amp_U_PR=Amp_U_PR
     & -anomhzzamp(i5,i6,i7,i1,1,s1756,s(i5,i6),s(i7,i1),za,zb)
     & /(prop17*prop56)*ZZ1756(jdu1,h17,h56)
     & +anomhzaamp(i5,i6,i7,i1,1,s1756,s(i5,i6),s(i7,i1),za,zb)
     & /(s(n1,n7)*prop56)*AZ1756(jdu1,h17,h56)
     & +anomhzaamp(i7,i1,i5,i6,1,s1756,s(i7,i1),s(i5,i6),za,zb)
     & /(prop17*s(n5,n6))*ZA1756(jdu1,h17,h56)
     & -anomhaaamp(i5,i6,i7,i1,1,s1756,s(i5,i6),s(i7,i1),za,zb)
     & /(s(n1,n7)*s(n5,n6))*AA1756(jdu1,h17,h56)
          else
      Amp_T_PR=Amp_T_PR
     & +za(i7,i3)*zb(i4,i1)/(prop17*prop34)*ZZ1734(jdu1,h17,h34)
      Amp_U_PR=Amp_U_PR
     & +za(i7,i5)*zb(i6,i1)/(prop17*prop56)*ZZ1734(jdu1,h17,h56)
          endif
          if( AnomalCouplDK.eq.1 ) then
      Amp_T_DK=Amp_T_DK
     & -anomhzzamp(i5,i6,i8,i2,1,s1734,s(i5,i6),s(i8,i2),za,zb)
     & /(prop28*prop56)*ZZ2856(jdu2,h28,h56)
     & +anomhzaamp(i5,i6,i8,i2,1,s1734,s(i5,i6),s(i8,i2),za,zb)
     & /(s(n2,n8)*prop56)*AZ2856(jdu2,h28,h56)
     & +anomhzaamp(i8,i2,i5,i6,1,s1734,s(i8,i2),s(i5,i6),za,zb)
     & /(prop28*s(n5,n6))*ZA2856(jdu2,h28,h56)
     & -anomhaaamp(i5,i6,i8,i2,1,s1734,s(i5,i6),s(i8,i2),za,zb)
     & /(s(n2,n8)*s(n5,n6))*AA2856(jdu2,h28,h56)
      Amp_U_DK=Amp_U_DK
     & -anomhzzamp(i3,i4,i8,i2,1,s1756,s(i3,i4),s(i8,i2),za,zb)
     & /(prop28*prop34)*ZZ2834(jdu2,h28,h34)
     & +anomhzaamp(i3,i4,i8,i2,1,s1756,s(i3,i4),s(i8,i2),za,zb)
     & /(s(n2,n8)*prop34)*AZ2834(jdu2,h28,h34)
     & +anomhzaamp(i8,i2,i3,i4,1,s1756,s(i8,i2),s(i3,i4),za,zb)
     & /(prop28*s(n3,n4))*ZA2834(jdu2,h28,h34)
     & -anomhaaamp(i3,i4,i8,i2,1,s1756,s(i3,i4),s(i8,i2),za,zb)
     & /(s(n2,n8)*s(n3,n4))*AA2834(jdu2,h28,h34)
          else
      Amp_T_DK=Amp_T_DK
     & +za(i8,i5)*zb(i6,i2)/(prop28*prop56)*ZZ2856(jdu2,h28,h56)
      Amp_U_DK=Amp_U_DK
     & +za(i8,i3)*zb(i4,i2)/(prop28*prop34)*ZZ2834(jdu2,h28,h34)
          endif
        endif

      endif
      ZZHamp(jdu1,jdu2,h17,h28,h34,h56)=
     & ZZHamp(jdu1,jdu2,h17,h28,h34,h56)
     & +fac*(
C---s-channel
     & Amp_S_DK*Amp_S_PR/prop3456
C---t-channel
     & +Amp_T_DK*Amp_T_PR/prop1734
C---u-channel
     & +Amp_U_DK*Amp_U_PR/prop1756
     & )

      Amp_S_PR=czip
      Amp_S_DK=czip
      Amp_T_PR=czip
      Amp_T_DK=czip
      Amp_U_PR=czip
      Amp_U_DK=czip
      if( h2mass.ge.zip ) then

        if(channeltoggle_stu.ne.1) then
          if( AnomalCouplPR.eq.1 ) then
      Amp_S_PR=Amp_S_PR
     & -anomhzzamp(i7,i1,i8,i2,2,s3456,s(i7,i1),s(i8,i2),za,zb)
     & /(prop17*prop28)*ZZ1728(jdu1,jdu2,h17,h28)
     & +anomhzaamp(i7,i1,i8,i2,2,s3456,s(i7,i1),s(i8,i2),za,zb)
     & /(prop17*s(n2,n8))*ZA1728(jdu1,jdu2,h17,h28)
     & +anomhzaamp(i8,i2,i7,i1,2,s3456,s(i8,i2),s(i7,i1),za,zb)
     & /(s(n1,n7)*prop28)*AZ1728(jdu1,jdu2,h17,h28)
     & -anomhaaamp(i7,i1,i8,i2,2,s3456,s(i7,i1),s(i8,i2),za,zb)
     & /(s(n1,n7)*s(n2,n8))*AA1728(jdu1,jdu2,h17,h28)
          else
      Amp_S_PR=Amp_S_PR
     & +za(i7,i8)*zb(i2,i1)/(prop17*prop28)*ZZ1728(jdu1,jdu2,h17,h28)
          endif
          if( AnomalCouplDK.eq.1 ) then
      Amp_S_DK=Amp_S_DK
     & -anomhzzamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*prop56)*ZZ3456(h34,h56)
     & +anomhzaamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(prop34*s(n5,n6))*ZA3456(h34,h56)
     & +anomhzaamp(i5,i6,i3,i4,2,s3456,s(i5,i6),s(i3,i4),za,zb)
     & /(s(n3,n4)*prop56)*AZ3456(h34,h56)
     & -anomhaaamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
     & /(s(n3,n4)*s(n5,n6))*AA3456(h34,h56)
          else
      Amp_S_DK=Amp_S_DK
     & +za(i3,i5)*zb(i6,i4)/(prop34*prop56)*ZZ3456(h34,h56)
          endif
        endif

        if(channeltoggle_stu.ne.0) then
          if( AnomalCouplPR.eq.1 ) then
      Amp_T_PR=Amp_T_PR
     & -anomhzzamp(i3,i4,i7,i1,2,s1734,s(i3,i4),s(i7,i1),za,zb)
     & /(prop17*prop34)*ZZ1734(jdu1,h17,h34)
     & +anomhzaamp(i3,i4,i7,i1,2,s1734,s(i3,i4),s(i7,i1),za,zb)
     & /(s(n1,n7)*prop34)*AZ1734(jdu1,h17,h34)
     & +anomhzaamp(i7,i1,i3,i4,2,s1734,s(i7,i1),s(i3,i4),za,zb)
     & /(prop17*s(n3,n4))*ZA1734(jdu1,h17,h34)
     & -anomhaaamp(i3,i4,i7,i1,2,s1734,s(i3,i4),s(i7,i1),za,zb)
     & /(s(n1,n7)*s(n3,n4))*AA1734(jdu1,h17,h34)
      Amp_U_PR=Amp_U_PR
     & -anomhzzamp(i5,i6,i7,i1,2,s1756,s(i5,i6),s(i7,i1),za,zb)
     & /(prop17*prop56)*ZZ1756(jdu1,h17,h56)
     & +anomhzaamp(i5,i6,i7,i1,2,s1756,s(i5,i6),s(i7,i1),za,zb)
     & /(s(n1,n7)*prop56)*AZ1756(jdu1,h17,h56)
     & +anomhzaamp(i7,i1,i5,i6,2,s1756,s(i7,i1),s(i5,i6),za,zb)
     & /(prop17*s(n5,n6))*ZA1756(jdu1,h17,h56)
     & -anomhaaamp(i5,i6,i7,i1,2,s1756,s(i5,i6),s(i7,i1),za,zb)
     & /(s(n1,n7)*s(n5,n6))*AA1756(jdu1,h17,h56)
          else
      Amp_T_PR=Amp_T_PR
     & +za(i7,i3)*zb(i4,i1)/(prop17*prop34)*ZZ1734(jdu1,h17,h34)
      Amp_U_PR=Amp_U_PR
     & +za(i7,i5)*zb(i6,i1)/(prop17*prop56)*ZZ1734(jdu1,h17,h56)
          endif
          if( AnomalCouplDK.eq.1 ) then
      Amp_T_DK=Amp_T_DK
     & -anomhzzamp(i5,i6,i8,i2,2,s1734,s(i5,i6),s(i8,i2),za,zb)
     & /(prop28*prop56)*ZZ2856(jdu2,h28,h56)
     & +anomhzaamp(i5,i6,i8,i2,2,s1734,s(i5,i6),s(i8,i2),za,zb)
     & /(s(n2,n8)*prop56)*AZ2856(jdu2,h28,h56)
     & +anomhzaamp(i8,i2,i5,i6,2,s1734,s(i8,i2),s(i5,i6),za,zb)
     & /(prop28*s(n5,n6))*ZA2856(jdu2,h28,h56)
     & -anomhaaamp(i5,i6,i8,i2,2,s1734,s(i5,i6),s(i8,i2),za,zb)
     & /(s(n2,n8)*s(n5,n6))*AA2856(jdu2,h28,h56)
      Amp_U_DK=Amp_U_DK
     & -anomhzzamp(i3,i4,i8,i2,2,s1756,s(i3,i4),s(i8,i2),za,zb)
     & /(prop28*prop34)*ZZ2834(jdu2,h28,h34)
     & +anomhzaamp(i3,i4,i8,i2,2,s1756,s(i3,i4),s(i8,i2),za,zb)
     & /(s(n2,n8)*prop34)*AZ2834(jdu2,h28,h34)
     & +anomhzaamp(i8,i2,i3,i4,2,s1756,s(i8,i2),s(i3,i4),za,zb)
     & /(prop28*s(n3,n4))*ZA2834(jdu2,h28,h34)
     & -anomhaaamp(i3,i4,i8,i2,2,s1756,s(i3,i4),s(i8,i2),za,zb)
     & /(s(n2,n8)*s(n3,n4))*AA2834(jdu2,h28,h34)
          else
      Amp_T_DK=Amp_T_DK
     & +za(i8,i5)*zb(i6,i2)/(prop28*prop56)*ZZ2856(jdu2,h28,h56)
      Amp_U_DK=Amp_U_DK
     & +za(i8,i3)*zb(i4,i2)/(prop28*prop34)*ZZ2834(jdu2,h28,h34)
          endif
        endif

      endif
      ZZHamp(jdu1,jdu2,h17,h28,h34,h56)=
     & ZZHamp(jdu1,jdu2,h17,h28,h34,h56)
     & +fac*(
C---s-channel
     & Amp_S_DK*Amp_S_PR/propX3456
C---t-channel
     & +Amp_T_DK*Amp_T_PR/propX1734
C---u-channel
     & +Amp_U_DK*Amp_U_PR/propX1756
     & )


      enddo
      enddo

      enddo
      enddo
      enddo
      enddo

      return
      end

