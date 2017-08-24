      subroutine amp2currentw(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'pid_pdg.f'
      include 'zacouplejk.f'
      logical isALepton,isANeutrino,isUpTypeLightQuark,isDnTypeQuark,
     & isAnUnknownJet
      integer jdu1,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      double complex zab2,zba2,amp(2,2),rxw,
     & propw34,propw56,propw28,propz3456,
     & gam3456f3(2,2,2),gam3456f4(2,2,2),
     & gam3456f6(2,2,2),gam3456f5(2,2,2),
     & gamv(2,2)
      double precision qn,t3,t4,s134,s128,s345,s347,s356,s156,s456,
     & s278,s346,s567,s3456,s28,s34,s56,
     & q3,q4,q5,q6,l3,l4,l5,l6,r3,r4,r5,r6
      parameter(qn=0d0)
C-----Begin statement functions
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s28=s(i2,i8)
      s34=s(i3,i4)
      s56=s(i5,i6)
      s128=t3(i1,i2,i8)
      s134=t3(i1,i3,i4)
      s156=t3(i1,i5,i6)
      s278=t3(i2,i7,i8)
      s345=t3(i3,i4,i5)
      s346=t3(i3,i4,i6)
      s347=t3(i3,i4,i7)
      s356=t3(i3,i5,i6)
      s456=t3(i4,i5,i6)
      s567=t3(i5,i6,i7)
      s3456=t4(i3,i4,i5,i6)

      propw34=s34-dcmplx(wmass**2,-wmass*wwidth)
      propw56=s56-dcmplx(wmass**2,-wmass*wwidth)
      propw28=s28-dcmplx(wmass**2,-wmass*wwidth)
      propz3456=s3456-dcmplx(zmass**2,-zmass*zwidth)

      ! Find the charges for the 34 line
      if(
     & (isUpTypeLightQuark(-pid_pdg(i4)).or.isAnUnknownJet(pid_pdg(i4)))
     & .and.
     & (isDnTypeQuark(pid_pdg(i3)) .or. isAnUnknownJet(pid_pdg(i3)))
     & ) then
         l3=L(1)
         l4=L(2)
         r3=R(1)
         r4=R(2)
         q3=Q(1)
         q4=Q(2)
      else if(
     & (isANeutrino(-pid_pdg(i4)) .and. isALepton(pid_pdg(i3)))
     & ) then
         l3=le
         l4=ln
         r3=re
         r4=rn
         q3=qe
         q4=qn
      else
         l3=0d0
         l4=0d0
         r3=0d0
         r4=0d0
         q3=0d0
         q4=0d0
         write(6,*)
     &   "Unrecognized id combination in amp2currentw(id3,id4)=",
     &   pid_pdg(i3),pid_pdg(i4)
      endif

      ! Find the charges for the 56 line
      if(
     & (isDnTypeQuark(-pid_pdg(i6)).or.isAnUnknownJet(pid_pdg(i6)))
     & .and.
     & (isUpTypeLightQuark(pid_pdg(i5)).or.isAnUnknownJet(pid_pdg(i5)))
     & ) then
         l5=L(2)
         l6=L(1)
         r5=R(2)
         r6=R(1)
         q5=Q(2)
         q6=Q(1)
      else if(
     & (isALepton(-pid_pdg(i6)) .and. isANeutrino(pid_pdg(i5)))
     & ) then
         l5=ln
         l6=le
         r5=rn
         r6=re
         q5=qn
         q6=qe
      else
         l5=0d0
         l6=0d0
         r5=0d0
         r6=0d0
         q5=0d0
         q6=0d0
         write(6,*)
     &   "Unrecognized id combination in amp2currentw(id5,id6)=",
     &   pid_pdg(i5),pid_pdg(i6)
      endif

      do jdu1=1,2
      gam3456f3(jdu1,1,1)=Q_jk(i1,i7,jdu1)*q3/s3456
     & +L_jk(i1,i7,jdu1)*l3/propZ3456
      gam3456f3(jdu1,1,2)=Q_jk(i1,i7,jdu1)*q3/s3456
     & +L_jk(i1,i7,jdu1)*r3/propZ3456
      gam3456f3(jdu1,2,1)=Q_jk(i1,i7,jdu1)*q3/s3456
     & +R_jk(i1,i7,jdu1)*l3/propZ3456
      gam3456f3(jdu1,2,2)=Q_jk(i1,i7,jdu1)*q3/s3456
     & +R_jk(i1,i7,jdu1)*r3/propZ3456
      gam3456f4(jdu1,1,1)=Q_jk(i1,i7,jdu1)*q4/s3456
     & +L_jk(i1,i7,jdu1)*l4/propZ3456
      gam3456f4(jdu1,1,2)=Q_jk(i1,i7,jdu1)*q4/s3456
     & +L_jk(i1,i7,jdu1)*r4/propZ3456
      gam3456f4(jdu1,2,1)=Q_jk(i1,i7,jdu1)*q4/s3456
     & +R_jk(i1,i7,jdu1)*l4/propZ3456
      gam3456f4(jdu1,2,2)=Q_jk(i1,i7,jdu1)*q4/s3456
     & +R_jk(i1,i7,jdu1)*r4/propZ3456

      gam3456f5(jdu1,1,1)=Q_jk(i1,i7,jdu1)*q5/s3456
     & +L_jk(i1,i7,jdu1)*l5/propZ3456
      gam3456f5(jdu1,1,2)=Q_jk(i1,i7,jdu1)*q5/s3456
     & +L_jk(i1,i7,jdu1)*r5/propZ3456
      gam3456f5(jdu1,2,1)=Q_jk(i1,i7,jdu1)*q5/s3456
     & +R_jk(i1,i7,jdu1)*l5/propZ3456
      gam3456f5(jdu1,2,2)=Q_jk(i1,i7,jdu1)*q5/s3456
     & +R_jk(i1,i7,jdu1)*r5/propZ3456
      gam3456f6(jdu1,1,1)=Q_jk(i1,i7,jdu1)*q6/s3456
     & +L_jk(i1,i7,jdu1)*l6/propZ3456
      gam3456f6(jdu1,1,2)=Q_jk(i1,i7,jdu1)*q6/s3456
     & +L_jk(i1,i7,jdu1)*r6/propZ3456
      gam3456f6(jdu1,2,1)=Q_jk(i1,i7,jdu1)*q6/s3456
     & +R_jk(i1,i7,jdu1)*l6/propZ3456
      gam3456f6(jdu1,2,2)=Q_jk(i1,i7,jdu1)*q6/s3456
     & +R_jk(i1,i7,jdu1)*r6/propZ3456

      gamV(jdu1,1)=Q_jk(i1,i7,jdu1)/s3456
     & +L_jk(i1,i7,jdu1)*rxw/propZ3456
      gamV(jdu1,2)=Q_jk(i1,i7,jdu1)/s3456
     & +R_jk(i1,i7,jdu1)*rxw/propZ3456
      enddo

      p1=i1
      p2=i2
      p3=i3
      p4=i4
      p5=i5
      p6=i6
      p7=i7
      p8=i8

      amp(:,:)=czip
      amp(1,2)=
     & + propw56**(-1)*propw34**(-1)*propw28**(-1)*(
     &    + s128**(-1)*s347**(-1)*cxw**(-1) * (
     & - za(p7,p5)*za(p7,p3)*zb(p7,p4)*zb(p1,p2)*zba2(p6,p1,p2,p8)
     & + za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2,p8)
     &    )
     &    + s278**(-1)*s134**(-1)*cxw**(-1) * (
     & - za(p7,p8)*za(p1,p3)*zb(p1,p6)*zb(p1,p4)*zba2(p2,p7,p8,p5)
     & - za(p7,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5)
     &    )
     &    + gamv(2,1)*s128**(-1) * (
     & - za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)
     & + za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1)
     & + za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2)
     & - za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p2)
     & - 2.D0*za(p7,p5)*zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8)
     & + 2.D0*za(p7,p3)*zb(p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8)
     &    )
     &    + gamv(1,1)*s278**(-1) * (
     & + za(p7,p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)
     & - za(p7,p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1)
     & + za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1)
     & - za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1)
     & + 2.D0*za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,p5)
     & - 2.D0*za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,p8,p3)
     &    )
     & )
     & + propw56**(-1)*propw28**(-1) * (
     &    + 2.D0*s128**(-1) * (
     &       + gam3456f3(2,1,1)*s356**(-1)*(
     & + za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p6)*zba2(p4,p1,p2,p8)
     & - za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p6,p3)*zba2(p4,p1,p2,p8)
     &       )
     &       + gam3456f4(2,1,1)*s456**(-1)*(
     & + za(p7,p3)*za(p5,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)
     & + za(p7,p3)*za(p5,p6)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)
     &       )
     &    )
     &    + 2.D0*s278**(-1) * (
     &       + gam3456f3(1,1,1)*s356**(-1)*(
     & - za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p5,p6)*zba2(p2,p7,p8,p5)
     & + za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p6,p3)*zba2(p2,p7,p8,p3)
     &       )
     &       + gam3456f4(1,1,1)*s456**(-1)*(
     & - za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*zba2(p2,p7,p8,p3)
     & - za(p7,p8)*za(p5,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p3)
     &       )
     &    )
     & )
     & + propw34**(-1)*propw28**(-1)*(
     &    + 2.D0*s128**(-1) * (
     &       + gam3456f6(2,1,1)*s346**(-1)*(
     & + za(p7,p5)*za(p6,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)
     & - za(p7,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)
     &       )
     &       + gam3456f5(2,1,1)*s345**(-1)*(
     & - za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p4)*zba2(p6,p1,p2,p8)
     & - za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2,p8)
     &       )
     &    )
     &    + 2.D0*s278**(-1) * (
     &       + gam3456f6(1,1,1)*s346**(-1)*(
     & - za(p7,p8)*za(p6,p3)*zb(p1,p6)*zb(p6,p4)*zba2(p2,p7,p8,p5)
     & + za(p7,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5)
     &       )
     &       + gam3456f5(1,1,1)*s345**(-1)*(
     & + za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p5,p4)*zba2(p2,p7,p8,p5)
     & + za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p3,p4)*zba2(p2,p7,p8,p3)
     &       )
     &    )
     & )

      amp(2,1)=
     & + propw56**(-1)*propw34**(-1)*propw28**(-1)*(
     &    + s278**(-1)*s156**(-1)*cxw**(-1) * (
     & - za(p7,p8)*za(p1,p5)*zb(p1,p6)*zb(p1,p4)*zba2(p2,p7,p8,p3)
     & + za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*zba2(p2,p7,p8,p3)
     &    )
     &    + s128**(-1)*s567**(-1)*cxw**(-1) * (
     & - za(p7,p5)*za(p7,p3)*zb(p7,p6)*zb(p1,p2)*zba2(p4,p1,p2,p8)
     & - za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p6)*zba2(p4,p1,p2,p8)
     &    )
     &    + gamv(1,1)*s128**(-1) * (
     & - za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)
     & + za(p1,p8)*za(p5,p3)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1)
     & + za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p5,p6,p2)
     & - za(p5,p3)*za(p8,p2)*zb(p1,p2)*zb(p6,p4)*zab2(p7,p3,p4,p2)
     & - 2.D0*za(p7,p5)*zb(p1,p2)*zab2(p3,p5,p6,p4)*zba2(p6,p1,p2,p8)
     & + 2.D0*za(p7,p3)*zb(p1,p2)*zab2(p5,p3,p4,p6)*zba2(p4,p1,p2,p8)
     &    )
     &    + gamv(2,1)*s278**(-1) * (
     & + za(p7,p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p5,p6,p1)
     & - za(p7,p8)*za(p5,p3)*zb(p7,p2)*zb(p6,p4)*zab2(p7,p3,p4,p1)
     & + za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p5,p6,p1)
     & - za(p7,p8)*za(p5,p3)*zb(p6,p4)*zb(p8,p2)*zab2(p8,p3,p4,p1)
     & + 2.D0*za(p7,p8)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zba2(p2,p7,p8,p5)
     & - 2.D0*za(p7,p8)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zba2(p2,p7,p8,p3)
     &    )
     & )
     & + propw56**(-1)*propw28**(-1)*(
     &    2.D0*s128**(-1) * (
     &       + gam3456f3(1,1,1)*s356**(-1)*(
     & + za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p6)*zba2(p4,p1,p2,p8)
     & - za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p6,p3)*zba2(p4,p1,p2,p8)
     &       )
     &       + gam3456f4(1,1,1)*s456**(-1)*(
     & + za(p7,p3)*za(p5,p6)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)
     & + za(p7,p3)*za(p5,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)
     &       )
     &    )
     &    + 2.D0**s278**(-1) * (
     &       + gam3456f3(2,1,1)*s356**(-1)*(
     & - za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p5,p6)*zba2(p2,p7,p8,p5)
     & + za(p7,p8)*za(p5,p3)*zb(p1,p4)*zb(p6,p3)*zba2(p2,p7,p8,p3)
     &       )
     &       + gam3456f4(2,1,1)*s456**(-1)*(
     & - za(p7,p8)*za(p5,p6)*zb(p1,p6)*zb(p6,p4)*zba2(p2,p7,p8,p3)
     & - za(p7,p8)*za(p5,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p3)
     &       )
     &    )
     & )
     & + propw34**(-1)*propw28**(-1)*(
     &    + 2.D0*s128**(-1) * (
     &       + gam3456f6(1,1,1)*s346**(-1)*(
     & + za(p7,p5)*za(p6,p3)*zb(p1,p2)*zb(p6,p4)*zba2(p6,p1,p2,p8)
     & - za(p7,p5)*za(p3,p4)*zb(p1,p2)*zb(p6,p4)*zba2(p4,p1,p2,p8)
     &       )
     &       + gam3456f5(1,1,1)*s345**(-1)*(
     & - za(p7,p5)*za(p5,p3)*zb(p1,p2)*zb(p5,p4)*zba2(p6,p1,p2,p8)
     & - za(p7,p3)*za(p5,p3)*zb(p1,p2)*zb(p3,p4)*zba2(p6,p1,p2,p8)
     &       )
     &    )
     &    + 2.D0*s278**(-1) * (
     &       + gam3456f6(2,1,1)*s346**(-1)*(
     & - za(p7,p8)*za(p6,p3)*zb(p1,p6)*zb(p6,p4)*zba2(p2,p7,p8,p5)
     & + za(p7,p8)*za(p3,p4)*zb(p1,p4)*zb(p6,p4)*zba2(p2,p7,p8,p5)
     &       )
     &       + gam3456f5(2,1,1)*s345**(-1)*(
     & + za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p5,p4)*zba2(p2,p7,p8,p5)
     & + za(p7,p8)*za(p5,p3)*zb(p1,p6)*zb(p3,p4)*zba2(p2,p7,p8,p3)
     &       )
     &    )
     & )

      amp(:,:)=amp(:,:)/cxw**2
      return
      end
