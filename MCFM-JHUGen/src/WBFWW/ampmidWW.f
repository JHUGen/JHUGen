      subroutine ampmidWW(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cmplxmass.f'
      include 'masses.f'
      include 'runstring.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'pid_pdg.f'
      include 'zacouplejk.f'
      include 'WWbits.f'
      include 'spinzerohiggs_anomcoupl.f'
      logical isALepton,isANeutrino,isUpTypeLightQuark,isDnTypeQuark,
     & isAnUnknownJet
      integer i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      double complex zab2,zba2,amp,sqwmass,rxw,facWWHWW,
     & propw34,propw56,propw28,propw17,anomhwwamp,
     & Amp_S_PR,Amp_S_DK,Amp_T_PR,Amp_T_DK,
     & propz3456,propa3456,proph3456,propX3456,
     & propz1347,propa1347,proph1347,propX1347,
     & gam3456f4,gam3456f3,gam3456f5,gam3456f6,
     & gam1347f1,gam1347f3,gam1347f4,gam1347f7,
     & gam1347f2,gam1347f5,gam1347f6,gam1347f8
      double precision t3,t4,s34,s56,s17,s28,s134,s137,s147,
     & s256,s258,s268,s456,s345,s347,s356,s346,s568,s3456,s1347,s1567,
     & twop17Dp3456,twop28Dp3456,twop34Dp3456,twop56Dp3456,
     & twop17Dp1347,twop28Dp1347,twop34Dp1347,twop56Dp1347,
     & twop17Dp2568,twop28Dp2568,twop34Dp2568,twop56Dp2568,
     & qn,
     & q3,q4,q5,q6,l3,l4,l5,l6,r3,r4,r5,r6,
     & line8q,line2q,line7q,line1q,
     & line8l,line2l,line7l,line1l,
     & line8r,line2r,line7r,line1r
      parameter(qn=0d0)
C-----Begin statement functions
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

c--- special fix for Madgraph check
      if (index(runstring,'mad') .gt. 0) then
        sqwmass=dcmplx(wmass**2,0d0)
      else
        sqwmass=cwmass2
      endif
      facWWHWW = sqwmass*cxw**(-3)*(4d0*cwmass2/vevsq*cxw/esq)

      rxw=sqrt((cone-cxw)/cxw)
      s3456=t4(i3,i4,i5,i6)
      s1347=t4(i1,i3,i4,i7)
      s1567=t4(i1,i5,i6,i7)
      s137=t3(i1,i3,i7)
      s134=t3(i1,i3,i4)
      s147=t3(i1,i4,i7)
      s258=t3(i2,i5,i8)
      s268=t3(i2,i6,i8)
      s345=t3(i3,i4,i5)
      s346=t3(i3,i4,i6)
      s347=t3(i3,i4,i7)
      s256=t3(i2,i5,i6)
      s356=t3(i3,i5,i6)
      s456=t3(i4,i5,i6)
      s568=t3(i5,i6,i8)
      s34=s(i3,i4)
      s56=s(i5,i6)
      s17=s(i1,i7)
      s28=s(i2,i8)

      twop17Dp3456=s28-s17-s3456
      twop28Dp3456=s17-s28-s3456
      twop34Dp3456=s34-s56+s3456
      twop56Dp3456=s56-s34+s3456

      twop17Dp1347=s17-s34+s1347
      twop28Dp1347=s56-s28-s1347
      twop34Dp1347=s34-s17+s1347
      twop56Dp1347=s28-s56-s1347

      twop17Dp2568=-twop17Dp1347
      twop28Dp2568=-twop28Dp1347
      twop34Dp2568=-twop34Dp1347
      twop56Dp2568=-twop56Dp1347

      propw34=s34-dcmplx(wmass**2,-wmass*wwidth)
      propw56=s56-dcmplx(wmass**2,-wmass*wwidth)
      propw17=s17-dcmplx(wmass**2,-wmass*wwidth)
      propw28=s28-dcmplx(wmass**2,-wmass*wwidth)

      propz3456=s3456-dcmplx(zmass**2,-zmass*zwidth)
      propa3456=s3456
      proph3456=dcmplx(s3456-hmass**2,hmass*hwidth)
      propX3456=dcmplx(s3456-h2mass**2,h2mass*h2width)

      propz1347=s1347-dcmplx(zmass**2,-zmass*zwidth)
      propa1347=s1347
      proph1347=dcmplx(s1347-hmass**2,hmass*hwidth)
      propX1347=dcmplx(s1347-h2mass**2,h2mass*h2width)

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
     &   "Unrecognized id combination in ampmidWW(id3,id4)=",
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
     &   "Unrecognized id combination in ampmidWW(id5,id6)=",
     &   pid_pdg(i5),pid_pdg(i6)
      endif
c      ! Find the charges for the 71 line
c      if(
c     & (isDnTypeQuark(-pid_pdg(i1)).or.isAnUnknownJet(pid_pdg(i1)))
c     & .and.
c     & (isUpTypeLightQuark(pid_pdg(i7)).or.isAnUnknownJet(pid_pdg(i7)))
c     & ) then
c         line7l=L(2)
c         line1l=L(1)
c         line7r=R(2)
c         line1r=R(1)
c         line7q=Q(2)
c         line1q=Q(1)
c      else if(
c     & (isALepton(-pid_pdg(i1)) .and. isANeutrino(pid_pdg(i7)))
c     & ) then
c         line7l=ln
c         line1l=le
c         line7r=rn
c         line1r=re
c         line7q=qn
c         line1q=qe
c      else
c         line7l=0d0
c         line1l=0d0
c         line7r=0d0
c         line1r=0d0
c         line7q=0d0
c         line1q=0d0
c         write(6,*)
c     &   "Unrecognized id combination in ampmidWW(id7,id1)=",
c     &   pid_pdg(i7),pid_pdg(i1)
c      endif
c      ! Find the charges for the 82 line
c      if(
c     & (isUpTypeLightQuark(-pid_pdg(i2)).or.isAnUnknownJet(pid_pdg(i2)))
c     & .and.
c     & (isDnTypeQuark(pid_pdg(i8)) .or. isAnUnknownJet(pid_pdg(i8)))
c     & ) then
c         line8l=L(1)
c         line2l=L(2)
c         line8r=R(1)
c         line2r=R(2)
c         line8q=Q(1)
c         line2q=Q(2)
c      else if(
c     & (isANeutrino(-pid_pdg(i2)) .and. isALepton(pid_pdg(i8)))
c     & ) then
c         line8l=le
c         line2l=ln
c         line8r=re
c         line2r=rn
c         line8q=qe
c         line2q=qn
c      else
c         line8l=0d0
c         line2l=0d0
c         line8r=0d0
c         line2r=0d0
c         line8q=0d0
c         line2q=0d0
c         write(6,*)
c     &   "Unrecognized id combination in ampmidWW(id8,id2)=",
c     &   pid_pdg(i8),pid_pdg(i2)
c      endif

      gam3456f3=q3/propa3456+l3*rxw/propz3456
      gam3456f4=q4/propa3456+l4*rxw/propz3456
      gam3456f5=q5/propa3456+l5*rxw/propz3456
      gam3456f6=q6/propa3456+l6*rxw/propz3456

      gam1347f1=line1q/propa1347+line1l*rxw/propz1347
      gam1347f2=line2q/propa1347+line2l*rxw/propz1347
      gam1347f3=q3/propa1347+l3*rxw/propz1347
      gam1347f4=q4/propa1347+l4*rxw/propz1347
      gam1347f3=q5/propa1347+l5*rxw/propz1347
      gam1347f4=q6/propa1347+l6*rxw/propz1347
      gam1347f7=line7q/propa1347+line7l*rxw/propz1347
      gam1347f8=line8q/propa1347+line8l*rxw/propz1347

      p1=i1
      p2=i2
      p3=i3
      p4=i4
      p5=i5
      p6=i6
      p7=i7
      p8=i8
      amp= + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propz3456**(-1) * ( 1.D0/4.D0*za(p3,p5)*za(p7,p8)*
     &    zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-3)*twop28Dp3456*
     &    twop56Dp3456 - 1.D0/4.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p6)*czmass2**(-1)*cxw**(-3)*twop28Dp3456*twop34Dp3456 - 1.D0/
     &    4.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*
     &    cxw**(-3)*twop17Dp3456*twop56Dp3456 + 1.D0/4.D0*za(p3,p5)*za(
     &    p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-3)*
     &    twop17Dp3456*twop34Dp3456 - 1.D0/4.D0*za(p3,p5)*za(p7,p8)*zb(
     &    p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*twop28Dp3456*
     &    twop56Dp3456 + 1.D0/4.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p6)*czmass2**(-1)*cxw**(-2)*twop28Dp3456*twop34Dp3456 + 1.D0/
     &    4.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*
     &    cxw**(-2)*twop17Dp3456*twop56Dp3456 - 1.D0/4.D0*za(p3,p5)*za(
     &    p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456*twop34Dp3456 )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propz3456**(-1) * ( za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*cxw**(-3)*s1567 - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(
     &    p4,p6)*cxw**(-3)*s1347 - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p6)*cxw**(-2)*s1567 + za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)
     &    *cxw**(-2)*s1347 + za(p3,p5)*zb(p4,p6)*zab2(p7,p2,p8,p1)*
     &    zab2(p8,p3,p4,p2)*cxw**(-3) - za(p3,p5)*zb(p4,p6)*zab2(p7,p2,
     &    p8,p1)*zab2(p8,p3,p4,p2)*cxw**(-2) - za(p3,p5)*zb(p4,p6)*
     &    zab2(p7,p2,p8,p1)*zab2(p8,p5,p6,p2)*cxw**(-3) + za(p3,p5)*zb(
     &    p4,p6)*zab2(p7,p2,p8,p1)*zab2(p8,p5,p6,p2)*cxw**(-2) - za(p3,
     &    p5)*zb(p4,p6)*zab2(p7,p3,p4,p1)*zab2(p8,p1,p7,p2)*cxw**(-3)
     &     + za(p3,p5)*zb(p4,p6)*zab2(p7,p3,p4,p1)*zab2(p8,p1,p7,p2)*
     &    cxw**(-2) + za(p3,p5)*zb(p4,p6)*zab2(p7,p5,p6,p1)*zab2(p8,p1,
     &    p7,p2)*cxw**(-3) - za(p3,p5)*zb(p4,p6)*zab2(p7,p5,p6,p1)*
     &    zab2(p8,p1,p7,p2)*cxw**(-2) - 2.D0*za(p3,p7)*zb(p1,p4)*zab2(
     &    p5,p3,p4,p6)*zab2(p8,p1,p7,p2)*cxw**(-3) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propz3456**(-1) * ( 2.D0*za(p3,p7)*zb(p1,p4)*zab2(
     &    p5,p3,p4,p6)*zab2(p8,p1,p7,p2)*cxw**(-2) + 2.D0*za(p3,p8)*zb(
     &    p2,p4)*zab2(p5,p3,p4,p6)*zab2(p7,p2,p8,p1)*cxw**(-3) - 2.D0*
     &    za(p3,p8)*zb(p2,p4)*zab2(p5,p3,p4,p6)*zab2(p7,p2,p8,p1)*
     &    cxw**(-2) + 2.D0*za(p5,p7)*zb(p1,p6)*zab2(p3,p5,p6,p4)*zab2(
     &    p8,p1,p7,p2)*cxw**(-3) - 2.D0*za(p5,p7)*zb(p1,p6)*zab2(p3,p5,
     &    p6,p4)*zab2(p8,p1,p7,p2)*cxw**(-2) - 2.D0*za(p5,p8)*zb(p2,p6)
     &    *zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1)*cxw**(-3) + 2.D0*za(p5,
     &    p8)*zb(p2,p6)*zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1)*cxw**(-2)
     &     - za(p7,p8)*zb(p1,p2)*zab2(p3,p1,p7,p4)*zab2(p5,p3,p4,p6)*
     &    cxw**(-3) + za(p7,p8)*zb(p1,p2)*zab2(p3,p1,p7,p4)*zab2(p5,p3,
     &    p4,p6)*cxw**(-2) + za(p7,p8)*zb(p1,p2)*zab2(p3,p2,p8,p4)*
     &    zab2(p5,p3,p4,p6)*cxw**(-3) - za(p7,p8)*zb(p1,p2)*zab2(p3,p2,
     &    p8,p4)*zab2(p5,p3,p4,p6)*cxw**(-2) + za(p7,p8)*zb(p1,p2)*
     &    zab2(p3,p5,p6,p4)*zab2(p5,p1,p7,p6)*cxw**(-3) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propz3456**(-1) * (  - za(p7,p8)*zb(p1,p2)*zab2(p3
     &    ,p5,p6,p4)*zab2(p5,p1,p7,p6)*cxw**(-2) - za(p7,p8)*zb(p1,p2)*
     &    zab2(p3,p5,p6,p4)*zab2(p5,p2,p8,p6)*cxw**(-3) + za(p7,p8)*zb(
     &    p1,p2)*zab2(p3,p5,p6,p4)*zab2(p5,p2,p8,p6)*cxw**(-2) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propa3456**(-1) * ( za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*cxw**(-2)*s1567 - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(
     &    p4,p6)*cxw**(-2)*s1347 + za(p3,p5)*zb(p4,p6)*zab2(p7,p2,p8,p1
     &    )*zab2(p8,p3,p4,p2)*cxw**(-2) - za(p3,p5)*zb(p4,p6)*zab2(p7,
     &    p2,p8,p1)*zab2(p8,p5,p6,p2)*cxw**(-2) - za(p3,p5)*zb(p4,p6)*
     &    zab2(p7,p3,p4,p1)*zab2(p8,p1,p7,p2)*cxw**(-2) + za(p3,p5)*zb(
     &    p4,p6)*zab2(p7,p5,p6,p1)*zab2(p8,p1,p7,p2)*cxw**(-2) - 2.D0*
     &    za(p3,p7)*zb(p1,p4)*zab2(p5,p3,p4,p6)*zab2(p8,p1,p7,p2)*
     &    cxw**(-2) + 2.D0*za(p3,p8)*zb(p2,p4)*zab2(p5,p3,p4,p6)*zab2(
     &    p7,p2,p8,p1)*cxw**(-2) + 2.D0*za(p5,p7)*zb(p1,p6)*zab2(p3,p5,
     &    p6,p4)*zab2(p8,p1,p7,p2)*cxw**(-2) - 2.D0*za(p5,p8)*zb(p2,p6)
     &    *zab2(p3,p5,p6,p4)*zab2(p7,p2,p8,p1)*cxw**(-2) - za(p7,p8)*
     &    zb(p1,p2)*zab2(p3,p1,p7,p4)*zab2(p5,p3,p4,p6)*cxw**(-2) + za(
     &    p7,p8)*zb(p1,p2)*zab2(p3,p2,p8,p4)*zab2(p5,p3,p4,p6)*
     &    cxw**(-2) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propa3456**(-1) * ( za(p7,p8)*zb(p1,p2)*zab2(p3,p5
     &    ,p6,p4)*zab2(p5,p1,p7,p6)*cxw**(-2) - za(p7,p8)*zb(p1,p2)*
     &    zab2(p3,p5,p6,p4)*zab2(p5,p2,p8,p6)*cxw**(-2) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propa3456**(-1) * (  - 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*czmass2**(-1)*cxw**(-2)*twop28Dp3456*(q5-q6)
     &    + 1.D0/2.D0
     &    *za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*
     &    cxw**(-2)*twop17Dp3456*(q5-q6) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & propa3456**(-1) * ( 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(
     &    p4,p6)*czmass2**(-1)*cxw**(-2)*twop28Dp3456*(q4-q3)
     &    - 1.D0/2.D0*
     &    za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*
     &    cxw**(-2)*twop17Dp3456*(q4-q3) )

      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1) * (  - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    cxw**(-3) - za(p3,p7)*za(p5,p8)*zb(p1,p4)*zb(p2,p6)*cxw**(-3)
     &     + 2.D0*za(p3,p8)*za(p5,p7)*zb(p1,p6)*zb(p2,p4)*cxw**(-3) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s346**(-1)*s258**(-1) * (  - za(p5,p8)*zb(p4,p6)*zba2(p1,p4,p6,
     &    p3)*zba2(p2,p5,p8,p7)*cxw**(-3) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s345**(-1)*s268**(-1) * ( za(p3,p5)*zb(p2,p6)*zba2(p1,p2,p6,p8)*
     &    zba2(p4,p3,p5,p7)*cxw**(-3) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s456**(-1)*s137**(-1) * ( + za(p3,p7)*zb(p4,p6)*zba2(p1,p3,p7,
     &    p8)*zba2(p2,p4,p6,p5)*cxw**(-3) )
      amp = amp + propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s356**(-1)*s147**(-1) * ( - za(p3,p5)*zb(p1,p4)*zba2(p2,p1,p4,
     &    p7)*zba2(p6,p3,p5,p8)*cxw**(-3) )

      amp = amp + gam3456f3*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s356**(-1) * (  - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*zba2(
     &    p6,p3,p5,p1)*cxw**(-2) + za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p2,
     &    p4)*zba2(p6,p3,p5,p2)*cxw**(-2) + za(p3,p5)*za(p7,p8)*zb(p1,
     &    p2)*zb(p4,p7)*zba2(p6,p3,p5,p7)*cxw**(-2) - za(p3,p5)*za(p7,
     &    p8)*zb(p1,p2)*zb(p4,p8)*zba2(p6,p3,p5,p8)*cxw**(-2) - 2.D0*
     &    za(p3,p5)*zb(p1,p4)*zab2(p8,p1,p7,p2)*zba2(p6,p3,p5,p7)*
     &    cxw**(-2) + 2.D0*za(p3,p5)*zb(p2,p4)*zab2(p7,p2,p8,p1)*zba2(
     &    p6,p3,p5,p8)*cxw**(-2) )
      amp = amp + gam3456f3*propw17**(-1)*propw28**(-1)*propw56**(-1)
     &  * ( 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 - 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      amp = amp + gam3456f4*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s456**(-1) * ( za(p1,p3)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p1,
     &    p4,p6,p5)*cxw**(-2) - za(p2,p3)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)
     &    *zba2(p2,p4,p6,p5)*cxw**(-2) - za(p3,p7)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*zba2(p7,p4,p6,p5)*cxw**(-2) - 2.D0*za(p3,p7)*zb(p4,
     &    p6)*zab2(p8,p1,p7,p2)*zba2(p1,p4,p6,p5)*cxw**(-2) + za(p3,p8)
     &    *za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p8,p4,p6,p5)*cxw**(-2) +
     &    2.D0*za(p3,p8)*zb(p4,p6)*zab2(p7,p2,p8,p1)*zba2(p2,p4,p6,p5)*
     &    cxw**(-2) )
      amp = amp + gam3456f4*propw17**(-1)*propw28**(-1)*propw56**(-1)
     &  * (  - 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 + 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      amp = amp + gam3456f5*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s345**(-1) * ( za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p6)*zba2(p4,
     &    p3,p5,p1)*cxw**(-2) - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p2,p6)
     &    *zba2(p4,p3,p5,p2)*cxw**(-2) - za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p6,p7)*zba2(p4,p3,p5,p7)*cxw**(-2) + za(p3,p5)*za(p7,p8)*
     &    zb(p1,p2)*zb(p6,p8)*zba2(p4,p3,p5,p8)*cxw**(-2) + 2.D0*za(p3,
     &    p5)*zb(p1,p6)*zab2(p8,p1,p7,p2)*zba2(p4,p3,p5,p7)*cxw**(-2)
     &     - 2.D0*za(p3,p5)*zb(p2,p6)*zab2(p7,p2,p8,p1)*zba2(p4,p3,p5,
     &    p8)*cxw**(-2) )
      amp = amp + gam3456f5*propw17**(-1)*propw28**(-1)*propw34**(-1)
     &  * ( 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 - 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      amp = amp + gam3456f6*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s346**(-1) * (  - za(p1,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(
     &    p1,p4,p6,p3)*cxw**(-2) + za(p2,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p6)*zba2(p2,p4,p6,p3)*cxw**(-2) + za(p5,p7)*za(p7,p8)*zb(p1,
     &    p2)*zb(p4,p6)*zba2(p7,p4,p6,p3)*cxw**(-2) + 2.D0*za(p5,p7)*
     &    zb(p4,p6)*zab2(p8,p1,p7,p2)*zba2(p1,p4,p6,p3)*cxw**(-2) - za(
     &    p5,p8)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p8,p4,p6,p3)*
     &    cxw**(-2) - 2.D0*za(p5,p8)*zb(p4,p6)*zab2(p7,p2,p8,p1)*zba2(
     &    p2,p4,p6,p3)*cxw**(-2) )
      amp = amp + gam3456f6*propw17**(-1)*propw28**(-1)*propw34**(-1)
     &  * (  - 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 + 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )

      amp = amp*Bbit

!    original MCFM code
!       amp = amp + Hbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
!      & propw56**(-1) * (
!      &    -za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*proph3456**(-1)   ! MARKUS: this is the WW-->H-->WW contribution (s-channel Higgs)
!      &   - za(p3,p7)*za(p5,p8)*zb(p1,p4)*zb(p2,p6)*proph1347**(-1)   ! MARKUS: this is the WW-->H-->WW contribution (t-channel Higgs)
!      &                 )*cxw**(-3)*sqwmass

!     new code with anomalous couplings
      Amp_S_PR=czip
      Amp_S_DK=czip
      Amp_T_PR=czip
      Amp_T_DK=czip
      if( hmass.ge.zip ) then
        if(channeltoggle_stu.ne.1) then
          if( AnomalCouplPR.eq.1 ) then
      Amp_S_PR=anomhwwamp(i7,i1,i8,i2,1,s3456,s(i7,i1),s(i8,i2),za,zb)
          else
      Amp_S_PR=za(i7,i8)*zb(i1,i2)
          endif
          if( AnomalCouplDK.eq.1 ) then
      ! MCFM uses W-W+!
      Amp_S_DK=anomhwwamp(i5,i6,i3,i4,1,s3456,s(i5,i6),s(i3,i4),za,zb)
          else
      Amp_S_DK=za(i3,i5)*zb(i4,i6)
          endif
        endif
        if(channeltoggle_stu.ne.0) then
      ! MCFM uses W-W+!
          if( AnomalCouplPR.eq.1 ) then
      Amp_T_PR=-anomhwwamp(i7,i1,i3,i4,1,s1347,s(i7,i1),s(i3,i4),za,zb)
          else
      Amp_T_PR=za(i7,i3)*zb(i4,i1)
          endif
          if( AnomalCouplDK.eq.1 ) then
      ! MCFM uses W-W+!
      Amp_T_DK=-anomhwwamp(i5,i6,i8,i2,1,s1347,s(i5,i6),s(i8,i2),za,zb)
          else
      Amp_T_DK=za(i5,i8)*zb(i2,i6)
          endif
        endif
      amp = amp   !*00000
     & + propw17**(-1)*propw28**(-1)*propw34**(-1)*propw56**(-1)*(
     &   - Amp_S_PR*Amp_S_DK*proph3456**(-1)
     &   - Amp_T_PR*Amp_T_DK*proph1347**(-1)
     &                 )*Hbit*facWWHWW
      endif

!     adding a second resonance
      Amp_S_PR=czip
      Amp_S_DK=czip
      Amp_T_PR=czip
      Amp_T_DK=czip
      if( h2mass.ge.zip ) then
        if(channeltoggle_stu.ne.1) then
          if( AnomalCouplPR.eq.1 ) then
      Amp_S_PR=anomhwwamp(i7,i1,i8,i2,2,s3456,s(i7,i1),s(i8,i2),za,zb)
          else
      Amp_S_PR=za(i7,i8)*zb(i1,i2)
          endif
          if( AnomalCouplDK.eq.1 ) then
      ! MCFM uses W-W+!
      Amp_S_DK=anomhwwamp(i5,i6,i3,i4,2,s3456,s(i5,i6),s(i3,i4),za,zb)
          else
      Amp_S_DK=za(i3,i5)*zb(i4,i6)
          endif
        endif
        if(channeltoggle_stu.ne.0) then
      ! MCFM uses W-W+!
          if( AnomalCouplPR.eq.1 ) then
      Amp_T_PR=-anomhwwamp(i7,i1,i3,i4,2,s1347,s(i7,i1),s(i3,i4),za,zb)
          else
      Amp_T_PR=za(i7,i3)*zb(i4,i1)
          endif
          if( AnomalCouplDK.eq.1 ) then
      ! MCFM uses W-W+!
      Amp_T_DK=-anomhwwamp(i5,i6,i8,i2,2,s1347,s(i5,i6),s(i8,i2),za,zb)
          else
      Amp_T_DK=za(i5,i8)*zb(i2,i6)
          endif
        endif
      amp = amp   !*00000
     & + propw17**(-1)*propw28**(-1)*propw34**(-1)*propw56**(-1)*(
     &   - Amp_S_PR*Amp_S_DK*propX3456**(-1)
     &   - Amp_T_PR*Amp_T_DK*propX1347**(-1)
     &                 )*Hbit*facWWHWW
      endif

      return
      end
