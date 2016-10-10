      subroutine ampmidWW(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cmplxmass.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zcouple.f'
      include 'WWbits.f'
      include 'spinzerohiggs_anomcoupl.f'            
      integer i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      double complex zab2,zba2,amp,game,gamn,sqwmass,rxw,
     & propw34,propw56,propw28,propw17,anomhwwamp,
     & Amp_S_PR,Amp_S_DK,Amp_T_PR,Amp_T_DK,propX3456,
     & propz3456,propa3456,proph3456,proph1347,propX1347
      double precision t3,t4,s34,s56,s17,s28,s137,s147,
     & s258,s268,s456,s345,s356,s346,s3456,s1347,s1567,
     & twop17Dp3456,twop28Dp3456,twop34Dp3456,twop56Dp3456
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

      rxw=sqrt((cone-cxw)/cxw)
      s3456=t4(i3,i4,i5,i6)
      s1347=t4(i1,i3,i4,i7)
      s1567=t4(i1,i5,i6,i7)
      s137=t3(i1,i3,i7)
      s147=t3(i1,i4,i7)
      s258=t3(i2,i5,i8)
      s268=t3(i2,i6,i8)
      s345=t3(i3,i4,i5)
      s346=t3(i3,i4,i6)
      s356=t3(i3,i5,i6)
      s456=t3(i4,i5,i6)
      s34=s(i3,i4)
      s56=s(i5,i6)
      s17=s(i1,i7)
      s28=s(i2,i8)

      twop17Dp3456=s28-s17-s3456
      twop28Dp3456=s17-s28-s3456
      twop34Dp3456=s34-s56+s3456
      twop56Dp3456=s56-s34+s3456

      proph3456=dcmplx(s3456-hmass**2,hmass*hwidth)
      proph1347=dcmplx(s1347-hmass**2,hmass*hwidth)
      propz3456=s3456-dcmplx(zmass**2,-zmass*zwidth)
      propa3456=s3456
      propw34=s34-dcmplx(wmass**2,-wmass*wwidth)
      propw56=s56-dcmplx(wmass**2,-wmass*wwidth)
      propw17=s17-dcmplx(wmass**2,-wmass*wwidth)
      propw28=s28-dcmplx(wmass**2,-wmass*wwidth)

      game=qe/propa3456+le*rxw/propz3456
      gamn=ln*rxw/propz3456

      p1=i1
      p2=i2
      p3=i3
      p4=i4
      p5=i5
      p6=i6
      p7=i7
      p8=i8
      amp= + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
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
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
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
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
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
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propz3456**(-1) * (  - za(p7,p8)*zb(p1,p2)*zab2(p3
     &    ,p5,p6,p4)*zab2(p5,p1,p7,p6)*cxw**(-2) - za(p7,p8)*zb(p1,p2)*
     &    zab2(p3,p5,p6,p4)*zab2(p5,p2,p8,p6)*cxw**(-3) + za(p7,p8)*zb(
     &    p1,p2)*zab2(p3,p5,p6,p4)*zab2(p5,p2,p8,p6)*cxw**(-2) )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
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
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1)*propa3456**(-1) * ( za(p7,p8)*zb(p1,p2)*zab2(p3,p5
     &    ,p6,p4)*zab2(p5,p1,p7,p6)*cxw**(-2) - za(p7,p8)*zb(p1,p2)*
     &    zab2(p3,p5,p6,p4)*zab2(p5,p2,p8,p6)*cxw**(-2) )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propw56**(-1) * (  - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    cxw**(-3) - za(p3,p7)*za(p5,p8)*zb(p1,p4)*zb(p2,p6)*cxw**(-3)
     &     + 2.D0*za(p3,p8)*za(p5,p7)*zb(p1,p6)*zb(p2,p4)*cxw**(-3) )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & propa3456**(-1) * (  - 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*czmass2**(-1)*cxw**(-2)*twop28Dp3456*qe + 1.D0/2.D0
     &    *za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*
     &    cxw**(-2)*twop17Dp3456*qe )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s346**(-1)*s258**(-1) * (  - za(p5,p8)*zb(p4,p6)*zba2(p1,p4,p6,
     &    p3)*zba2(p2,p5,p8,p7)*cxw**(-3) )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s345**(-1)*s268**(-1) * ( za(p3,p5)*zb(p2,p6)*zba2(p1,p2,p6,p8)*
     &    zba2(p4,p3,p5,p7)*cxw**(-3) )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & propa3456**(-1) * ( 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(
     &    p4,p6)*czmass2**(-1)*cxw**(-2)*twop28Dp3456*qe - 1.D0/2.D0*
     &    za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*
     &    cxw**(-2)*twop17Dp3456*qe )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s456**(-1) * ( za(p3,p7)*zb(p4,p6)*zba2(p1,p3,p7,p8)*zba2(p2,p4,
     &    p6,p5)*cxw**(-3)*s137**(-1) )
      amp = amp + Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s356**(-1)*s147**(-1) * (  - za(p3,p5)*zb(p1,p4)*zba2(p2,p1,p4,
     &    p7)*zba2(p6,p3,p5,p8)*cxw**(-3) )
      
!    original MCFM code
!       amp = amp + Hbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
!      & propw56**(-1) * (  
!      &    -za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*proph3456**(-1)   ! MARKUS: this is the WW-->H-->WW contribution (s-channel Higgs)
!      &   - za(p3,p7)*za(p5,p8)*zb(p1,p4)*zb(p2,p6)*proph1347**(-1)   ! MARKUS: this is the WW-->H-->WW contribution (t-channel Higgs)
!      &                 )*cxw**(-3)*sqwmass
!       print *, "MARKUS check: set Bbit=0!"
!       print *, "MARKUS check: old WW-->H-->WW:",amp


!     new code with anomalous couplings
      Amp_S_PR=-anomhwwamp(i7,i1,i8,i2,1,s3456,s(i7,i1),s(i8,i2),za,zb)
      Amp_S_DK=-anomhwwamp(i3,i4,i5,i6,1,s3456,s(i3,i4),s(i5,i6),za,zb)
      Amp_T_PR=-anomhwwamp(i7,i1,i3,i4,1,s1347,s(i3,i4),s(i7,i1),za,zb)
      Amp_T_DK=-anomhwwamp(i8,i2,i5,i6,1,s1347,s(i5,i6),s(i8,i2),za,zb)
      amp = amp   !*00000 
     & + Hbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*propw56**(-1)*(  
     &   - Amp_S_PR*Amp_S_DK*proph3456**(-1)   ! MARKUS: this is the WW-->H-->WW contribution (s-channel Higgs)
     &   - Amp_T_PR*Amp_T_DK*proph1347**(-1)   ! MARKUS: this is the WW-->H-->WW contribution (t-channel Higgs)
     &                 )*cxw**(-3)*sqwmass
!       print *, "MARKUS check: remove *00000 above when checked"
!       print *, "MARKUS check: new WW-->H-->WW:",amp
!       pause
     
!     adding a second resonance
      if( h2mass.ge.zip ) then      
      propX3456=dcmplx(s3456-h2mass**2,h2mass*h2width)
      propX1347=dcmplx(s1347-h2mass**2,h2mass*h2width)
      Amp_S_PR=-anomhwwamp(i7,i1,i8,i2,2,s3456,s(i7,i1),s(i8,i2),za,zb)
      Amp_S_DK=-anomhwwamp(i3,i4,i5,i6,2,s3456,s(i3,i4),s(i5,i6),za,zb)
      Amp_T_PR=-anomhwwamp(i7,i1,i3,i4,2,s1347,s(i3,i4),s(i7,i1),za,zb)
      Amp_T_DK=-anomhwwamp(i8,i2,i5,i6,2,s1347,s(i5,i6),s(i8,i2),za,zb)
      amp = amp   !*00000 
     & + Hbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*propw56**(-1)*(  
     &   - Amp_S_PR*Amp_S_DK*propX3456**(-1)
     &   - Amp_T_PR*Amp_T_DK*propX1347**(-1)
     &                 )*cxw**(-3)*sqwmass
      endif
     
     
      amp = amp + gamn*Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s346**(-1) * (  - za(p1,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(
     &    p1,p4,p6,p3)*cxw**(-2) + za(p2,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,
     &    p6)*zba2(p2,p4,p6,p3)*cxw**(-2) + za(p5,p7)*za(p7,p8)*zb(p1,
     &    p2)*zb(p4,p6)*zba2(p7,p4,p6,p3)*cxw**(-2) + 2.D0*za(p5,p7)*
     &    zb(p4,p6)*zab2(p8,p1,p7,p2)*zba2(p1,p4,p6,p3)*cxw**(-2) - za(
     &    p5,p8)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p8,p4,p6,p3)*
     &    cxw**(-2) - 2.D0*za(p5,p8)*zb(p4,p6)*zab2(p7,p2,p8,p1)*zba2(
     &    p2,p4,p6,p3)*cxw**(-2) )
      amp = amp + gamn*Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)
     &  * (  - 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 + 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      amp = amp + gamn*Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s356**(-1) * (  - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*zba2(
     &    p6,p3,p5,p1)*cxw**(-2) + za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p2,
     &    p4)*zba2(p6,p3,p5,p2)*cxw**(-2) + za(p3,p5)*za(p7,p8)*zb(p1,
     &    p2)*zb(p4,p7)*zba2(p6,p3,p5,p7)*cxw**(-2) - za(p3,p5)*za(p7,
     &    p8)*zb(p1,p2)*zb(p4,p8)*zba2(p6,p3,p5,p8)*cxw**(-2) - 2.D0*
     &    za(p3,p5)*zb(p1,p4)*zab2(p8,p1,p7,p2)*zba2(p6,p3,p5,p7)*
     &    cxw**(-2) + 2.D0*za(p3,p5)*zb(p2,p4)*zab2(p7,p2,p8,p1)*zba2(
     &    p6,p3,p5,p8)*cxw**(-2) )
      amp = amp + gamn*Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)
     &  * ( 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 - 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      amp = amp + game*Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)*
     & s345**(-1) * ( za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p6)*zba2(p4,
     &    p3,p5,p1)*cxw**(-2) - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p2,p6)
     &    *zba2(p4,p3,p5,p2)*cxw**(-2) - za(p3,p5)*za(p7,p8)*zb(p1,p2)*
     &    zb(p6,p7)*zba2(p4,p3,p5,p7)*cxw**(-2) + za(p3,p5)*za(p7,p8)*
     &    zb(p1,p2)*zb(p6,p8)*zba2(p4,p3,p5,p8)*cxw**(-2) + 2.D0*za(p3,
     &    p5)*zb(p1,p6)*zab2(p8,p1,p7,p2)*zba2(p4,p3,p5,p7)*cxw**(-2)
     &     - 2.D0*za(p3,p5)*zb(p2,p6)*zab2(p7,p2,p8,p1)*zba2(p4,p3,p5,
     &    p8)*cxw**(-2) )
      amp = amp + game*Bbit*propw17**(-1)*propw28**(-1)*propw34**(-1)
     &  * ( 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 - 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      amp = amp + game*Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)*
     & s456**(-1) * ( za(p1,p3)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p1,
     &    p4,p6,p5)*cxw**(-2) - za(p2,p3)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)
     &    *zba2(p2,p4,p6,p5)*cxw**(-2) - za(p3,p7)*za(p7,p8)*zb(p1,p2)*
     &    zb(p4,p6)*zba2(p7,p4,p6,p5)*cxw**(-2) - 2.D0*za(p3,p7)*zb(p4,
     &    p6)*zab2(p8,p1,p7,p2)*zba2(p1,p4,p6,p5)*cxw**(-2) + za(p3,p8)
     &    *za(p7,p8)*zb(p1,p2)*zb(p4,p6)*zba2(p8,p4,p6,p5)*cxw**(-2) + 
     &    2.D0*za(p3,p8)*zb(p4,p6)*zab2(p7,p2,p8,p1)*zba2(p2,p4,p6,p5)*
     &    cxw**(-2) )
      amp = amp + game*Bbit*propw17**(-1)*propw28**(-1)*propw56**(-1)
     &  * (  - 1.D0/2.D0*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p4,p6)*
     &    czmass2**(-1)*cxw**(-2)*twop28Dp3456 + 1.D0/2.D0*za(p3,p5)*
     &    za(p7,p8)*zb(p1,p2)*zb(p4,p6)*czmass2**(-1)*cxw**(-2)*
     &    twop17Dp3456 )
      return
      end
