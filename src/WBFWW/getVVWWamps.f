      subroutine getVVWWamps(amp,ampa,ampb,p,za,zb,zab,zba,
     & j1,j2,j3,j4,j5,j6,j7,j8,doHO,doBO)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      include 'WWbits.f'
      integer nmax,jmax
      parameter(nmax=10)
      integer k,l,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq,
     & nfinc
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2,j1,j2,j3,j4,j5,j6,j7,j8
      double precision p(mxpart,4)
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),cdotpr,
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2),
     & k7341(4),k8341(4),k8342(4),
     & s7341,s8341,s8342
      double complex
     & jmid17(2,2,2,2),jvbf17(2,2,2,2),jtwodiags17(2,2,2,2),
     &jtwo17(2,2,2,2),jtwo28(2,2,2,2),jZWZa17(2,2,2,2),jZWZb17(2,2,2,2),
     & jmid18(2,2,2,2),jvbf18(2,2,2,2),jtwodiags18(2,2,2,2),
     &jtwo18(2,2,2,2),jtwo27(2,2,2,2),jZWZa18(2,2,2,2),jZWZb18(2,2,2,2),
     & jtwoWexch17(2,2),jtwoWexch28(2,2),
     & j7_34_1z(2,4),j7_34_1g(2,4),j8_56_2z(2,4),j8_56_2g(2,4),
     & jtwoWexch18(2,2),jtwoWexch27(2,2),
     & j8_34_1z(2,4),j8_34_1g(2,4),j7_56_2z(2,4),j7_56_2g(2,4),
     & jmidWW17,jmidWW18,jmidWW28,
     & j8_34_2z(2,4),j8_34_2g(2,4),j7_56_1z(2,4),j7_56_1g(2,4)
      logical doHO,doBO

      amp(:,:,:)=czip
      ampa(:,:,:)=czip
      ampb(:,:,:)=czip

      k7341(:)=0.5d0*(zab(j1,:,j1)+zab(j3,:,j3)
     & +zab(j4,:,j4)+zab(j7,:,j7))
      k8341(:)=0.5d0*(zab(j1,:,j1)+zab(j3,:,j3)
     & +zab(j4,:,j4)+zab(j8,:,j8))
      k8342(:)=0.5d0*(zab(j2,:,j2)+zab(j3,:,j3)
     & +zab(j4,:,j4)+zab(j8,:,j8))
      s7341=cdotpr(k7341,k7341)
      s8341=cdotpr(k8341,k8341)
      s8342=cdotpr(k8342,k8342)

c--- These contributions contain Hbit and Bbit
c--- contribution from jVBF
      call ampvbf(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,jvbf17)
      call ampvbf(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,jvbf18)

c--- W mid diagrams: contribution from jmidWW
      call ampmidWW(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,jmidWW17)
      call ampmidWW(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,jmidWW18)
      call ampmidWW(j2,j1,j3,j4,j5,j6,j8,j7,za,zb,jmidWW28)

c--- these are not used in calculation of Higgs contribution
      if (doHO .eqv. .false.) then

c--- mid diagrams: contribution from jcentre
      call ampmid(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,jmid17)
      call ampmid(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,jmid18)
c--- contribution from jtwoWW
      call amp2current(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,jtwo17)
      call amp2current(j2,j1,j3,j4,j5,j6,j8,j7,za,zb,jtwo28)
      call amp2current(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,jtwo18)
      call amp2current(j2,j1,j3,j4,j5,j6,j7,j8,za,zb,jtwo27)

c--- Z/W/Z diagrams: contribution from jZWZ
      call ampZWZ(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,jZWZa17)
      call ampZWZ(j1,j2,j5,j6,j3,j4,j7,j8,za,zb,jZWZb17)
      call ampZWZ(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,jZWZa18)
      call ampZWZ(j1,j2,j5,j6,j3,j4,j8,j7,za,zb,jZWZb18)

c--- W-exchange diagrams: contribution from jtwodiags
      call amptwodiags(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,
     & jtwodiags17)
      call amptwodiags(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,
     & jtwodiags18)

c--- W exchange diagrams for flavor-changing contributions: contribution from jtwoWexch
      call amp2currentw(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,
     & jtwoWexch17)
      call amp2currentw(j2,j1,j3,j4,j5,j6,j8,j7,za,zb,
     & jtwoWexch28)
      call amp2currentw(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,
     & jtwoWexch18)
      call amp2currentw(j2,j1,j3,j4,j5,j6,j7,j8,za,zb,
     & jtwoWexch27)

      call jonew(j7,j3,j4,j1,za,zb,zab,j7_34_1z,j7_34_1g)
      call jonew(j8,j5,j6,j2,za,zb,zab,j8_56_2z,j8_56_2g)
      call jonew(j8,j3,j4,j1,za,zb,zab,j8_34_1z,j8_34_1g)
      call jonew(j7,j5,j6,j2,za,zb,zab,j7_56_2z,j7_56_2g)
      call jonew(j8,j3,j4,j2,za,zb,zab,j8_34_2z,j8_34_2g)
      call jonew(j7,j5,j6,j1,za,zb,zab,j7_56_1z,j7_56_1g)

      else

      jmid17=czip
      jmid18=czip
      jtwo17=czip
      jtwo28=czip
      jtwo18=czip
      jtwo27=czip
      jZWZa17=czip
      jZWZb17=czip
      jZWZa18=czip
      jZWZb18=czip
      jtwodiags17=czip
      jtwodiags18=czip
      jtwoWexch17=czip
      jtwoWexch28=czip
      jtwoWexch18=czip
      jtwoWexch27=czip
      j7_34_1z=czip
      j7_34_1g=czip
      j8_56_2z=czip
      j8_56_2g=czip
      j8_34_1z=czip
      j8_34_1g=czip
      j7_56_2z=czip
      j7_56_2g=czip
      j8_34_2z=czip
      j8_34_2g=czip
      j7_56_1z=czip
      j7_56_1g=czip

      endif


C-----setup for (dqcq_dqcq)
      do h1=1,2
      do h2=1,2
      amp(dqcq_dqcq,h1,h2)=
     & +jmid17(1,2,h1,h2)
     & +jvbf17(1,2,h1,h2)
     & +jtwo17(1,2,h1,h2)+jtwo28(2,1,h2,h1)
     & +jZWZa17(1,2,h1,h2)+jZWZb17(1,2,h1,h2)
     & +jtwodiags17(1,2,h1,h2)
      enddo
      enddo

C-----setup for (dqcq_uqsq)
      amp(dqcq_uqsq,1,1)=
     & +jtwoWexch17(1,2)+jtwoWexch28(2,1)
     & +jmidWW17
     & +cdotpr(j7_34_1g(1,:),j8_56_2g(2,:))/s7341
     & +(cdotpr(j7_34_1z(1,:),j8_56_2z(2,:))
     &  -cdotpr(j7_34_1z(1,:),k7341(:))
     &  *cdotpr(k7341(:),j8_56_2z(2,:))/czmass2)
     & /(s7341-dcmplx(zmass**2,-zmass*zwidth))

C-----setup for (uqcq_uqcq)
      do h1=1,2
      do h2=1,2
c--- contribution from jcentre
      amp(uqcq_uqcq,h1,h2)=
     & +jmid17(2,2,h1,h2)
     & +jvbf17(2,2,h1,h2)
     & +jtwo17(2,2,h1,h2)+jtwo28(2,2,h2,h1)
     & +jZWZa17(2,2,h1,h2)+jZWZb17(2,2,h1,h2)
     & +jtwodiags17(2,2,h1,h2)
      enddo
      enddo


C-----setup for (dqsq_dqsq)
      do h1=1,2
      do h2=1,2
c--- contribution from jcentre
      amp(dqsq_dqsq,h1,h2)=
     & +jmid17(1,1,h1,h2)
     & +jvbf17(1,1,h1,h2)
     & +jtwo17(1,1,h1,h2)+jtwo28(1,1,h2,h1)
     & +jZWZa17(1,1,h1,h2)+jZWZb17(1,1,h1,h2)
     & +jtwodiags17(1,1,h1,h2)
      enddo
      enddo

C-----setup for (dqdq_dqdq)
      do h1=1,2
      do h2=1,2
c--- contribution from jcentre
c-------- ampa
      ampa(dqdq_dqdq,h1,h2)=amp(dqsq_dqsq,h1,h2)
c-------- ampb
      ampb(dqdq_dqdq,h1,h2)=
     & +jmid18(1,1,h1,h2)
     & +jvbf18(1,1,h1,h2)
     & +jtwo18(1,1,h1,h2)+jtwo27(1,1,h2,h1)
     & +jZWZa18(1,1,h1,h2)+jZWZb18(1,1,h1,h2)
     & +jtwodiags18(1,1,h1,h2)
      enddo
      enddo

C-----setup for (uquq_uquq)
      do h1=1,2
      do h2=1,2
c--- contribution from jcentre
c-------- ampa
      ampa(uquq_uquq,h1,h2)=amp(uqcq_uqcq,h1,h2)
c-------- ampb
      ampb(uquq_uquq,h1,h2)=
     & +jmid18(2,2,h1,h2)
     & +jvbf18(2,2,h1,h2)
     & +jtwo18(2,2,h1,h2)+jtwo27(2,2,h2,h1)
     & +jZWZa18(2,2,h1,h2)+jZWZb18(2,2,h1,h2)
     & +jtwodiags18(2,2,h1,h2)
      enddo
      enddo

C-----setup for (dquq_dquq)
c-------- ampa
      do h1=1,2
      do h2=1,2
      ampa(dquq_dquq,h1,h2)=amp(dqcq_dqcq,h1,h2)
      enddo
      enddo
c-------- ampb
      ampb(dquq_dquq,1,1)=
     & +jtwoWexch18(1,2)+jtwoWexch27(2,1)
     & +jmidWW18
     & +cdotpr(j8_34_1g(1,:),j7_56_2g(2,:))/s8341
     & +(cdotpr(j8_34_1z(1,:),j7_56_2z(2,:))
     &  -cdotpr(j8_34_1z(1,:),k8341(:))
     &  *cdotpr(k8341(:),j7_56_2z(2,:))/czmass2)
     & /(s8341-dcmplx(zmass**2,-zmass*zwidth))

C-----setup for (uqbq_uqbq)
      do h1=1,2
      do h2=1,2
      amp(uqbq_uqbq,h1,h2)=
     & +jmid17(2,1,h1,h2)
     & +jvbf17(2,1,h1,h2)
     & +jtwo17(2,1,h1,h2)+jtwo28(1,2,h2,h1)
     & +jZWZa17(2,1,h1,h2)+jZWZb17(2,1,h1,h2)
     & +jtwodiags17(2,1,h1,h2)
      enddo
      enddo

C-----setup for (uqsq_dqcq)
      amp(uqsq_dqcq,1,1)=
     & +jtwoWexch17(2,1)+jtwoWexch28(1,2)
     & +jmidWW28
     & +cdotpr(j8_34_2g(1,:),j7_56_1g(2,:))/s8342
     & +(cdotpr(j8_34_2z(1,:),j7_56_1z(2,:))
     &  -cdotpr(j8_34_2z(1,:),k8342(:))
     &  *cdotpr(k8342(:),j7_56_1z(2,:))/czmass2)
     & /(s8342-dcmplx(zmass**2,-zmass*zwidth))


      return

      end

