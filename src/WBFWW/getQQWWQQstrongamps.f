      subroutine getQQWWQQstrongamps(amp,ampa,ampb,p,za,zb,zab,zba,
     & j1,j2,j3,j4,j5,j6,j7,j8)
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
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq
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
     & jtwo17(2,2,2,2),jtwo28(2,2,2,2),
     & jtwo18(2,2,2,2),jtwo27(2,2,2,2),
     & j7_34_1g(2,4),j8_56_2g(2,4),
     & j8_34_1g(2,4),j7_56_2g(2,4),
     & j8_34_2g(2,4),j7_56_1g(2,4)

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

c--- contribution from jtwoWW
      call amp2currentstrong(j1,j2,j3,j4,j5,j6,j7,j8,za,zb,
     & jtwo17)
      call amp2currentstrong(j2,j1,j3,j4,j5,j6,j8,j7,za,zb,
     &jtwo28)
      call amp2currentstrong(j1,j2,j3,j4,j5,j6,j8,j7,za,zb,
     & jtwo18)
      call amp2currentstrong(j2,j1,j3,j4,j5,j6,j7,j8,za,zb,
     & jtwo27)

c--- flavor-changing contributions
      call jonewstrong(j7,j3,j4,j1,za,zb,zab,j7_34_1g)
      call jonewstrong(j8,j5,j6,j2,za,zb,zab,j8_56_2g)
      call jonewstrong(j8,j3,j4,j1,za,zb,zab,j8_34_1g)
      call jonewstrong(j7,j5,j6,j2,za,zb,zab,j7_56_2g)
      call jonewstrong(j8,j3,j4,j2,za,zb,zab,j8_34_2g)
      call jonewstrong(j7,j5,j6,j1,za,zb,zab,j7_56_1g)


C-----setup for (dqcq_dqcq)
      do h1=1,2
      do h2=1,2
      amp(dqcq_dqcq,h1,h2)=
     & +jtwo17(1,2,h1,h2)+jtwo28(2,1,h2,h1)
      enddo
      enddo

C-----setup for (dqcq_uqsq)
      amp(dqcq_uqsq,1,1)=
     & +cdotpr(j7_34_1g(1,:),j8_56_2g(2,:))/s7341

C-----setup for (uqcq_uqcq)
      do h1=1,2
      do h2=1,2
c--- contribution from jcentre
      amp(uqcq_uqcq,h1,h2)=
     & +jtwo17(2,2,h1,h2)+jtwo28(2,2,h2,h1)
      enddo
      enddo

C-----setup for (dqsq_dqsq)
      do h1=1,2
      do h2=1,2
c--- contribution from jcentre
      amp(dqsq_dqsq,h1,h2)=
     & +jtwo17(1,1,h1,h2)+jtwo28(1,1,h2,h1)
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
     & +jtwo18(1,1,h1,h2)+jtwo27(1,1,h2,h1)
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
     & +jtwo18(2,2,h1,h2)+jtwo27(2,2,h2,h1)
      enddo
      enddo

C-----setup for (dquq_dquq)
c-------- ampb
      ampb(dquq_dquq,1,1)=
     & +cdotpr(j8_34_1g(1,:),j7_56_2g(2,:))/s8341

      do h1=1,2
      do h2=1,2
c-------- ampa
      ampa(dquq_dquq,h1,h2)=amp(dqcq_dqcq,h1,h2)
      enddo
      enddo

C-----setup for (uqbq_uqbq)
      do h1=1,2
      do h2=1,2
      amp(uqbq_uqbq,h1,h2)=
     & +jtwo17(2,1,h1,h2)+jtwo28(1,2,h2,h1)
      enddo
      enddo

C-----setup for (uqsq_dqcq)
      amp(uqsq_dqcq,1,1)=
     & +cdotpr(j8_34_2g(1,:),j7_56_1g(2,:))/s8342


      return

      end

