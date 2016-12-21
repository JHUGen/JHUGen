      subroutine getQQZZQQstrongamps(amp,ampa,ampb,p,za,zb,zab,zba,
     & j1,j2,j3,j4,j5,j6,j7,j8)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'WWbits.f'
      integer nmax
      parameter(nmax=10)
      integer k,l,i1,i2,i3,i4,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq
c     & dquq_dquq,dqcq_uqsq,uqsq_dqcq
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6)
c     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2,h3,h5,j1,j2,j3,j4,j5,j6,j7,j8
      double precision p(mxpart,4),
     & t4,s17,s28,s18,s27,s7341,s7561,s7342,s7562
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),cdotpr,
     & j7_1(4,2),j7_2(4,2),j8_1(4,2),j8_2(4,2),
     & j7_34_1(4,2,2,2),j7_34_2(4,2,2,2),
     & j7_56_1(4,2,2,2),j7_56_2(4,2,2,2),
     & j8_34_1(4,2,2,2),j8_34_2(4,2,2,2),
     & j8_56_1(4,2,2,2),j8_56_2(4,2,2,2),
     & jl7_34_1(4,2,2,2),jl7_34_2(4,2,2,2),
     & jl7_56_1(4,2,2,2),jl7_56_2(4,2,2,2),
     & jl8_34_1(4,2,2,2),jl8_34_2(4,2,2,2),
     & jl8_56_1(4,2,2,2),jl8_56_2(4,2,2,2),
     & jw7_34_1(4,2,2),jw7_34_2(4,2,2),
     & jw7_56_1(4,2,2),jw7_56_2(4,2,2),
     & jw8_34_1(4,2,2),jw8_34_2(4,2,2),
     & jw8_56_1(4,2,2),jw8_56_2(4,2,2),
     & j7_3456_1(4,2,2,2,2),j8_3456_2(4,2,2,2,2),
     & j7_3456_2(4,2,2,2,2),j8_3456_1(4,2,2,2,2),
     & jw7_3456_1(4,2,2,2),jw8_3456_2(4,2,2,2),
     & jw7_3456_2(4,2,2,2),jw8_3456_1(4,2,2,2),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2)
c--- Begin statement functions
      t4(i1,i2,i3,i4)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)
     & +s(i2,i3)+s(i2,i4)+s(i3,i4)
c--- End statement functions

      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

      s17=s(j1,j7)
      s28=s(j2,j8)
      s27=s(j2,j7)
      s18=s(j1,j8)
      s7341=t4(j7,j3,j4,j1)
      s7342=t4(j7,j3,j4,j2)
      s7561=t4(j7,j5,j6,j1)
      s7562=t4(j7,j5,j6,j2)

      call jzero(j7,j1,zab,zba,j7_1)
      call jzero(j7,j2,zab,zba,j7_2)
      call jzero(j8,j1,zab,zba,j8_1)
      call jzero(j8,j2,zab,zba,j8_2)

      call jone(j7,j3,j4,j1,za,zb,zab,zba,j7_34_1,jw7_34_1,jl7_34_1)
      call jone(j7,j3,j4,j2,za,zb,zab,zba,j7_34_2,jw7_34_2,jl7_34_2)
      call jone(j7,j5,j6,j1,za,zb,zab,zba,j7_56_1,jw7_56_1,jl7_56_1)
      call jone(j7,j5,j6,j2,za,zb,zab,zba,j7_56_2,jw7_56_2,jl7_56_2)
      call jone(j8,j3,j4,j1,za,zb,zab,zba,j8_34_1,jw8_34_1,jl8_34_1)
      call jone(j8,j3,j4,j2,za,zb,zab,zba,j8_34_2,jw8_34_2,jl8_34_2)
      call jone(j8,j5,j6,j1,za,zb,zab,zba,j8_56_1,jw8_56_1,jl8_56_1)
      call jone(j8,j5,j6,j2,za,zb,zab,zba,j8_56_2,jw8_56_2,jl8_56_2)

      call jtwo(j7,j3,j4,j5,j6,j1,za,zb,zab,zba,j7_3456_1,jw7_3456_1)
      call jtwo(j7,j3,j4,j5,j6,j2,za,zb,zab,zba,j7_3456_2,jw7_3456_2)
      call jtwo(j8,j3,j4,j5,j6,j1,za,zb,zab,zba,j8_3456_1,jw8_3456_1)
      call jtwo(j8,j3,j4,j5,j6,j2,za,zb,zab,zba,j8_3456_2,jw8_3456_2)

C-----setup for (uqbq_uqbq) (2,5)->(2,5)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C---one-one currents
      amp(uqbq_uqbq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,1,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,1,h2,h3))/s7561

C---two-one currents
      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,1,h2,h3,h5))/s17
      enddo
      enddo
      enddo
      enddo

C-----setup for (uqcq_uqcq) (2,4)->(2,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C---one-one currents
      amp(uqcq_uqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,2,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,2,h2,h3))/s7561

C---two-one currents
      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))/s17
      enddo
      enddo
      enddo
      enddo


C-----setup for (dqcq_dqcq) (1,4)-->(1,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      amp(dqcq_dqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,2,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,2,h2,h3))/s7561

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))/s17
      enddo
      enddo
      enddo
      enddo

C-----setup for (dqsq_dqsq) (1,3)-->(1,3)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      amp(dqsq_dqsq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,1,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,1,h2,h3))/s7561

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j8_3456_2(:,1,h2,h3,h5),j7_1(:,h1))/s17
      enddo
      enddo
      enddo
      enddo

C-----setup for ((uquq_uquq)  (2,2)-->(2,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C-----------------ampa
      ampa(uquq_uquq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)

C-----------------ampb
      ampb(uquq_uquq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,2,h1,h3),j7_56_2(:,2,h2,h5))/s7562
     & +cdotpr(j8_56_1(:,2,h1,h5),j7_34_2(:,2,h2,h3))/s7342

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,2,h1,h3,h5),j7_2(:,h2))/s27
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,2,h2,h3,h5))/s18
      enddo
      enddo
      enddo
      enddo

C-----setup for ((dqdq_dqdq)  (1,1)-->(1,1)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C-----------------ampa
      ampa(dqdq_dqdq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)

C-----------------ampb
      ampb(dqdq_dqdq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,1,h1,h3),j7_56_2(:,1,h2,h5))/s7562
     & +cdotpr(j8_56_1(:,1,h1,h5),j7_34_2(:,1,h2,h3))/s7342

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,1,h1,h3,h5),j7_2(:,h2))/s27
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,1,h2,h3,h5))/s18
      enddo
      enddo
      enddo
      enddo

      ! U.Sarica: interference terms receive extra (-) in amp**2
      ampb(:,:,:,:,:)=-ampb(:,:,:,:,:)

      return

      end
