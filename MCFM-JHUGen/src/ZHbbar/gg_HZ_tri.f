!----- fills amplitudes for gg=>HZ process
!===== these are the pieces which come from triangle topologies
      subroutine gg_HZ_tri(p,amp,mt2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'cplx.h'
      real(kind=dp):: p(mxpart,4),mt2
      complex(kind=dp):: amp(2,2,2),test
      complex(kind=dp):: ggHZ_pp_tri
      external ggHZ_pp_tri
      complex(kind=dp):: prop_34,prop_12

!----- debug
!----- KC phase space point check
!      p(1,4)=  -3.0000000000000000_dp
!      p(1,1)=   2.1213203435596424_dp
!      p(1,2)=   1.0606601717798212_dp
!      p(1,3)=   1.8371173070873839_dp
!      p(2,4)=  -3.0000000000000000_dp
!      p(2,1)=  -1.3016512173526746_dp
!      p(2,2)=   2.6346884688170253_dp
!      p(2,3)= -0.60342421284441206_dp
!      p(3,4)=  0.85714285714285710_dp
!      p(3,1)= -0.31578947368421051_dp
!      p(3,2)=  0.79685060448070799_dp
!      p(3,3)=   0.0000000000000000_dp
!     p(4,4)=  -1.50000000000000000_dp
!      p(4,1)= -0.88167787843870971_dp
!      p(4,2)=   1.1895840590472970_dp
!      p(4,3)= -0.23986222114450001_dp


!      scale=one
!      musq=one
!      mt=two
!      mt=2.32one
!      mt=0.4255266775_dp

      call spinoru(4,p,za,zb)
!====== z propagtors
      prop_34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      prop_12=s(1,2)/cplx2(s(1,2)-zmass**2,zmass*zwidth)
!---- for checking total xs use the formula below
!      prop_12=s(1,2)/cplx2(s(1,2)-zmass**2,zip)

 !------ left lepton amplitudes
      amp(2,2,1)=ggHZ_pp_tri(1,2,3,4,za,zb,mt2)
      amp(2,1,1)=czip
      amp(1,2,1)=czip
      amp(1,1,1)=-ggHZ_pp_tri(1,2,4,3,zb,za,mt2)

!------ right lepton amplitudes
      amp(2,2,2)=ggHZ_pp_tri(1,2,4,3,za,zb,mt2)
      amp(2,1,2)=czip
      amp(1,2,2)=czip
      amp(1,1,2)=-ggHZ_pp_tri(1,2,3,4,zb,za,mt2)


!---- debug KC check
 !     write(6,*) '************ L *********'
 !     write(6,*) 'amp(2,2,1)*im ', amp(2,2,1)*im
 !     write(6,*) 'amp(1,2,1)*im ', amp(1,2,1)*im
 !     write(6,*) 'amp(2,1,1)*im ', amp(2,1,1)*im
 !     write(6,*) 'amp(1,1,1)*im ', amp(1,1,1)*im
 !     write(6,*) '************************'
 !     write(6,*) '************ R *********'
 !     write(6,*) 'amp(2,2,2)*im ', amp(2,2,2)*im
 !     write(6,*) 'amp(1,2,2)*im ', amp(1,2,2)*im
 !     write(6,*) 'amp(2,1,2)*im ', amp(2,1,2)*im
 !     write(6,*) 'amp(1,1,2) *im', amp(1,1,2)*im
!      write(6,*) '************************'
!      stop

      amp(:,:,:)=prop_12*prop_34*amp(:,:,:)
      return
      end


      function ggHZ_pp_tri(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(kind=dp)::ggHZ_pp_tri
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'masses.f'
      integer:: i1,i2,i3,i4
      real(kind=dp):: mt2
      real(kind=dp):: mZsq,mHsq,s15,s25,s12,t
      integer:: Nbox,Ntri,i
      parameter(Nbox=3,Ntri=6)

      complex(kind=dp):: qlI4,qlI3,qlI2
      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5)
      complex(kind=dp):: D0(Nbox),C0(Ntri)
      complex(kind=dp):: di(Nbox),ci(Ntri)
      complex(kind=dp):: rat
      complex(kind=dp):: cmzsq,cmz
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)

      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_pp_tri=czip

!--- debug
!      zwidth=zip
!      zmass=2.5_dp
!      zwidth=zip

!----- note I'm using the complex mass scheme to ensure Gauge invariance,
      cmzsq=zmass**2-im*zmass*zwidth


!======= kinematic configurations
      s12=s(i1,i2)
      mZsq=s(i3,i4)
      s25=t(i1,i3,i4)
      s15=t(i2,i3,i4)
      mHsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      di(:)=czip
      ci(:)=czip



!      C0(c12)=qlI3(s12,zip,zip,mt2,mt2,mt2,musq,0)

      ci(c12)=  (2*(-((mt2*za(i1,i3)*zb(i2,i1)*zb(i4,i1))/za(i1,i2)) -
     -      (mt2*za(i2,i3)*zb(i2,i1)*zb(i4,i2))/za(i1,i2)))/
     -  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))

!------ longitudinal pieces
      ci(c12)=ci(c12)+
     &(2*mt2*zb(i2,i1)**2*(za(i1,i3)*zb(i4,i1) + za(i2,i3)*zb(i4,i2)))/
     -  (cmzsq*s(i1,i2)*s(i3,i4))

      rat= (-2*(-((za(i1,i3)*zb(i2,i1)*zb(i4,i1))/za(i1,i2)) -
     -      (za(i2,i3)*zb(i2,i1)*zb(i4,i2))/za(i1,i2)))/
     -  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))

!------ longitudinal pieces

      rat=rat+
     &  (-2*zb(i2,i1)**2*
     -    (za(i1,i3)*zb(i4,i1) + za(i2,i3)*zb(i4,i2)))/
     -  (cmzsq*s(i1,i2)*s(i3,i4))


!----- debug
!      write(6,*) 'coefficient of tri*im ',ci(c12)*im
!     write(6,*) 'rational term*im ',rat*im

      ggHZ_pp_tri=C0(c12)*ci(c12)-rat/two

 !     write(6,*) dsqrt(mt2),C0(c12),ci(c12),ggHZ_pp_tri
 !     pause
      return
      end
