      function WZccmsq(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      real(dp):: WZccmsq

C-----Matrix elements squared for
C     u(p1)+c(p8)-->d(p2)+c(p7)+W(p3,p4)+Z(p5,p6) (nwz=+1)
C     d(p1)+c(p8)-->u(p2)+c(p7)+W(p3,p4)+Z(p5,p6) (nwz=-1)
C     Author:R.K.Ellis, February 2013
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'nwz.f'
      complex(dp):: TWZbbab,TWZbbnab,TWbbpZab,TWbbmZab,d123(2,2),
     & d456(2,2),d78(2,2),d910(2,2),prop34,prop56,prop3456,amp(2,2),
     & d11(2,2),d12(2,2),d13(2,2),TWZbbnr1,TWZbbnr2
      real(dp):: s3456,s4,fac,q3,q4,l3,l4
      integer:: jd,ju,p1,p2,p3,p4,p5,p6,p7,p8
C     statement functions
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
C     end statement functions

C     First index of d123 etc is helicity of Z (56) decay line
C     Second index of d123 etc is helicity of (78) bbar-line
      prop56=s(p5,p6)/cplx2(s(p5,p6)-zmass**2,zmass*zwidth)
      prop34=s(p3,p4)/cplx2(s(p3,p4)-wmass**2,wmass*wwidth)
      s3456=s4(p3,p4,p5,p6)
      prop3456=s3456/cplx2(s3456-wmass**2,wmass*wwidth)
      fac=aveqq*4._dp*V*gwsq**2*esq**2*gsq**2

      d123(1,1)=TWZbbab(p1,p2,p3,p4,p5,p6,p7,p8)
      d123(2,1)=TWZbbab(p1,p2,p3,p4,p6,p5,p7,p8)
      d123(1,2)=TWZbbab(p1,p2,p3,p4,p5,p6,p8,p7)
      d123(2,2)=TWZbbab(p1,p2,p3,p4,p6,p5,p8,p7)

      d456(1,1)=TWZbbab(p1,p2,p5,p6,p3,p4,p7,p8)
      d456(2,1)=TWZbbab(p1,p2,p6,p5,p3,p4,p7,p8)
      d456(1,2)=TWZbbab(p1,p2,p5,p6,p3,p4,p8,p7)
      d456(2,2)=TWZbbab(p1,p2,p6,p5,p3,p4,p8,p7)

      d78(1,1)=TWZbbnab(p1,p2,p3,p4,p5,p6,p7,p8)
      d78(2,1)=TWZbbnab(p1,p2,p3,p4,p6,p5,p7,p8)
      d78(1,2)=TWZbbnab(p1,p2,p3,p4,p5,p6,p8,p7)
      d78(2,2)=TWZbbnab(p1,p2,p3,p4,p6,p5,p8,p7)

      d910(1,1)=TWbbmZab(p1,p2,p3,p4,p5,p6,p7,p8)
      d910(2,1)=TWbbmZab(p1,p2,p3,p4,p6,p5,p7,p8)
      d910(1,2)=TWbbpZab(p1,p2,p3,p4,p5,p6,p7,p8)
      d910(2,2)=TWbbpZab(p1,p2,p3,p4,p6,p5,p7,p8)

      d11(1,1)=TWZbbnr1(p1,p2,p3,p4,p5,p6,p7,p8)
      d11(2,1)=TWZbbnr1(p1,p2,p3,p4,p6,p5,p7,p8)
      d11(1,2)=TWZbbnr1(p1,p2,p3,p4,p5,p6,p8,p7)
      d11(2,2)=TWZbbnr1(p1,p2,p3,p4,p6,p5,p8,p7)

      d12(1,1)=TWZbbnr2(p1,p2,p3,p4,p5,p6,p7,p8)
      d12(2,1)=TWZbbnr2(p1,p2,p3,p4,p6,p5,p7,p8)
      d12(1,2)=TWZbbnr2(p1,p2,p3,p4,p5,p6,p8,p7)
      d12(2,2)=TWZbbnr2(p1,p2,p3,p4,p6,p5,p8,p7)


C     For u->dW^+ process Z is emitted after W for 123 and before W for 456
      if (nwz == 1) then
      ju=2
      jd=1
      q3=0._dp
      l3=ln
      q4=qe
      l4=le
      d13(1,1)=TWZbbnr2(p1,p2,p5,p6,p3,p4,p7,p8)
      d13(1,2)=TWZbbnr2(p1,p2,p5,p6,p3,p4,p8,p7)
      d13(2,1)=0._dp
      d13(2,2)=0._dp
      elseif (nwz == -1) then
      jd=2
      ju=1
      q3=qe
      l3=le
      q4=0._dp
      l4=ln
      d13(1,1)=TWZbbnr1(p1,p2,p5,p6,p3,p4,p7,p8)
      d13(1,2)=TWZbbnr1(p1,p2,p5,p6,p3,p4,p8,p7)
      d13(2,1)=0._dp
      d13(2,2)=0._dp
      endif

!     First index of d123 etc is helicity of Z (56) decay line
!     Second index of d123 etc is helicity of (78) bbar-line
      amp(1,1)=
     &  ((Q(jd)*qe+L(jd)*le*prop56)*d123(1,1)
     &  +(Q(ju)*qe+L(ju)*le*prop56)*d456(1,1)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*le*prop56)*d78(1,1)*prop3456
     &  +(Q(2)*qe+L(2)*le*prop56)*d910(1,1)
     &  +0.5_dp/xw*d13(1,1)*prop3456)*prop34
     &  +(q4*qe+l4*le*prop56)*d11(1,1)*prop3456
     &  +(q3*qe+l3*le*prop56)*d12(1,1)*prop3456
     &
      amp(2,1)=
     &  ((Q(jd)*qe+L(jd)*re*prop56)*d123(2,1)
     &  +(Q(ju)*qe+L(ju)*re*prop56)*d456(2,1)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*re*prop56)*d78(2,1)*prop3456
     &  +(Q(2)*qe+L(2)*re*prop56)*d910(2,1)
     &  +0.5_dp/xw*d13(2,1)*prop3456)*prop34
     &  +(q4*qe+l4*re*prop56)*d11(2,1)*prop3456
     &  +(q3*qe+l3*re*prop56)*d12(2,1)*prop3456

      amp(1,2)=
     &  ((Q(jd)*qe+L(jd)*le*prop56)*d123(1,2)
     &  +(Q(ju)*qe+L(ju)*le*prop56)*d456(1,2)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*le*prop56)*d78(1,2)*prop3456
     &  +(Q(2)*qe+R(2)*le*prop56)*d910(1,2)
     &  +0.5_dp/xw*d13(1,2)*prop3456)*prop34
     &  +(q4*qe+l4*le*prop56)*d11(1,2)*prop3456
     &  +(q3*qe+l3*le*prop56)*d12(1,2)*prop3456

      amp(2,2)=
     &  ((Q(jd)*qe+L(jd)*re*prop56)*d123(2,2)
     &  +(Q(ju)*qe+L(ju)*re*prop56)*d456(2,2)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*re*prop56)*d78(2,2)*prop3456
     &  +(Q(2)*qe+R(2)*re*prop56)*d910(2,2)
     &  +0.5_dp/xw*d13(2,2)*prop3456)*prop34
     &  +(q4*qe+l4*re*prop56)*d11(2,2)*prop3456
     &  +(q3*qe+l3*re*prop56)*d12(2,2)*prop3456

      WZccmsq=fac*(abs(amp(1,1))**2+abs(amp(2,1))**2
     &            +abs(amp(1,2))**2+abs(amp(2,2))**2)
      return
      end
