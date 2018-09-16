      subroutine setupzprops(i1,i2,i3,i4,i5,i6,i7,i8,
     & gmZ7341,gmZ7561,gmZ17,gmZ28,
     & gmZ7342,gmZ7562,gmZ27,gmZ18,
     & ggWW,propw17,propw18,propw27,propw28,
     & propw7341,propw7561,propw7342,propw7562,
     & ll7341,ll7561,ll7342,ll7562,
     & gmZl7341,gmZl7561,gmZl7342,gmZl7562,
     & gmZl8562,gmZl8342,gmZl8561,gmZl8341)
      implicit none
      include 'types.f'
c--- Author: R.K. Ellis, November 2014
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
c      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
c      include 'runstring.f'
c      include 'masses.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7,i8,jdu1,jdu2
      real(dp):: t4,bit,s7341,s7561,s7342,s7562,
     & s34,s56,s17,s27,s18,s28
      complex(dp):: 
     & prop17,prop27,prop18,prop28,prop7341,prop7561,prop7342,prop7562,
     & prop34,prop56,
     & propw17,propw18,propw27,propw28,
     & propw7341,propw7561,propw7342,propw7562,
     & gmZ7341(2,2,2,2),gmZ7561(2,2,2,2),gmZ17(2,2,2,2),gmZ28(2,2,2,2),
     & gmZ7342(2,2,2,2),gmZ7562(2,2,2,2),gmZ27(2,2,2,2),gmZ18(2,2,2,2),
     & ll7341(2,2),ll7561(2,2),ll7342(2,2),ll7562(2,2),
     & gmZl7341(2,2,2),gmZl7561(2,2,2),gmZl7342(2,2,2),gmZl7562(2,2,2),
     & gmZl8562(2,2,2),gmZl8342(2,2,2),gmZl8561(2,2,2),gmZl8341(2,2,2),
     & ggWW(2,2),rxw

C---begin statement functions
      t4(i1,i2,i3,i4)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)
     & +s(i2,i3)+s(i2,i4)+s(i3,i4)
C---end statement functions

c      if (runstring(1:5) == 'check') then
c      bit=0d0
c      else
      bit=1d0
c      endif

      s7341=t4(i1,i7,i3,i4)
      s7561=t4(i1,i7,i5,i6)
      s7342=t4(i2,i7,i3,i4)
      s7562=t4(i2,i7,i5,i6)

      s34=s(i3,i4)
      s56=s(i5,i6)
      s17=s(i1,i7)
      s27=s(i2,i7)
      s18=s(i1,i8)
      s28=s(i2,i8)
      prop34=cplx1(s34)-czmass2
      prop56=cplx1(s56)-czmass2
      prop17=cplx1(s17)-czmass2
      prop27=cplx1(s27)-czmass2
      prop18=cplx1(s18)-czmass2
      prop28=cplx1(s28)-czmass2
      prop7341=cplx1(s7341)-czmass2
      prop7561=cplx1(s7561)-czmass2
      prop7342=cplx1(s7342)-czmass2
      prop7562=cplx1(s7562)-czmass2

      propw17=cplx1(s17)-cwmass2
      propw18=cplx1(s18)-cwmass2
      propw27=cplx1(s27)-cwmass2
      propw28=cplx1(s28)-cwmass2
      propw7341=cplx1(s7341)-cwmass2
      propw7561=cplx1(s7561)-cwmass2
      propw7342=cplx1(s7342)-cwmass2
      propw7562=cplx1(s7562)-cwmass2

C----setup couplings and propagators
      rxw=(cone-cxw)/cxw

      ggWW(1,1)=
     & cplx1(bit*qe**2/(s34*s56))+rxw*cplx1(le**2)/prop34/prop56
      ggWW(1,2)=
     & cplx1(bit*qe**2/(s34*s56))+rxw*cplx1(le*re)/prop34/prop56
      ggWW(2,1)=
     & cplx1(bit*qe**2/(s34*s56))+rxw*cplx1(re*le)/prop34/prop56
      ggWW(2,2)=
     & cplx1(bit*qe**2/(s34*s56))+rxw*cplx1(re**2)/prop34/prop56

      do jdu1=1,2
      do jdu2=1,2
      gmZ7341(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7341)+cplx1(L(jdu1)*L(jdu2))/prop7341
      gmZ7341(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7341)+cplx1(L(jdu1)*R(jdu2))/prop7341
      gmZ7341(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7341)+cplx1(R(jdu1)*L(jdu2))/prop7341
      gmZ7341(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7341)+cplx1(R(jdu1)*R(jdu2))/prop7341

      gmZ7561(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7561)+cplx1(L(jdu1)*L(jdu2))/prop7561
      gmZ7561(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7561)+cplx1(L(jdu1)*R(jdu2))/prop7561
      gmZ7561(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7561)+cplx1(R(jdu1)*L(jdu2))/prop7561
      gmZ7561(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7561)+cplx1(R(jdu1)*R(jdu2))/prop7561

      gmZ7342(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7342)+cplx1(L(jdu1)*L(jdu2))/prop7342
      gmZ7342(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7342)+cplx1(L(jdu1)*R(jdu2))/prop7342
      gmZ7342(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7342)+cplx1(R(jdu1)*L(jdu2))/prop7342
      gmZ7342(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7342)+cplx1(R(jdu1)*R(jdu2))/prop7342

      gmZ7562(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7562)+cplx1(L(jdu1)*L(jdu2))/prop7562
      gmZ7562(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7562)+cplx1(L(jdu1)*R(jdu2))/prop7562
      gmZ7562(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7562)+cplx1(R(jdu1)*L(jdu2))/prop7562
      gmZ7562(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s7562)+cplx1(R(jdu1)*R(jdu2))/prop7562

      gmZ17(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s17)+cplx1(L(jdu1)*L(jdu2))/prop17
      gmZ17(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s17)+cplx1(L(jdu1)*R(jdu2))/prop17
      gmZ17(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s17)+cplx1(R(jdu1)*L(jdu2))/prop17
      gmZ17(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s17)+cplx1(R(jdu1)*R(jdu2))/prop17

      gmZ27(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s27)+cplx1(L(jdu1)*L(jdu2))/prop27
      gmZ27(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s27)+cplx1(L(jdu1)*R(jdu2))/prop27
      gmZ27(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s27)+cplx1(R(jdu1)*L(jdu2))/prop27
      gmZ27(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s27)+cplx1(R(jdu1)*R(jdu2))/prop27

      gmZ18(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s18)+cplx1(L(jdu1)*L(jdu2))/prop18
      gmZ18(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s18)+cplx1(L(jdu1)*R(jdu2))/prop18
      gmZ18(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s18)+cplx1(R(jdu1)*L(jdu2))/prop18
      gmZ18(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s18)+cplx1(R(jdu1)*R(jdu2))/prop18

      gmZ28(jdu1,jdu2,1,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s28)+cplx1(L(jdu1)*L(jdu2))/prop28
      gmZ28(jdu1,jdu2,1,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s28)+cplx1(L(jdu1)*R(jdu2))/prop28
      gmZ28(jdu1,jdu2,2,1)=
     & cplx1(Q(jdu1)*Q(jdu2)/s28)+cplx1(R(jdu1)*L(jdu2))/prop28
      gmZ28(jdu1,jdu2,2,2)=
     & cplx1(Q(jdu1)*Q(jdu2)/s28)+cplx1(R(jdu1)*R(jdu2))/prop28
      enddo
      enddo

c--- lepton-lepton couplings
      ll7341(1,1)=
     & cplx1(q1*q2/s7341)+cplx1(l1*l2)/prop7341
      ll7341(1,2)=
     & cplx1(q1*q2/s7341)+cplx1(l1*r2)/prop7341
      ll7341(2,1)=
     & cplx1(q1*q2/s7341)+cplx1(r1*l2)/prop7341
      ll7341(2,2)=
     & cplx1(q1*q2/s7341)+cplx1(r1*r2)/prop7341

      ll7561(1,1)=
     & cplx1(q1*q2/s7561)+cplx1(l1*l2)/prop7561
      ll7561(1,2)=
     & cplx1(q1*q2/s7561)+cplx1(l1*r2)/prop7561
      ll7561(2,1)=
     & cplx1(q1*q2/s7561)+cplx1(r1*l2)/prop7561
      ll7561(2,2)=
     & cplx1(q1*q2/s7561)+cplx1(r1*r2)/prop7561

      ll7342(1,1)=
     & cplx1(q1*q2/s7342)+cplx1(l1*l2)/prop7342
      ll7342(1,2)=
     & cplx1(q1*q2/s7342)+cplx1(l1*r2)/prop7342
      ll7342(2,1)=
     & cplx1(q1*q2/s7342)+cplx1(r1*l2)/prop7342
      ll7342(2,2)=
     & cplx1(q1*q2/s7342)+cplx1(r1*r2)/prop7342

      ll7562(1,1)=
     & cplx1(q1*q2/s7562)+cplx1(l1*l2)/prop7562
      ll7562(1,2)=
     & cplx1(q1*q2/s7562)+cplx1(l1*r2)/prop7562
      ll7562(2,1)=
     & cplx1(q1*q2/s7562)+cplx1(r1*l2)/prop7562
      ll7562(2,2)=
     & cplx1(q1*q2/s7562)+cplx1(r1*r2)/prop7562

c--- gamma/Z-lepton couplings
      do jdu1=1,2
      gmZl7341(jdu1,1,1)=
     & cplx1(Q(jdu1)*q1/s7341)+cplx1(L(jdu1)*l1)/prop7341
      gmZl7341(jdu1,1,2)=
     & cplx1(Q(jdu1)*q1/s7341)+cplx1(L(jdu1)*r1)/prop7341
      gmZl7341(jdu1,2,1)=
     & cplx1(Q(jdu1)*q1/s7341)+cplx1(R(jdu1)*l1)/prop7341
      gmZl7341(jdu1,2,2)=
     & cplx1(Q(jdu1)*q1/s7341)+cplx1(R(jdu1)*r1)/prop7341

      gmZl7561(jdu1,1,1)=
     & cplx1(Q(jdu1)*q2/s7561)+cplx1(L(jdu1)*l2)/prop7561
      gmZl7561(jdu1,1,2)=
     & cplx1(Q(jdu1)*q2/s7561)+cplx1(L(jdu1)*r2)/prop7561
      gmZl7561(jdu1,2,1)=
     & cplx1(Q(jdu1)*q2/s7561)+cplx1(R(jdu1)*l2)/prop7561
      gmZl7561(jdu1,2,2)=
     & cplx1(Q(jdu1)*q2/s7561)+cplx1(R(jdu1)*r2)/prop7561

      gmZl7342(jdu1,1,1)=
     & cplx1(Q(jdu1)*q1/s7342)+cplx1(L(jdu1)*l1)/prop7342
      gmZl7342(jdu1,1,2)=
     & cplx1(Q(jdu1)*q1/s7342)+cplx1(L(jdu1)*r1)/prop7342
      gmZl7342(jdu1,2,1)=
     & cplx1(Q(jdu1)*q1/s7342)+cplx1(R(jdu1)*l1)/prop7342
      gmZl7342(jdu1,2,2)=
     & cplx1(Q(jdu1)*q1/s7342)+cplx1(R(jdu1)*r1)/prop7342

      gmZl7562(jdu1,1,1)=
     & cplx1(Q(jdu1)*q2/s7562)+cplx1(L(jdu1)*l2)/prop7562
      gmZl7562(jdu1,1,2)=
     & cplx1(Q(jdu1)*q2/s7562)+cplx1(L(jdu1)*r2)/prop7562
      gmZl7562(jdu1,2,1)=
     & cplx1(Q(jdu1)*q2/s7562)+cplx1(R(jdu1)*l2)/prop7562
      gmZl7562(jdu1,2,2)=
     & cplx1(Q(jdu1)*q2/s7562)+cplx1(R(jdu1)*r2)/prop7562
      enddo
      
c--- gamma/Z-lepton couplings
      do jdu1=1,2
      gmZl8562(jdu1,1,1)=
     & cplx1(Q(jdu1)*q2/s7341)+cplx1(L(jdu1)*l2)/prop7341
      gmZl8562(jdu1,1,2)=
     & cplx1(Q(jdu1)*q2/s7341)+cplx1(L(jdu1)*r2)/prop7341
      gmZl8562(jdu1,2,1)=
     & cplx1(Q(jdu1)*q2/s7341)+cplx1(R(jdu1)*l2)/prop7341
      gmZl8562(jdu1,2,2)=
     & cplx1(Q(jdu1)*q2/s7341)+cplx1(R(jdu1)*r2)/prop7341

      gmZl8342(jdu1,1,1)=
     & cplx1(Q(jdu1)*q1/s7561)+cplx1(L(jdu1)*l1)/prop7561
      gmZl8342(jdu1,1,2)=
     & cplx1(Q(jdu1)*q1/s7561)+cplx1(L(jdu1)*r1)/prop7561
      gmZl8342(jdu1,2,1)=
     & cplx1(Q(jdu1)*q1/s7561)+cplx1(R(jdu1)*l1)/prop7561
      gmZl8342(jdu1,2,2)=
     & cplx1(Q(jdu1)*q1/s7561)+cplx1(R(jdu1)*r1)/prop7561

      gmZl8561(jdu1,1,1)=
     & cplx1(Q(jdu1)*q2/s7342)+cplx1(L(jdu1)*l2)/prop7342
      gmZl8561(jdu1,1,2)=
     & cplx1(Q(jdu1)*q2/s7342)+cplx1(L(jdu1)*r2)/prop7342
      gmZl8561(jdu1,2,1)=
     & cplx1(Q(jdu1)*q2/s7342)+cplx1(R(jdu1)*l2)/prop7342
      gmZl8561(jdu1,2,2)=
     & cplx1(Q(jdu1)*q2/s7342)+cplx1(R(jdu1)*r2)/prop7342

      gmZl8341(jdu1,1,1)=
     & cplx1(Q(jdu1)*q1/s7562)+cplx1(L(jdu1)*l1)/prop7562
      gmZl8341(jdu1,1,2)=
     & cplx1(Q(jdu1)*q1/s7562)+cplx1(L(jdu1)*r1)/prop7562
      gmZl8341(jdu1,2,1)=
     & cplx1(Q(jdu1)*q1/s7562)+cplx1(R(jdu1)*l1)/prop7562
      gmZl8341(jdu1,2,2)=
     & cplx1(Q(jdu1)*q1/s7562)+cplx1(R(jdu1)*r1)/prop7562
      enddo
      
      return
      end
