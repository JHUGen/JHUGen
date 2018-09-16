      subroutine ampZWZ(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer:: jdu1,jdu2,h17,h28,i1,i2,i3,i4,i5,i6,i7,i8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: amp(2,2,2,2),propZ17,propZ28,propW1347,
     & gamZ17lept3(2,2),gamZ17anti4(2,2),
     & gamZ28lept5(2,2),gamZ28anti6(2,2)
      real(dp):: qn,t3,t4,s17,s28,s137,s147,s258,s268,
     & q3,l3,q4,l4,q5,l5,q6,l6
C     amp(jdu1,jdu2,h17,h28)
      parameter(qn=0d0)
C-----Begin statement functions
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &               +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      s17=s(i1,i7)
      s28=s(i2,i8)
      s137=t3(i1,i3,i7)
      s147=t3(i1,i4,i7)
      s258=t3(i2,i5,i8)
      s268=t3(i2,i6,i8)

c--- set up lepton couplings according to call
      if (i3 == 3) then
        q3=qe
        l3=le
        q4=qn
        l4=ln
        q5=qn
        l5=ln
        q6=qe
        l6=le
      else
        q3=qn
        l3=ln
        q4=qe
        l4=le
        q5=qe
        l5=le
        q6=qn
        l6=ln
      endif

      propz17=s17-czmass2
      propz28=s28-czmass2
      propw1347=t4(i1,i3,i4,i7)-cwmass2

      do jdu1=1,2
      gamZ17lept3(jdu1,1)=Q(jdu1)*q3/s17+L(jdu1)*l3/propZ17
      gamZ17lept3(jdu1,2)=Q(jdu1)*q3/s17+R(jdu1)*l3/propZ17

      gamZ17anti4(jdu1,1)=Q(jdu1)*q4/s17+L(jdu1)*l4/propZ17
      gamZ17anti4(jdu1,2)=Q(jdu1)*q4/s17+R(jdu1)*l4/propZ17

      gamZ28lept5(jdu1,1)=Q(jdu1)*q5/s28+L(jdu1)*l5/propZ28
      gamZ28lept5(jdu1,2)=Q(jdu1)*q5/s28+R(jdu1)*l5/propZ28

      gamZ28anti6(jdu1,1)=Q(jdu1)*q6/s28+L(jdu1)*l6/propZ28
      gamZ28anti6(jdu1,2)=Q(jdu1)*q6/s28+R(jdu1)*l6/propZ28
      enddo

      p3=i3
      p4=i4
      p5=i5
      p6=i6

      do h17=1,2
      do h28=1,2
      if (h17 == 1) then
      p1=i1
      p7=i7
      elseif (h17 == 2) then
      p1=i7
      p7=i1
      endif
      if (h28 == 1) then
      p2=i2
      p8=i8
      elseif (h28 == 2) then
      p2=i8
      p8=i2
      endif
      do jdu1=1,2
      do jdu2=1,2
      amp(jdu1,jdu2,h17,h28)= + gamZ17lept3(jdu1,h17)*gamZ28lept5(jdu2,
     & h28)*propW1347**(-1)*s137**(-1)*s258**(-1) * (  - 4.D0*za(p7,p3)
     &    *za(p7,p8)*za(p8,p5)*zb(p7,p1)*zb(p4,p6)*zb(p8,p2) + 4.D0*za(
     &    p7,p3)*za(p7,p5)*za(p8,p5)*zb(p7,p1)*zb(p4,p6)*zb(p2,p5) + 4.D
     &    0*za(p7,p3)*za(p3,p8)*za(p8,p5)*zb(p1,p3)*zb(p4,p6)*zb(p8,p2)
     &     - 4.D0*za(p7,p3)*za(p3,p5)*za(p8,p5)*zb(p1,p3)*zb(p4,p6)*zb(
     &    p2,p5) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17lept3(
     & jdu1,h17)*gamZ28lept5(jdu2,h28)*propW1347**(-1) * ( 2.D0*za(p7,
     &    p3)*za(p8,p5)*zb(p1,p4)*zb(p2,p6)*cwmass2**(-1) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17lept3(
     & jdu1,h17)*gamZ28anti6(jdu2,h28)*propW1347**(-1)*s137**(-1)*
     & s268**(-1) * (  - 4.D0*za(p7,p3)*za(p7,p5)*za(p8,p2)*zb(p7,p1)*
     &    zb(p4,p2)*zb(p2,p6) - 4.D0*za(p7,p3)*za(p7,p5)*za(p8,p6)*zb(
     &    p7,p1)*zb(p4,p6)*zb(p2,p6) + 4.D0*za(p7,p3)*za(p3,p5)*za(p8,
     &    p2)*zb(p1,p3)*zb(p4,p2)*zb(p2,p6) + 4.D0*za(p7,p3)*za(p3,p5)*
     &    za(p8,p6)*zb(p1,p3)*zb(p4,p6)*zb(p2,p6) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17lept3(
     & jdu1,h17)*gamZ28anti6(jdu2,h28)*propW1347**(-1) * (  - 2.D0*za(
     &    p7,p3)*za(p8,p5)*zb(p1,p4)*zb(p2,p6)*cwmass2**(-1) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17anti4(
     & jdu1,h17)*gamZ28lept5(jdu2,h28)*propW1347**(-1)*s147**(-1)*
     & s258**(-1) * (  - 4.D0*za(p7,p1)*za(p3,p8)*za(p8,p5)*zb(p1,p4)*
     &    zb(p1,p6)*zb(p8,p2) + 4.D0*za(p7,p1)*za(p3,p5)*za(p8,p5)*zb(
     &    p1,p4)*zb(p1,p6)*zb(p2,p5) - 4.D0*za(p7,p4)*za(p3,p8)*za(p8,
     &    p5)*zb(p1,p4)*zb(p4,p6)*zb(p8,p2) + 4.D0*za(p7,p4)*za(p3,p5)*
     &    za(p8,p5)*zb(p1,p4)*zb(p4,p6)*zb(p2,p5) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17anti4(
     & jdu1,h17)*gamZ28lept5(jdu2,h28)*propW1347**(-1) * (  - 2.D0*za(
     &    p7,p3)*za(p8,p5)*zb(p1,p4)*zb(p2,p6)*cwmass2**(-1) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17anti4(
     & jdu1,h17)*gamZ28anti6(jdu2,h28)*propW1347**(-1)*s147**(-1)*
     & s268**(-1) * (  - 4.D0*za(p7,p1)*za(p3,p5)*za(p8,p2)*zb(p1,p4)*
     &    zb(p1,p2)*zb(p2,p6) - 4.D0*za(p7,p1)*za(p3,p5)*za(p8,p6)*zb(
     &    p1,p4)*zb(p1,p6)*zb(p2,p6) - 4.D0*za(p7,p4)*za(p3,p5)*za(p8,
     &    p2)*zb(p1,p4)*zb(p4,p2)*zb(p2,p6) - 4.D0*za(p7,p4)*za(p3,p5)*
     &    za(p8,p6)*zb(p1,p4)*zb(p4,p6)*zb(p2,p6) )
      amp(jdu1,jdu2,h17,h28) = amp(jdu1,jdu2,h17,h28) + gamZ17anti4(
     & jdu1,h17)*gamZ28anti6(jdu2,h28)*propW1347**(-1) * ( 2.D0*za(p7,
     &    p3)*za(p8,p5)*zb(p1,p4)*zb(p2,p6)*cwmass2**(-1) )
      enddo
      enddo
      enddo
      enddo
      amp(:,:,:,:)=amp(:,:,:,:)/cxw
      return
      end
