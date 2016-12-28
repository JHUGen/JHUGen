      subroutine jonew(p7,p3,p4,p1,za,zb,zab,jZ,jg)
      implicit none
      include 'constants.f'
      include 'cmplxmass.f'
      include 'masses.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zcouple.f'
      include 'pid_pdg.f'
      include 'zacouplejk.f'
      integer p1,p2,p3,p4,p7,ro
      double complex zab(mxpart,4,mxpart),zab2,jZ(2,4),jg(2,4),
     & rxw,propw34,propw17
      double precision p17(4),p34(4),t3,s134,s347,
     & s34,s137,s147,s17,q3,q4,l3,l4,
     & Ldn,Lup,Qdn,Qup
C-----Begin statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s17=s(p1,p7)
      s34=s(p3,p4)
      s134=t3(p1,p3,p4)
      s137=t3(p1,p3,p7)
      s147=t3(p1,p4,p7)
      s347=t3(p3,p4,p7)
      p17(:)=0.5d0*dble(zab(p1,:,p1)+zab(p7,:,p7))
      p34(:)=0.5d0*dble(zab(p3,:,p3)+zab(p4,:,p4))

      propw17=s17-dcmplx(wmass**2,-wmass*wwidth)
      propw34=s34-dcmplx(wmass**2,-wmass*wwidth)

c--- determine correct couplings from calling parameters
      Ldn = L_jk(p1,p7,1)
      Lup = L_jk(p1,p7,2)
      Qdn = Q_jk(p1,p7,1)
      Qup = Q_jk(p1,p7,2)
      if (p3 .eq. 3) then
        q3=q1
        l3=l1
      else
        q3=q2
        l3=l2
      endif
      ! Find the complementary charges for the n3-n4 line
      if(pid_pdg(p3).eq.0 .and. pid_pdg(p4).eq.0) then
         ! Assign same couplings if both particles are unknown (jets).
         l4=l3
         q4=q3
      else if(
     & abs(pid_pdg(p3)).eq.11 .or.
     & abs(pid_pdg(p3)).eq.13 .or.
     & abs(pid_pdg(p3)).eq.15) then
         ! p3-p4 is an e current
         l4=ln
         q4=0d0
      else if(
     & abs(pid_pdg(p3)).eq.12 .or.
     & abs(pid_pdg(p3)).eq.14 .or.
     & abs(pid_pdg(p3)).eq.16) then
         ! p3-p4 is a nu current
         l4=le
         q4=qe
      else if(
     & abs(pid_pdg(p3)).eq.1 .or.
     & abs(pid_pdg(p3)).eq.3 .or.
     & abs(pid_pdg(p3)).eq.5) then
         ! p3-p4 is a d quark current
         l4=L(2)
         q4=Q(2)
      else if(
     & abs(pid_pdg(p3)).eq.2 .or.
     & abs(pid_pdg(p3)).eq.4) then
         ! p3-p4 is a u quark current
         l4=L(1)
         q4=Q(1)
      else
         ! p3-p4 case is unknown, kill the relevant amplitude.
         l4=0d0
         q4=0d0
      endif

      do ro=1,4
      jZ(1,ro)= + propw34**(-1) * ( Ldn*za(p7,p3)*zb(p7,p4)*zab(p7,ro,
     &    p1)*s347**(-1) + Ldn*za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*
     &    s347**(-1) - Lup*za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)*
     &    s134**(-1) + Lup*za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*
     &    s134**(-1) + za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1)*rxw -
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1)*rxw - zab2(p7,p3,p4
     &    ,p1)*zab(p3,ro,p4)*propw17**(-1)*rxw + zab2(p3,p1,p7,p4)*zab(
     &    p7,ro,p1)*propw17**(-1)*rxw )
      jZ(1,ro) = jZ(1,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*l3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*l4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*l4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*l3*s147**(-1)
      jZ(2,ro)= + propw34**(-1) * (  - Ldn*za(p1,p3)*zb(p1,p4)*zab(p7,
     &    ro,p1)*s134**(-1) + Ldn*za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*
     &    s134**(-1) + Lup*za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)*
     &    s347**(-1) + Lup*za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*
     &    s347**(-1) - za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1)*rxw +
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1)*rxw + zab2(p7,p3,p4
     &    ,p1)*zab(p3,ro,p4)*propw17**(-1)*rxw - zab2(p3,p1,p7,p4)*zab(
     &    p7,ro,p1)*propw17**(-1)*rxw )
      jZ(2,ro) = jZ(2,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*l3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*l4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*l4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*l3*s147**(-1)
      jg(1,ro)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)*
     &    Qdn*s347**(-1) + za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1) -
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1) + za(p7,p3)*zb(p3,
     &    p4)*zab(p3,ro,p1)*Qdn*s347**(-1) - za(p1,p3)*zb(p1,p4)*zab(p7,
     &    ro,p1)*Qup*s134**(-1) + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*Qup*
     &    s134**(-1) - zab2(p7,p3,p4,p1)*zab(p3,ro,p4)*propw17**(-1) +
     &    zab2(p3,p1,p7,p4)*zab(p7,ro,p1)*propw17**(-1) )
      jg(1,ro) = jg(1,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*q3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*q4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*q4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*q3*s147**(-1)
      jg(2,ro)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)*
     &    Qup*s347**(-1) - za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1) +
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1) + za(p7,p3)*zb(p3,
     &    p4)*zab(p3,ro,p1)*Qup*s347**(-1) - za(p1,p3)*zb(p1,p4)*zab(p7,
     &    ro,p1)*Qdn*s134**(-1) + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*Qdn*
     &    s134**(-1) + zab2(p7,p3,p4,p1)*zab(p3,ro,p4)*propw17**(-1) -
     &    zab2(p3,p1,p7,p4)*zab(p7,ro,p1)*propw17**(-1) )
      jg(2,ro) = jg(2,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*q3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*q4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*q4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*q3*s147**(-1)
      enddo
      jZ(:,:)=jZ(:,:)/cxw
      jg(:,:)=jg(:,:)/cxw
      return
      end
