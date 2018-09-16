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
      logical isALepton,isANeutrino,isUpTypeLightQuark,isDnTypeQuark,
     & isAnUnknownJet
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
      ! Find the charges for the n3-n4 line
      if(
     & (
     & isDnTypeQuark(abs(pid_pdg(p3))) .or.
     & isUpTypeLightQuark(abs(pid_pdg(p4)))
     & ) .or. isALepton(abs(pid_pdg(p3)))
     & ) then
         ! p3-p4 is a d/e current with negative charge
         l3=L_jk(p3,p4,1)
         l4=L_jk(p3,p4,2)
         q3=Q_jk(p3,p4,1)
         q4=Q_jk(p3,p4,2)
      else if(
     & (
     & isUpTypeLightQuark(abs(pid_pdg(p3))) .or.
     & isDnTypeQuark(abs(pid_pdg(p4)))
     & ) .or. isANeutrino(abs(pid_pdg(p3)))
     & ) then
         ! p3-p4 is a u/nu current with non-negative charge
         l3=L_jk(p3,p4,2)
         l4=L_jk(p3,p4,1)
         q3=Q_jk(p3,p4,2)
         q4=Q_jk(p3,p4,1)
      else if(
     & isAnUnknownJet(pid_pdg(p3)) .and.
     & isAnUnknownJet(pid_pdg(p4))
     & ) then
         ! Assign same couplings if both particles are unknown (jets).
         l3=0d0
         q3=0d0
         if(p3.eq.3) then
            l3=l1
            q3=q1
         else if(p3.eq.5) then
            l3=l2
            q3=q2
         endif
         l4=l3
         q4=q3
      else
         ! p3-p4 case is unknown, kill the relevant amplitude.
         l3=0d0
         q3=0d0
         l4=0d0
         q4=0d0
      endif

      do ro=1,4
      jZ(1,ro)= propw34**(-1) * (
     &   + Ldn*s347**(-1)*(
     &   + za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &   + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)
     &   )
     &   + Lup*s134**(-1)*(
     &   - za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)
     &   + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)
     &   )
     &   + propw17**(-1)*rxw*(
     &   + za(p7,p3)*zb(p1,p4)*p17(ro)
     &   - za(p7,p3)*zb(p1,p4)*p34(ro)
     &   - zab2(p7,p3,p4,p1)*zab(p3,ro,p4)
     &   + zab2(p3,p1,p7,p4)*zab(p7,ro,p1)
     &   )
     & )
      jZ(1,ro) = jZ(1,ro) + propw17**(-1)*(
     &   + l3*s147**(-1)*(
     &   - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)
     &   - za(p7,p4)*zb(p1,p4)*zab(p3,ro,p4)
     &   )
     &   + l4*s137**(-1)*(
     &   - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)
     &   + za(p7,p3)*zb(p1,p3)*zab(p3,ro,p4)
     &   )
     & )
      jZ(2,ro)= + propw34**(-1) * (
     &   + Lup*s347**(-1)*(
     &   + za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &   + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)
     &   )
     &   + Ldn*s134**(-1)*(
     &   - za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)
     &   + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)
     &   )
     &   + propw17**(-1)*rxw*(
     &   - za(p7,p3)*zb(p1,p4)*p17(ro)
     &   + za(p7,p3)*zb(p1,p4)*p34(ro)
     &   + zab2(p7,p3,p4,p1)*zab(p3,ro,p4)
     &   - zab2(p3,p1,p7,p4)*zab(p7,ro,p1)
     &   )
     & )
      jZ(2,ro) = jZ(2,ro) + propw17**(-1)*(
     &   + l3*s147**(-1)*(
     &   - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)
     &   - za(p7,p4)*zb(p1,p4)*zab(p3,ro,p4)
     &   )
     &   + l4*s137**(-1)*(
     &   - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)
     &   + za(p7,p3)*zb(p1,p3)*zab(p3,ro,p4)
     &   )
     & )

      jg(1,ro)= propw34**(-1) * (
     &   + Qdn*s347**(-1)*(
     &   + za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &   + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)
     &   )
     &   + Qup*s134**(-1)*(
     &   - za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)
     &   + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)
     &   )
     &   + propw17**(-1)*(
     &   + za(p7,p3)*zb(p1,p4)*p17(ro)
     &   - za(p7,p3)*zb(p1,p4)*p34(ro)
     &   - zab2(p7,p3,p4,p1)*zab(p3,ro,p4)
     &   + zab2(p3,p1,p7,p4)*zab(p7,ro,p1)
     &   )
     & )
      jg(1,ro) = jg(1,ro) + propw17**(-1)*(
     &   + q3*s147**(-1)*(
     &   - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)
     &   - za(p7,p4)*zb(p1,p4)*zab(p3,ro,p4)
     &   )
     &   + q4*s137**(-1)*(
     &   - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)
     &   + za(p7,p3)*zb(p1,p3)*zab(p3,ro,p4)
     &   )
     & )
      jg(2,ro)= + propw34**(-1) * (
     &   + Qup*s347**(-1)*(
     &   + za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &   + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)
     &   )
     &   + Qdn*s134**(-1)*(
     &   - za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)
     &   + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)
     &   )
     &   + propw17**(-1)*(
     &   - za(p7,p3)*zb(p1,p4)*p17(ro)
     &   + za(p7,p3)*zb(p1,p4)*p34(ro)
     &   + zab2(p7,p3,p4,p1)*zab(p3,ro,p4)
     &   - zab2(p3,p1,p7,p4)*zab(p7,ro,p1)
     &   )
     & )
      jg(2,ro) = jg(2,ro) + propw17**(-1)*(
     &   + q3*s147**(-1)*(
     &   - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)
     &   - za(p7,p4)*zb(p1,p4)*zab(p3,ro,p4)
     &   )
     &   + q4*s137**(-1)*(
     &   - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)
     &   + za(p7,p3)*zb(p1,p3)*zab(p3,ro,p4)
     &   )
     & )

      enddo
      jZ(:,:)=jZ(:,:)/cxw
      jg(:,:)=jg(:,:)/cxw
      return
      end
