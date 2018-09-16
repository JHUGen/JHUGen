      subroutine jww(p7,p3,p4,p1,za,zb,zab,jw)
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
      integer:: jdu1,p1,p2,p3,p4,p7,ro
      complex(dp):: zab(mxpart,4,mxpart),zab2,zba2,jw(4,2,2),
     & rxw,propw34,propz17,gamv(2,2),gamvl3(2,2),gamvl4(2,2)
      real(dp):: qn,p17(4),p34(4),t3,s134,s347,
     & s34,s137,s147,s17,q3,q4,l3,l4
C     amp(jdu1,jdu2,h17,h28)
      parameter(qn=0d0)
C-----Begin statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s17=s(p1,p7)
      s34=s(p3,p4)
      s134=t3(p1,p3,p4)
      s137=t3(p1,p3,p7)
      s147=t3(p1,p4,p7)
      s347=t3(p3,p4,p7)
      p17(:)=0.5d0*real(zab(p1,:,p1)+zab(p7,:,p7))
      p34(:)=0.5d0*real(zab(p3,:,p3)+zab(p4,:,p4))

      propz17=s17-czmass2
      propw34=s34-cwmass2

      q3=qn
      l3=ln
      q4=qe
      l4=le
      do jdu1=1,2
      gamV(jdu1,1)=Q(jdu1)/s17+L(jdu1)*rxw/propz17
      gamV(jdu1,2)=Q(jdu1)/s17+R(jdu1)*rxw/propz17
      gamvl3(jdu1,1)=Q(jdu1)*q3/s17+L(jdu1)*l3/propz17
      gamvl3(jdu1,2)=Q(jdu1)*q3/s17+R(jdu1)*l3/propz17
      gamvl4(jdu1,1)=Q(jdu1)*q4/s17+L(jdu1)*l4/propz17
      gamvl4(jdu1,2)=Q(jdu1)*q4/s17+R(jdu1)*l4/propz17
      enddo

      do ro=1,4
      jw(ro,1,1)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &    *cxw**(-1)*s347**(-1) + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*
     &    cxw**(-1)*s347**(-1) )
      jw(ro,1,1) = jw(ro,1,1) + gamv(1,1)*propw34**(-1) * ( 2.D0*za(p7,
     &    p3)*zb(p1,p4)*p17(ro) - 2.D0*za(p7,p3)*zb(p1,p4)*p34(ro) - 2.D
     &    0*zab2(p7,p3,p4,p1)*zab(p3,ro,p4) + 2.D0*zab2(p3,p1,p7,p4)*
     &    zab(p7,ro,p1) )
      jw(ro,1,1) = jw(ro,1,1) - 2.D0*gamvl3(1,1)*za(p7,p3)*zb(p7,p1)*
     & zab(p7,ro,p4)*s137**(-1) + 2.D0*gamvl3(1,1)*za(p7,p3)*zb(p1,p3)*
     &    zab(p3,ro,p4)*s137**(-1) - 2.D0*gamvl4(1,1)*za(p7,p1)*zb(p1,
     &    p4)*zab(p3,ro,p1)*s147**(-1) - 2.D0*gamvl4(1,1)*za(p7,p4)*zb(
     &    p1,p4)*zab(p3,ro,p4)*s147**(-1)
      jw(ro,1,2)= + gamv(1,2)*propw34**(-1) * ( 2.D0*za(p1,p3)*zb(p7,p4
     &    )*p17(ro) - 2.D0*za(p1,p3)*zb(p7,p4)*p34(ro) + 2.D0*zab2(p3,
     &    p1,p7,p4)*zab(p1,ro,p7) - 2.D0*zba2(p7,p3,p4,p1)*zab(p3,ro,p4
     &    ) )
      jw(ro,1,2) = jw(ro,1,2) + 2.D0*gamvl3(1,2)*za(p1,p3)*zb(p7,p1)*
     & zab(p1,ro,p4)*s137**(-1) + 2.D0*gamvl3(1,2)*za(p1,p3)*zb(p7,p3)*
     &    zab(p3,ro,p4)*s137**(-1) + 2.D0*gamvl4(1,2)*za(p7,p1)*zb(p7,
     &    p4)*zab(p3,ro,p7)*s147**(-1) - 2.D0*gamvl4(1,2)*za(p1,p4)*zb(
     &    p7,p4)*zab(p3,ro,p4)*s147**(-1)
      jw(ro,2,1)= + propw34**(-1) * (  - za(p1,p3)*zb(p1,p4)*zab(p7,ro,
     &    p1)*cxw**(-1)*s134**(-1) + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*
     &    cxw**(-1)*s134**(-1) )
      jw(ro,2,1) = jw(ro,2,1) + gamv(2,1)*propw34**(-1) * ( 2.D0*za(p7,
     &    p3)*zb(p1,p4)*p17(ro) - 2.D0*za(p7,p3)*zb(p1,p4)*p34(ro) - 2.D
     &    0*zab2(p7,p3,p4,p1)*zab(p3,ro,p4) + 2.D0*zab2(p3,p1,p7,p4)*
     &    zab(p7,ro,p1) )
      jw(ro,2,1) = jw(ro,2,1) - 2.D0*gamvl3(2,1)*za(p7,p3)*zb(p7,p1)*
     & zab(p7,ro,p4)*s137**(-1) + 2.D0*gamvl3(2,1)*za(p7,p3)*zb(p1,p3)*
     &    zab(p3,ro,p4)*s137**(-1) - 2.D0*gamvl4(2,1)*za(p7,p1)*zb(p1,
     &    p4)*zab(p3,ro,p1)*s147**(-1) - 2.D0*gamvl4(2,1)*za(p7,p4)*zb(
     &    p1,p4)*zab(p3,ro,p4)*s147**(-1)
      jw(ro,2,2)= + gamv(2,2)*propw34**(-1) * ( 2.D0*za(p1,p3)*zb(p7,p4
     &    )*p17(ro) - 2.D0*za(p1,p3)*zb(p7,p4)*p34(ro) + 2.D0*zab2(p3,
     &    p1,p7,p4)*zab(p1,ro,p7) - 2.D0*zba2(p7,p3,p4,p1)*zab(p3,ro,p4
     &    ) )
      jw(ro,2,2) = jw(ro,2,2) + 2.D0*gamvl3(2,2)*za(p1,p3)*zb(p7,p1)*
     & zab(p1,ro,p4)*s137**(-1) + 2.D0*gamvl3(2,2)*za(p1,p3)*zb(p7,p3)*
     &    zab(p3,ro,p4)*s137**(-1) + 2.D0*gamvl4(2,2)*za(p7,p1)*zb(p7,
     &    p4)*zab(p3,ro,p7)*s147**(-1) - 2.D0*gamvl4(2,2)*za(p1,p4)*zb(
     &    p7,p4)*zab(p3,ro,p4)*s147**(-1)
      enddo
      jw(:,:,:)=jw(:,:,:)/sqrt(2d0*cxw)
      return
      end
