      function A51(j1,j2,j3,j4,j5,za,zb)

!  Amplitudes taken from Appendix IV of
!  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
!  %``One loop amplitudes for e+ e- to four partons,''
!  Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
!  [hep-ph/9708239].
!  Modified to remove momentum conservation relations
      implicit none
      include 'types.f'
      complex(dp):: A51

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: Vcc,Fcc,Vsc,Fsc,l12,l23,L0,L1,Lsm1,A5lom
      complex(dp):: lnrat,zab2
      real(dp):: s123

      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

!    -i * A5tree  Eq.(IV.1)
      s123=s(j1,j2)+s(j2,j3)+s(j3,j1)
      A5lom=-za(j3,j4)*zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)

!--leading N  Eq. (IV.2)
      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))
      Vcc=
     & -(epinv**2+epinv*l12+half*l12**2)
     & -(epinv**2+epinv*l23+half*l23**2)
     & -two*(epinv+l23)-four

!--Eq. (IV.3)
      Fcc=zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
     & *(za(j3,j4)*Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & +two*za(j3,j1)*zab2(j4,j2,j3,j1)
     &   *L0(-s(j2,j3),-s123)/s123)

!--Eq. (IV.4)
      Vsc =half*(epinv+l23)+one
!--Eq. (IV.5)
      Fsc =-za(j3,j4)*za(j3,j1)*zb(j1,j5)
     & /(za(j1,j2)*za(j2,j3))*L0(-s(j2,j3),-s123)/s123
     & +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/s123**2

      A51=(Vcc+Vsc)*A5Lom+Fcc+Fsc

      return
      end
