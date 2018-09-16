      function A51_VH(j1,j2,j3,j4,j5,za,zb)

!  Amplitudes taken from Appendix IV of
!  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
!  %``One loop amplitudes for e+ e- to four partons,''
!  Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
!  [hep-ph/9708239].
!  Modified to remove momentum conservation relations
!===== C.Williams added missing terms which go like <3|P_125|4] in
!===== W + 1 jet and vanish for that case but are non zero for us

      implicit none
      include 'types.f'
      complex(dp):: A51_VH
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: Vcc,Fcc,Vsc,Fsc,l12,l23,L0,L1,Lsm1,A5lom
      complex(dp):: lnrat,zab2,zab3
      real(dp):: s123
      complex(dp):: MPC
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)
     &     +za(j1,j4)*zb(j4,j5)

!    -i * A5tree
!==== debug
!      epinv=0_dp
!      musq=1_dp

      s123=s(j1,j2)+s(j2,j3)+s(j3,j1)
      A5lom=-za(j3,j4)*zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))


!--leading N
      Vcc=
     & -(epinv**2+epinv*l12+half*l12**2)
     & -(epinv**2+epinv*l23+half*l23**2)
     & -two*(epinv+l23)-four

      Fcc=zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
     & *(za(j3,j4)*Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & +two*za(j3,j1)*zab2(j4,j2,j3,j1)
     &   *L0(-s(j2,j3),-s123)/s123)

      Vsc =half*(epinv+l23)+one
      Fsc =-za(j3,j4)*za(j3,j1)*zb(j1,j5)
     & /(za(j1,j2)*za(j2,j3))*L0(-s(j2,j3),-s123)/s123
     & +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/s123**2

!===== correction for piece
      MPC=-(3._dp*L0(-s(j2,j3),-s123)*za(j1,j3)*
     &   zab3(j4,j1,j2,j3,j5)*zb(j2,j1))/(2.*s123**2*za(j1,j2))

      A51_VH=(Vcc+Vsc)*A5Lom+Fcc+Fsc+MPC

!      write(6,*) 'LO LC ',A5LOm
!      write(6,*) 'NLO, epinv (2) ',-2_dp*A5LOm
!      write(6,*) 'NLO, epinv(1) ', (-l12-l23-2_dp+half)*A5LOm
!      write(6,*) 'NLO   ',A51_VH

!      write(6,*) 'Box coefficient ',zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123!)
!     & *(za(j3,j4))
!      write(6,*) 'rational coefficient', +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
!     & /(za(j1,j2)*za(j2,j3))/s123**2*(cone/cplx1(one-s(j2,j3)/s123))
!      pause


      return
      end
