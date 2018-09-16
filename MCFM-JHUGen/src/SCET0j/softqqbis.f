!#################################################################
!##### O(as) soft contribution to tau_1, calculates soft1(-1:1)
!##### contributions to delta(tau1), L0(tau1) and L1(tau1)
!##### returns integrated contribution.
!##### a factor of alphas/2/pi has been extracted out
!##### see Eq. (2.22) of 1012.4480,
!#################################################################
      subroutine softqqbis(order,soft1,soft2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'nf.f'
      include 'scet_const.f'
      integer :: order
      real(dp) :: soft1(-1:1),soft2(-1:3)

      if (order < 1) return
c--- coefficient of (alphas/2/pi)
!---- contribution to delta(tau1)
      soft1(-1) = CF*pisqo6
!---- contribution to L0(tau1)
      soft1(0) = 0.0_dp
!---- contribution to L1(tau1)
      soft1(1) =-8*CF


!---- L taucut
c      Lt1cut = log(taucut/scale)
c      soft(1)=soft1(-1)
c     & +soft1(0)*Lt1cut
c     & +soft1(1)*Lt1cut**2/two

      if (order < 2) return
c--- coefficient of (alphas/2/pi)^2
c--- note overall factor of 4*CF from Eq. (2.22) compared with Eq. (2.23)

!---- contribution to delta(tau1)
      soft2(-1) =
     & +CF**2*( - 27/ten*zeta2**2)
     & +CF*(- 5/27._dp*be0 - 160/27._dp*CA - 37/12._dp*be0*zeta2
     &      + 2*CA*zeta2 + 22/five*CA*zeta2**2 + 58/12._dp*be0*zeta3)
!---- contribution to L0(tau1)
      soft2(0) = 
     & + CF**2*64*zeta3
     & + CF*(28/nine*be0 - 2*be0*zeta2 + 32/nine*CA - 14*CA*zeta3)
!---- contribution to L1(tau1)
      soft2(1) = 
     & - CF**2*72*zeta2
     & + CF*(- 20/three*be0 - 16/three*CA + 8*CA*zeta2)
!---- contribution to L2(tau1)
      soft2(2) = 4*CF*be0
!---- contribution to L3(tau1)
      soft2(3) = 32*CF**2

c      soft(2)=soft2(-1)
c     & +soft2(0)*Lt1cut
c     & +soft2(1)*Lt1cut**2/two
c     & +soft2(2)*Lt1cut**3/three
c     & +soft2(3)*Lt1cut**4/four
      return
      end


