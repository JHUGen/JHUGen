!===== T. Dennen, May 2014
!===== Initialise array of triangle integral values.
!===== For use with qqb->4gamma and qqb->2j2gamma
!===== coefficients by same author.
      subroutine Triint_init(i1,i2,i3,i4,i5,i6,Triint,ord)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'scale.f'
      integer:: i1,i2,i3,i4,i5,i6
      integer:: ord
      complex(dp):: Triint(90), qlI3
      real(dp):: t

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      Triint(:) = czip

      Triint(1) = qlI3(0d0,0d0,s(i1,i2),0d0,0d0,0d0,musq,ord)
      Triint(2) = qlI3(0d0,s(i1,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(3) = qlI3(s(i5,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(4) = qlI3(0d0,s(i1,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(5) = qlI3(s(i4,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(6) = qlI3(0d0,s(i1,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(7) = qlI3(s(i4,i5),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(8) = qlI3(s(i3,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(9) = qlI3(0d0,s(i1,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(10) = qlI3(s(i3,i5),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(11) = qlI3(s(i3,i4),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(12) = qlI3(s(i2,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(13) = qlI3(s(i2,i5),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(14) = qlI3(s(i2,i4),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(15) = qlI3(s(i2,i3),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(16) = qlI3(0d0,s(i2,i3),t(i4,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(17) = qlI3(0d0,t(i2,i3,i4),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(18) = qlI3(s(i1,i2),0d0,t(i4,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(19) = qlI3(s(i1,i2),t(i3,i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(20) = qlI3(t(i1,i2,i3),0d0,s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(21) = qlI3(t(i1,i2,i3),s(i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(22) = qlI3(s(i1,i2),t(i3,i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(23) = qlI3(t(i1,i2,i3),s(i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(24) = qlI3(0d0,t(i2,i3,i5),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(25) = qlI3(s(i1,i2),t(i3,i5,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(26) = qlI3(0d0,t(i2,i3,i6),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(27) = qlI3(0d0,s(i2,i4),t(i3,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(28) = qlI3(t(i1,i2,i4),0d0,s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(29) = qlI3(t(i1,i2,i4),s(i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(30) = qlI3(t(i1,i2,i4),s(i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(31) = qlI3(0d0,t(i2,i4,i5),s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(32) = qlI3(0d0,t(i2,i4,i6),s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(33) = qlI3(0d0,s(i2,i5),t(i3,i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(34) = qlI3(t(i1,i2,i5),0d0,s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(35) = qlI3(t(i1,i2,i5),s(i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(36) = qlI3(t(i1,i2,i5),s(i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(37) = qlI3(0d0,t(i2,i5,i6),s(i3,i4),0d0,0d0,0d0,musq,ord)
      Triint(38) = qlI3(0d0,s(i2,i6),t(i3,i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(39) = qlI3(t(i1,i2,i6),0d0,s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(40) = qlI3(t(i1,i2,i6),s(i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(41) = qlI3(t(i1,i2,i6),s(i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(42) = qlI3(s(i1,i3),0d0,t(i4,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(43) = qlI3(s(i1,i3),t(i2,i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(44) = qlI3(s(i1,i3),t(i2,i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(45) = qlI3(s(i1,i3),t(i2,i5,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(46) = qlI3(t(i1,i3,i4),0d0,s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(47) = qlI3(t(i1,i3,i4),s(i2,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(48) = qlI3(t(i1,i3,i4),s(i2,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(49) = qlI3(t(i1,i3,i5),0d0,s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(50) = qlI3(t(i1,i3,i5),s(i2,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(51) = qlI3(t(i1,i3,i5),s(i2,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(52) = qlI3(t(i1,i3,i6),0d0,s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(53) = qlI3(t(i1,i3,i6),s(i2,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(54) = qlI3(t(i1,i3,i6),s(i2,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(55) = qlI3(s(i1,i4),0d0,t(i3,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(56) = qlI3(s(i1,i4),t(i2,i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(57) = qlI3(s(i1,i4),t(i2,i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(58) = qlI3(s(i1,i4),t(i2,i5,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(59) = qlI3(t(i1,i4,i5),0d0,s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(60) = qlI3(t(i1,i4,i5),s(i2,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(61) = qlI3(t(i1,i4,i5),s(i2,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(62) = qlI3(t(i1,i4,i6),0d0,s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(63) = qlI3(t(i1,i4,i6),s(i2,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(64) = qlI3(t(i1,i4,i6),s(i2,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(65) = qlI3(s(i1,i5),0d0,t(i3,i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(66) = qlI3(s(i1,i5),t(i2,i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(67) = qlI3(s(i1,i5),t(i2,i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(68) = qlI3(s(i1,i5),t(i2,i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(69) = qlI3(t(i1,i5,i6),0d0,s(i3,i4),0d0,0d0,0d0,musq,ord)
      Triint(70) = qlI3(t(i1,i5,i6),s(i2,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(71) = qlI3(t(i1,i5,i6),s(i2,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(72) = qlI3(s(i1,i6),0d0,t(i3,i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(73) = qlI3(s(i1,i6),t(i2,i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(74) = qlI3(s(i1,i6),t(i2,i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(75) = qlI3(s(i1,i6),t(i2,i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(76) = qlI3(s(i1,i2),s(i3,i4),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(77) = qlI3(s(i1,i2),s(i3,i5),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(78) = qlI3(s(i1,i2),s(i3,i6),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(79) = qlI3(s(i1,i3),s(i2,i4),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(80) = qlI3(s(i1,i3),s(i2,i5),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(81) = qlI3(s(i1,i3),s(i2,i6),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(82) = qlI3(s(i1,i4),s(i2,i3),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(83) = qlI3(s(i1,i4),s(i2,i5),s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(84) = qlI3(s(i1,i4),s(i2,i6),s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(85) = qlI3(s(i1,i5),s(i2,i3),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(86) = qlI3(s(i1,i5),s(i2,i4),s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(87) = qlI3(s(i1,i5),s(i2,i6),s(i3,i4),0d0,0d0,0d0,musq,ord)
      Triint(88) = qlI3(s(i1,i6),s(i2,i3),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(89) = qlI3(s(i1,i6),s(i2,i4),s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(90) = qlI3(s(i1,i6),s(i2,i5),s(i3,i4),0d0,0d0,0d0,musq,ord)

      return
      end
