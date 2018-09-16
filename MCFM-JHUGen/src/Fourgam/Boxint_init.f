!===== T. Dennen, May 2014
!===== Initialise array of box integral values.
!===== For use with qqb->4gamma and qqb->2j2gamma
!===== coefficients by same author.
      subroutine Boxint_init(i1,i2,i3,i4,i5,i6,Boxint,ord)
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
      complex(dp):: Boxint(195), qlI4
      real(dp):: t

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      
      Boxint(:) = czip

      Boxint(1) = 
     &  qlI4(0d0,0d0,0d0,t(i4,i5,i6),s(i1,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

!      write(6,*) Boxint(1),musq,t(i4,i5,i6)
 !     pause
      Boxint(2) = 
     &  qlI4(0d0,0d0,t(i3,i4,i5),0d0,s(i1,i2),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(3) = 
     &  qlI4(0d0,t(i2,i3,i4),0d0,0d0,s(i5,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(4) = 
     &  qlI4(t(i1,i2,i3),0d0,0d0,0d0,s(i5,i6),s(i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(5) = 
     &  qlI4(0d0,0d0,t(i3,i4,i6),0d0,s(i1,i2),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(6) = 
     &  qlI4(0d0,t(i2,i3,i4),0d0,0d0,s(i5,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(7) = 
     &  qlI4(t(i1,i2,i3),0d0,0d0,0d0,s(i5,i6),s(i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(8) = 
     &  qlI4(0d0,t(i2,i3,i5),0d0,0d0,s(i4,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(9) = 
     &  qlI4(t(i1,i2,i3),0d0,0d0,0d0,s(i4,i6),s(i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(10) = 
     &  qlI4(0d0,0d0,t(i3,i5,i6),0d0,s(i1,i2),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(11) = 
     &  qlI4(0d0,t(i2,i3,i5),0d0,0d0,s(i4,i6),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(12) = 
     &  qlI4(0d0,t(i2,i3,i6),0d0,0d0,s(i4,i5),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(13) = 
     &  qlI4(0d0,t(i2,i3,i6),0d0,0d0,s(i4,i5),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(14) = 
     &  qlI4(0d0,0d0,0d0,t(i3,i5,i6),s(i1,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(15) = 
     &  qlI4(t(i1,i2,i4),0d0,0d0,0d0,s(i5,i6),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(16) = 
     &  qlI4(t(i1,i2,i4),0d0,0d0,0d0,s(i5,i6),s(i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(17) = 
     &  qlI4(0d0,t(i2,i4,i5),0d0,0d0,s(i3,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(18) = 
     &  qlI4(t(i1,i2,i4),0d0,0d0,0d0,s(i3,i6),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(19) = 
     &  qlI4(0d0,0d0,t(i4,i5,i6),0d0,s(i1,i2),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(20) = 
     &  qlI4(0d0,t(i2,i4,i5),0d0,0d0,s(i3,i6),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(21) = 
     &  qlI4(0d0,t(i2,i4,i6),0d0,0d0,s(i3,i5),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(22) = 
     &  qlI4(0d0,t(i2,i4,i6),0d0,0d0,s(i3,i5),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(23) = 
     &  qlI4(0d0,0d0,0d0,t(i3,i4,i6),s(i1,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(24) = 
     &  qlI4(t(i1,i2,i5),0d0,0d0,0d0,s(i4,i6),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(25) = 
     &  qlI4(t(i1,i2,i5),0d0,0d0,0d0,s(i4,i6),s(i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(26) = 
     &  qlI4(t(i1,i2,i5),0d0,0d0,0d0,s(i3,i6),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(27) = 
     &  qlI4(0d0,t(i2,i5,i6),0d0,0d0,s(i3,i4),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(28) = 
     &  qlI4(0d0,t(i2,i5,i6),0d0,0d0,s(i3,i4),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(29) = 
     &  qlI4(0d0,0d0,0d0,t(i3,i4,i5),s(i1,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(30) = 
     &  qlI4(t(i1,i2,i6),0d0,0d0,0d0,s(i4,i5),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(31) = 
     &  qlI4(t(i1,i2,i6),0d0,0d0,0d0,s(i4,i5),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(32) = 
     &  qlI4(t(i1,i2,i6),0d0,0d0,0d0,s(i3,i5),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(33) = 
     &  qlI4(0d0,0d0,0d0,t(i4,i5,i6),s(i1,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(34) = 
     &  qlI4(0d0,0d0,t(i2,i4,i5),0d0,s(i1,i3),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(35) = 
     &  qlI4(0d0,0d0,t(i2,i4,i6),0d0,s(i1,i3),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(36) = 
     &  qlI4(0d0,0d0,t(i2,i5,i6),0d0,s(i1,i3),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(37) = 
     &  qlI4(t(i1,i3,i4),0d0,0d0,0d0,s(i5,i6),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(38) = 
     &  qlI4(t(i1,i3,i4),0d0,0d0,0d0,s(i5,i6),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(39) = 
     &  qlI4(0d0,t(i3,i4,i5),0d0,0d0,s(i2,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(40) = 
     &  qlI4(t(i1,i3,i4),0d0,0d0,0d0,s(i2,i6),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(41) = 
     &  qlI4(0d0,t(i3,i4,i6),0d0,0d0,s(i2,i5),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(42) = 
     &  qlI4(t(i1,i3,i5),0d0,0d0,0d0,s(i4,i6),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(43) = 
     &  qlI4(t(i1,i3,i5),0d0,0d0,0d0,s(i4,i6),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(44) = 
     &  qlI4(t(i1,i3,i5),0d0,0d0,0d0,s(i2,i6),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(45) = 
     &  qlI4(0d0,t(i3,i5,i6),0d0,0d0,s(i2,i4),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(46) = 
     &  qlI4(t(i1,i3,i6),0d0,0d0,0d0,s(i4,i5),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(47) = 
     &  qlI4(t(i1,i3,i6),0d0,0d0,0d0,s(i4,i5),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(48) = 
     &  qlI4(t(i1,i3,i6),0d0,0d0,0d0,s(i2,i5),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(49) = 
     &  qlI4(0d0,0d0,t(i2,i3,i5),0d0,s(i1,i4),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(50) = 
     &  qlI4(0d0,0d0,t(i2,i3,i6),0d0,s(i1,i4),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(51) = 
     &  qlI4(t(i1,i4,i5),0d0,0d0,0d0,s(i3,i6),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(52) = 
     &  qlI4(t(i1,i4,i5),0d0,0d0,0d0,s(i3,i6),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(53) = 
     &  qlI4(t(i1,i4,i5),0d0,0d0,0d0,s(i2,i6),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(54) = 
     &  qlI4(t(i1,i4,i6),0d0,0d0,0d0,s(i3,i5),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(55) = 
     &  qlI4(t(i1,i4,i6),0d0,0d0,0d0,s(i3,i5),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(56) = 
     &  qlI4(t(i1,i4,i6),0d0,0d0,0d0,s(i2,i5),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(57) = 
     &  qlI4(0d0,0d0,t(i2,i3,i4),0d0,s(i1,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(58) = 
     &  qlI4(t(i1,i5,i6),0d0,0d0,0d0,s(i3,i4),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(59) = 
     &  qlI4(t(i1,i5,i6),0d0,0d0,0d0,s(i3,i4),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(60) = 
     &  qlI4(t(i1,i5,i6),0d0,0d0,0d0,s(i2,i4),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(61) = 
     &  qlI4(0d0,s(i2,i3),0d0,s(i5,i6),t(i1,i2,i3),t(i2,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(62) = 
     &  qlI4(s(i1,i2),0d0,s(i4,i5),0d0,t(i1,i2,i3),t(i3,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(63) = 
     &  qlI4(s(i1,i2),0d0,s(i4,i6),0d0,t(i1,i2,i3),t(i3,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(64) = 
     &  qlI4(0d0,s(i2,i3),0d0,s(i4,i6),t(i1,i2,i3),t(i2,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(65) = 
     &  qlI4(s(i1,i2),0d0,s(i5,i6),0d0,t(i1,i2,i3),t(i3,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(66) = 
     &  qlI4(0d0,s(i2,i3),0d0,s(i4,i5),t(i1,i2,i3),t(i2,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(67) = 
     &  qlI4(0d0,s(i2,i4),0d0,s(i5,i6),t(i1,i2,i4),t(i2,i4,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(68) = 
     &  qlI4(s(i1,i2),0d0,s(i3,i5),0d0,t(i1,i2,i4),t(i4,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(69) = 
     &  qlI4(s(i1,i2),0d0,s(i3,i6),0d0,t(i1,i2,i4),t(i4,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(70) = 
     &  qlI4(0d0,s(i2,i4),0d0,s(i3,i6),t(i1,i2,i4),t(i2,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(71) = 
     &  qlI4(0d0,s(i2,i4),0d0,s(i3,i5),t(i1,i2,i4),t(i2,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(72) = 
     &  qlI4(0d0,s(i2,i5),0d0,s(i4,i6),t(i1,i2,i5),t(i2,i5,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(73) = 
     &  qlI4(s(i1,i2),0d0,s(i3,i4),0d0,t(i1,i2,i5),t(i5,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(74) = 
     &  qlI4(0d0,s(i2,i5),0d0,s(i3,i6),t(i1,i2,i5),t(i2,i5,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(75) = 
     &  qlI4(0d0,s(i2,i5),0d0,s(i3,i4),t(i1,i2,i5),t(i2,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(76) = 
     &  qlI4(0d0,s(i2,i6),0d0,s(i4,i5),t(i1,i2,i6),t(i2,i6,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(77) = 
     &  qlI4(0d0,s(i2,i6),0d0,s(i3,i5),t(i1,i2,i6),t(i2,i6,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(78) = 
     &  qlI4(0d0,s(i2,i6),0d0,s(i3,i4),t(i1,i2,i6),t(i2,i6,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(79) = 
     &  qlI4(s(i1,i3),0d0,s(i4,i5),0d0,t(i1,i3,i2),t(i2,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(80) = 
     &  qlI4(s(i1,i3),0d0,s(i4,i6),0d0,t(i1,i3,i2),t(i2,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(81) = 
     &  qlI4(s(i1,i3),0d0,s(i5,i6),0d0,t(i1,i3,i2),t(i2,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(82) = 
     &  qlI4(0d0,s(i3,i4),0d0,s(i5,i6),t(i1,i3,i4),t(i3,i4,i2),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(83) = 
     &  qlI4(s(i1,i3),0d0,s(i2,i5),0d0,t(i1,i3,i4),t(i4,i2,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(84) = 
     &  qlI4(s(i1,i3),0d0,s(i2,i6),0d0,t(i1,i3,i4),t(i4,i2,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(85) = 
     &  qlI4(0d0,s(i3,i5),0d0,s(i4,i6),t(i1,i3,i5),t(i3,i5,i2),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(86) = 
     &  qlI4(s(i1,i3),0d0,s(i2,i4),0d0,t(i1,i3,i5),t(i5,i2,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(87) = 
     &  qlI4(0d0,s(i3,i6),0d0,s(i4,i5),t(i1,i3,i6),t(i3,i6,i2),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(88) = 
     &  qlI4(s(i1,i4),0d0,s(i3,i5),0d0,t(i1,i4,i2),t(i2,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(89) = 
     &  qlI4(s(i1,i4),0d0,s(i3,i6),0d0,t(i1,i4,i2),t(i2,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(90) = 
     &  qlI4(s(i1,i4),0d0,s(i5,i6),0d0,t(i1,i4,i2),t(i2,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(91) = 
     &  qlI4(s(i1,i4),0d0,s(i2,i5),0d0,t(i1,i4,i3),t(i3,i2,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(92) = 
     &  qlI4(s(i1,i4),0d0,s(i2,i6),0d0,t(i1,i4,i3),t(i3,i2,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(93) = 
     &  qlI4(s(i1,i4),0d0,s(i2,i3),0d0,t(i1,i4,i5),t(i5,i2,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(94) = 
     &  qlI4(s(i1,i5),0d0,s(i3,i4),0d0,t(i1,i5,i2),t(i2,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(95) = 
     &  qlI4(s(i1,i5),0d0,s(i3,i6),0d0,t(i1,i5,i2),t(i2,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(96) = 
     &  qlI4(s(i1,i5),0d0,s(i4,i6),0d0,t(i1,i5,i2),t(i2,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(97) = 
     &  qlI4(s(i1,i5),0d0,s(i2,i4),0d0,t(i1,i5,i3),t(i3,i2,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(98) = 
     &  qlI4(s(i1,i5),0d0,s(i2,i6),0d0,t(i1,i5,i3),t(i3,i2,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(99) = 
     &  qlI4(s(i1,i5),0d0,s(i2,i3),0d0,t(i1,i5,i4),t(i4,i2,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(100) = 
     &  qlI4(s(i1,i6),0d0,s(i3,i4),0d0,t(i1,i6,i2),t(i2,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(101) = 
     &  qlI4(s(i1,i6),0d0,s(i3,i5),0d0,t(i1,i6,i2),t(i2,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(102) = 
     &  qlI4(s(i1,i6),0d0,s(i4,i5),0d0,t(i1,i6,i2),t(i2,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(103) = 
     &  qlI4(s(i1,i6),0d0,s(i2,i4),0d0,t(i1,i6,i3),t(i3,i2,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(104) = 
     &  qlI4(s(i1,i6),0d0,s(i2,i5),0d0,t(i1,i6,i3),t(i3,i2,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(105) = 
     &  qlI4(s(i1,i6),0d0,s(i2,i3),0d0,t(i1,i6,i4),t(i4,i2,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(106) = 
     &  qlI4(0d0,0d0,s(i3,i4),s(i5,i6),s(i1,i2),t(i2,i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(107) = 
     &  qlI4(0d0,s(i2,i3),s(i4,i5),0d0,t(i1,i2,i3),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(108) = 
     &  qlI4(s(i1,i2),0d0,0d0,s(i5,i6),t(i1,i2,i3),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(109) = 
     &  qlI4(s(i1,i2),s(i3,i4),0d0,0d0,s(i5,i6),t(i3,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(110) = 
     &  qlI4(0d0,s(i2,i3),s(i4,i6),0d0,t(i1,i2,i3),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(111) = 
     &  qlI4(s(i1,i2),s(i3,i4),0d0,0d0,s(i5,i6),t(i3,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(112) = 
     &  qlI4(0d0,0d0,s(i3,i5),s(i4,i6),s(i1,i2),t(i2,i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(113) = 
     &  qlI4(s(i1,i2),0d0,0d0,s(i4,i6),t(i1,i2,i3),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(114) = 
     &  qlI4(s(i1,i2),s(i3,i5),0d0,0d0,s(i4,i6),t(i3,i5,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(115) = 
     &  qlI4(0d0,s(i2,i3),s(i5,i6),0d0,t(i1,i2,i3),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(116) = 
     &  qlI4(s(i1,i2),s(i3,i5),0d0,0d0,s(i4,i6),t(i3,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(117) = 
     &  qlI4(0d0,0d0,s(i3,i6),s(i4,i5),s(i1,i2),t(i2,i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(118) = 
     &  qlI4(s(i1,i2),0d0,0d0,s(i4,i5),t(i1,i2,i3),s(i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(119) = 
     &  qlI4(s(i1,i2),s(i3,i6),0d0,0d0,s(i4,i5),t(i3,i6,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(120) = 
     &  qlI4(s(i1,i2),s(i3,i6),0d0,0d0,s(i4,i5),t(i3,i6,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(121) = 
     &  qlI4(0d0,s(i2,i4),s(i3,i5),0d0,t(i1,i2,i4),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(122) = 
     &  qlI4(s(i1,i2),0d0,0d0,s(i5,i6),t(i1,i2,i4),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(123) = 
     &  qlI4(0d0,s(i2,i4),s(i3,i6),0d0,t(i1,i2,i4),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(124) = 
     &  qlI4(0d0,0d0,s(i4,i5),s(i3,i6),s(i1,i2),t(i2,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(125) = 
     &  qlI4(s(i1,i2),s(i4,i5),0d0,0d0,s(i3,i6),t(i4,i5,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(126) = 
     &  qlI4(0d0,s(i2,i4),s(i5,i6),0d0,t(i1,i2,i4),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(127) = 
     &  qlI4(0d0,0d0,s(i4,i6),s(i3,i5),s(i1,i2),t(i2,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(128) = 
     &  qlI4(s(i1,i2),s(i4,i6),0d0,0d0,s(i3,i5),t(i4,i6,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(129) = 
     &  qlI4(0d0,s(i2,i5),s(i3,i4),0d0,t(i1,i2,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(130) = 
     &  qlI4(0d0,s(i2,i5),s(i3,i6),0d0,t(i1,i2,i5),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(131) = 
     &  qlI4(0d0,s(i2,i5),s(i4,i6),0d0,t(i1,i2,i5),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(132) = 
     &  qlI4(0d0,0d0,s(i5,i6),s(i3,i4),s(i1,i2),t(i2,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(133) = 
     &  qlI4(0d0,s(i2,i6),s(i3,i4),0d0,t(i1,i2,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(134) = 
     &  qlI4(0d0,s(i2,i6),s(i3,i5),0d0,t(i1,i2,i6),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(135) = 
     &  qlI4(0d0,s(i2,i6),s(i4,i5),0d0,t(i1,i2,i6),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(136) = 
     &  qlI4(0d0,0d0,s(i2,i4),s(i5,i6),s(i1,i3),t(i3,i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(137) = 
     &  qlI4(s(i1,i3),0d0,0d0,s(i5,i6),t(i1,i3,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(138) = 
     &  qlI4(s(i1,i3),s(i2,i4),0d0,0d0,s(i5,i6),t(i2,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(139) = 
     &  qlI4(s(i1,i3),s(i2,i4),0d0,0d0,s(i5,i6),t(i2,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(140) = 
     &  qlI4(0d0,0d0,s(i2,i5),s(i4,i6),s(i1,i3),t(i3,i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(141) = 
     &  qlI4(s(i1,i3),0d0,0d0,s(i4,i6),t(i1,i3,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(142) = 
     &  qlI4(s(i1,i3),s(i2,i5),0d0,0d0,s(i4,i6),t(i2,i5,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(143) = 
     &  qlI4(s(i1,i3),s(i2,i5),0d0,0d0,s(i4,i6),t(i2,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(144) = 
     &  qlI4(0d0,0d0,s(i2,i6),s(i4,i5),s(i1,i3),t(i3,i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(145) = 
     &  qlI4(s(i1,i3),0d0,0d0,s(i4,i5),t(i1,i3,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(146) = 
     &  qlI4(s(i1,i3),s(i2,i6),0d0,0d0,s(i4,i5),t(i2,i6,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(147) = 
     &  qlI4(s(i1,i3),s(i2,i6),0d0,0d0,s(i4,i5),t(i2,i6,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(148) = 
     &  qlI4(0d0,s(i3,i4),s(i2,i5),0d0,t(i1,i3,i4),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(149) = 
     &  qlI4(s(i1,i3),0d0,0d0,s(i5,i6),t(i1,i3,i4),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(150) = 
     &  qlI4(0d0,s(i3,i4),s(i2,i6),0d0,t(i1,i3,i4),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(151) = 
     &  qlI4(s(i1,i3),s(i4,i5),0d0,0d0,s(i2,i6),t(i4,i5,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(152) = 
     &  qlI4(s(i1,i3),s(i4,i6),0d0,0d0,s(i2,i5),t(i4,i6,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(153) = 
     &  qlI4(0d0,s(i3,i5),s(i2,i4),0d0,t(i1,i3,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(154) = 
     &  qlI4(0d0,s(i3,i5),s(i2,i6),0d0,t(i1,i3,i5),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(155) = 
     &  qlI4(0d0,s(i3,i6),s(i2,i4),0d0,t(i1,i3,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(156) = 
     &  qlI4(0d0,s(i3,i6),s(i2,i5),0d0,t(i1,i3,i6),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(157) = 
     &  qlI4(0d0,0d0,s(i2,i3),s(i5,i6),s(i1,i4),t(i4,i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(158) = 
     &  qlI4(s(i1,i4),0d0,0d0,s(i5,i6),t(i1,i4,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(159) = 
     &  qlI4(s(i1,i4),s(i2,i3),0d0,0d0,s(i5,i6),t(i2,i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(160) = 
     &  qlI4(s(i1,i4),s(i2,i3),0d0,0d0,s(i5,i6),t(i2,i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(161) = 
     &  qlI4(s(i1,i4),0d0,0d0,s(i3,i6),t(i1,i4,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(162) = 
     &  qlI4(s(i1,i4),s(i2,i5),0d0,0d0,s(i3,i6),t(i2,i5,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(163) = 
     &  qlI4(s(i1,i4),s(i2,i5),0d0,0d0,s(i3,i6),t(i2,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(164) = 
     &  qlI4(s(i1,i4),0d0,0d0,s(i3,i5),t(i1,i4,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(165) = 
     &  qlI4(s(i1,i4),s(i2,i6),0d0,0d0,s(i3,i5),t(i2,i6,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(166) = 
     &  qlI4(s(i1,i4),s(i2,i6),0d0,0d0,s(i3,i5),t(i2,i6,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(167) = 
     &  qlI4(s(i1,i4),0d0,0d0,s(i5,i6),t(i1,i4,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(168) = 
     &  qlI4(s(i1,i4),s(i3,i5),0d0,0d0,s(i2,i6),t(i3,i5,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(169) = 
     &  qlI4(s(i1,i4),s(i3,i6),0d0,0d0,s(i2,i5),t(i3,i6,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(170) = 
     &  qlI4(0d0,s(i4,i5),s(i2,i3),0d0,t(i1,i4,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(171) = 
     &  qlI4(0d0,s(i4,i6),s(i2,i3),0d0,t(i1,i4,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(172) = 
     &  qlI4(s(i1,i5),0d0,0d0,s(i4,i6),t(i1,i5,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(173) = 
     &  qlI4(s(i1,i5),s(i2,i3),0d0,0d0,s(i4,i6),t(i2,i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(174) = 
     &  qlI4(s(i1,i5),s(i2,i3),0d0,0d0,s(i4,i6),t(i2,i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(175) = 
     &  qlI4(s(i1,i5),0d0,0d0,s(i3,i6),t(i1,i5,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(176) = 
     &  qlI4(s(i1,i5),s(i2,i4),0d0,0d0,s(i3,i6),t(i2,i4,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(177) = 
     &  qlI4(s(i1,i5),s(i2,i4),0d0,0d0,s(i3,i6),t(i2,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(178) = 
     &  qlI4(s(i1,i5),0d0,0d0,s(i3,i4),t(i1,i5,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(179) = 
     &  qlI4(s(i1,i5),s(i2,i6),0d0,0d0,s(i3,i4),t(i2,i6,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(180) = 
     &  qlI4(s(i1,i5),s(i2,i6),0d0,0d0,s(i3,i4),t(i2,i6,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(181) = 
     &  qlI4(s(i1,i5),0d0,0d0,s(i4,i6),t(i1,i5,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(182) = 
     &  qlI4(s(i1,i5),s(i3,i4),0d0,0d0,s(i2,i6),t(i3,i4,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(183) = 
     &  qlI4(s(i1,i5),s(i3,i6),0d0,0d0,s(i2,i4),t(i3,i6,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(184) = 
     &  qlI4(s(i1,i6),0d0,0d0,s(i4,i5),t(i1,i6,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(185) = 
     &  qlI4(s(i1,i6),s(i2,i3),0d0,0d0,s(i4,i5),t(i2,i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(186) = 
     &  qlI4(s(i1,i6),s(i2,i3),0d0,0d0,s(i4,i5),t(i2,i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(187) = 
     &  qlI4(s(i1,i6),0d0,0d0,s(i3,i5),t(i1,i6,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(188) = 
     &  qlI4(s(i1,i6),s(i2,i4),0d0,0d0,s(i3,i5),t(i2,i4,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(189) = 
     &  qlI4(s(i1,i6),s(i2,i4),0d0,0d0,s(i3,i5),t(i2,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(190) = 
     &  qlI4(s(i1,i6),0d0,0d0,s(i3,i4),t(i1,i6,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(191) = 
     &  qlI4(s(i1,i6),s(i2,i5),0d0,0d0,s(i3,i4),t(i2,i5,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(192) = 
     &  qlI4(s(i1,i6),s(i2,i5),0d0,0d0,s(i3,i4),t(i2,i5,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(193) = 
     &  qlI4(s(i1,i6),0d0,0d0,s(i4,i5),t(i1,i6,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(194) = 
     &  qlI4(s(i1,i6),s(i3,i4),0d0,0d0,s(i2,i5),t(i3,i4,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(195) = 
     &  qlI4(s(i1,i6),s(i3,i5),0d0,0d0,s(i2,i4),t(i3,i5,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)


      return
      end
