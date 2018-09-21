MODULE ModHashCollection
use ModParameters
implicit none
save
public

! MCFM qq_VVqq hash - JHU conventions
integer, parameter :: Hash_MCFM_qqVVqq_Size = 164
integer, target :: Hash_MCFM_qqVVqq(1:Hash_MCFM_qqVVqq_Size,1:4)

! MCFM qq_VVqq hash (generation) - PDF/PDG conventions
integer, parameter :: Hash_MCFM_qqVVqq_Gen_Size = 121 ! Maximum size
integer, target :: Hash_MCFM_qqVVqq_Gen(1:Hash_MCFM_qqVVqq_Gen_Size,1:4)

! JHUGen onshell 0-jet hash - PDF conventions
integer, parameter :: Hash_PPXchannel_Size = 121 ! Maximum size
integer, target :: Hash_PPXchannel(1:Hash_PPXchannel_Size,1:3)

! JHUGen onshell VBF hash (approximate) - PDF/PDG conventions
integer, parameter :: Hash_OnshellVBF_Size = 121 ! Maximum size
integer, target :: Hash_OnshellVBF(1:Hash_OnshellVBF_Size,1:3)

! JHUGen onshell VBF hash (exact) - PDF/PDG conventions
integer, parameter :: Hash_OnshellVBF_nosplit_Size = 121 ! Maximum size
integer, target :: Hash_OnshellVBF_nosplit(1:Hash_OnshellVBF_nosplit_Size,1:3)

! JHUGen onshell HJJ hash (approximate) - PDF conventions
integer, parameter :: Hash_OnshellHJJ_Size = 121 ! Maximum size
integer, target :: Hash_OnshellHJJ(1:Hash_OnshellHJJ_Size,1:3)

! JHUGen onshell HJJ hash (exact) - PDF conventions
integer, parameter :: Hash_OnshellHJJ_nosplit_Size = 121 ! Maximum size
integer, target :: Hash_OnshellHJJ_nosplit(1:Hash_OnshellHJJ_nosplit_Size,1:3)

! JHUGen onshell GEN channel hash - PDF conventions
integer, parameter :: Hash_GENchannel_Size = 121 ! Maximum size
integer, target :: Hash_GENchannel(1:Hash_GENchannel_Size,1:3)

! JHUGen onshell TH hash - PDF conventions
integer, parameter :: Hash_THchannel_Size = 121 ! Maximum size
integer, target :: Hash_THchannel(1:Hash_THchannel_Size,1:3)

! Variable to control initialization
logical, private :: hashcoll_hashes_initialized = .false.

contains

subroutine Init_Hash_MCFM_qqVVqq()
implicit none
   Hash_MCFM_qqVVqq(1,:) = (/ Up_ , Chm_ , Up_ , Chm_ /)
   Hash_MCFM_qqVVqq(2,:) = (/ Dn_ , Str_ , Dn_ , Str_ /)
   Hash_MCFM_qqVVqq(3,:) = (/ Dn_ , Bot_ , Dn_ , Bot_ /)
   Hash_MCFM_qqVVqq(4,:) = (/ Str_ , Bot_ , Str_ , Bot_ /)
   Hash_MCFM_qqVVqq(5,:) = (/ Up_ , Str_ , Up_ , Str_ /)
   Hash_MCFM_qqVVqq(6,:) = (/ Up_ , Bot_ , Up_ , Bot_ /)
   Hash_MCFM_qqVVqq(7,:) = (/ Chm_ , Bot_ , Chm_ , Bot_ /)
   Hash_MCFM_qqVVqq(8,:) = (/ Dn_ , Chm_ , Dn_ , Chm_ /)
   Hash_MCFM_qqVVqq(9,:) = (/ Dn_ , Up_ , Dn_ , Up_ /)
   Hash_MCFM_qqVVqq(10,:) = (/ Str_ , Chm_ , Str_ , Chm_ /)
   Hash_MCFM_qqVVqq(11,:) = (/ Dn_ , Chm_ , Up_ , Str_ /)
   Hash_MCFM_qqVVqq(12,:) = (/ Up_ , Str_ , Dn_ , Chm_ /)
   Hash_MCFM_qqVVqq(13,:) = (/ Up_ , Up_ , Up_ , Up_ /)
   Hash_MCFM_qqVVqq(14,:) = (/ Chm_ , Chm_ , Chm_ , Chm_ /)
   Hash_MCFM_qqVVqq(15,:) = (/ Dn_ , Dn_ , Dn_ , Dn_ /)
   Hash_MCFM_qqVVqq(16,:) = (/ Str_ , Str_ , Str_ , Str_ /)
   Hash_MCFM_qqVVqq(17,:) = (/ Bot_ , Bot_ , Bot_ , Bot_ /)
   Hash_MCFM_qqVVqq(18,:) = (/ Chm_ , Up_ , Up_ , Chm_ /)
   Hash_MCFM_qqVVqq(19,:) = (/ Str_ , Dn_ , Dn_ , Str_ /)
   Hash_MCFM_qqVVqq(20,:) = (/ Bot_ , Dn_ , Dn_ , Bot_ /)
   Hash_MCFM_qqVVqq(21,:) = (/ Bot_ , Str_ , Str_ , Bot_ /)
   Hash_MCFM_qqVVqq(22,:) = (/ Str_ , Up_ , Up_ , Str_ /)
   Hash_MCFM_qqVVqq(23,:) = (/ Bot_ , Up_ , Up_ , Bot_ /)
   Hash_MCFM_qqVVqq(24,:) = (/ Bot_ , Chm_ , Chm_ , Bot_ /)

   Hash_MCFM_qqVVqq(25,:) = (/ Chm_ , Dn_ , Dn_ , Chm_ /)
   Hash_MCFM_qqVVqq(26,:) = (/ Up_ , Dn_ , Dn_ , Up_ /)
   Hash_MCFM_qqVVqq(27,:) = (/ Chm_ , Str_ , Str_ , Chm_ /)
   Hash_MCFM_qqVVqq(28,:) = (/ Chm_ , Dn_ , Up_ , Str_ /)
   Hash_MCFM_qqVVqq(29,:) = (/ Str_ , Up_ , Dn_ , Chm_ /)

   Hash_MCFM_qqVVqq(30,:) = (/ AChm_ , AUp_ , AChm_ , AUp_ /)
   Hash_MCFM_qqVVqq(31,:) = (/ AStr_ , ADn_ , AStr_ , ADn_ /)
   Hash_MCFM_qqVVqq(32,:) = (/ ABot_ , ADn_ , ABot_ , ADn_ /)
   Hash_MCFM_qqVVqq(33,:) = (/ ABot_ , AStr_ , ABot_ , AStr_ /)
   Hash_MCFM_qqVVqq(34,:) = (/ AStr_ , AUp_ , AStr_ , AUp_ /)
   Hash_MCFM_qqVVqq(35,:) = (/ ABot_ , AUp_ , ABot_ , AUp_ /)
   Hash_MCFM_qqVVqq(36,:) = (/ ABot_ , AChm_ , ABot_ , AChm_ /)
   Hash_MCFM_qqVVqq(37,:) = (/ AChm_ , ADn_ , AChm_ , ADn_ /)
   Hash_MCFM_qqVVqq(38,:) = (/ AUp_ , ADn_ , AUp_ , ADn_ /)
   Hash_MCFM_qqVVqq(39,:) = (/ AChm_ , AStr_ , AChm_ , AStr_ /)
   Hash_MCFM_qqVVqq(40,:) = (/ AStr_ , AUp_ , AChm_ , ADn_ /)
   Hash_MCFM_qqVVqq(41,:) = (/ AChm_ , ADn_ , AStr_ , AUp_ /)
   Hash_MCFM_qqVVqq(42,:) = (/ AUp_ , AUp_ , AUp_ , AUp_ /)
   Hash_MCFM_qqVVqq(43,:) = (/ AChm_ , AChm_ , AChm_ , AChm_ /)
   Hash_MCFM_qqVVqq(44,:) = (/ ADn_ , ADn_ , ADn_ , ADn_ /)
   Hash_MCFM_qqVVqq(45,:) = (/ AStr_ , AStr_ , AStr_ , AStr_ /)
   Hash_MCFM_qqVVqq(46,:) = (/ ABot_ , ABot_ , ABot_ , ABot_ /)

   Hash_MCFM_qqVVqq(47,:) = (/ AUp_ , AChm_ , AChm_ , AUp_ /)
   Hash_MCFM_qqVVqq(48,:) = (/ ADn_ , AStr_ , AStr_ , ADn_ /)
   Hash_MCFM_qqVVqq(49,:) = (/ ADn_ , ABot_ , ABot_ , ADn_ /)
   Hash_MCFM_qqVVqq(50,:) = (/ AStr_ , ABot_ , ABot_ , AStr_ /)
   Hash_MCFM_qqVVqq(51,:) = (/ AUp_ , AStr_ , AStr_ , AUp_ /)
   Hash_MCFM_qqVVqq(52,:) = (/ AUp_ , ABot_ , ABot_ , AUp_ /)
   Hash_MCFM_qqVVqq(53,:) = (/ AChm_ , ABot_ , ABot_ , AChm_ /)
   Hash_MCFM_qqVVqq(54,:) = (/ ADn_ , AChm_ , AChm_ , ADn_ /)
   Hash_MCFM_qqVVqq(55,:) = (/ ADn_ , AUp_ , AUp_ , ADn_ /)
   Hash_MCFM_qqVVqq(56,:) = (/ AStr_ , AChm_ , AChm_ , AStr_ /)
   Hash_MCFM_qqVVqq(57,:) = (/ AUp_ , AStr_ , AChm_ , ADn_ /)
   Hash_MCFM_qqVVqq(58,:) = (/ ADn_ , AChm_ , AStr_ , AUp_ /)

   Hash_MCFM_qqVVqq(59,:) = (/ AUp_ , Chm_ , AUp_ , Chm_ /)
   Hash_MCFM_qqVVqq(60,:) = (/ ADn_ , Str_ , ADn_ , Str_ /)
   Hash_MCFM_qqVVqq(61,:) = (/ ADn_ , Bot_ , ADn_ , Bot_ /)
   Hash_MCFM_qqVVqq(62,:) = (/ AStr_ , Bot_ , AStr_ , Bot_ /)
   Hash_MCFM_qqVVqq(63,:) = (/ AUp_ , Str_ , AUp_ , Str_ /)
   Hash_MCFM_qqVVqq(64,:) = (/ AUp_ , Bot_ , AUp_ , Bot_ /)
   Hash_MCFM_qqVVqq(65,:) = (/ AChm_ , Bot_ , AChm_ , Bot_ /)
   Hash_MCFM_qqVVqq(66,:) = (/ ADn_ , Chm_ , ADn_ , Chm_ /)
   Hash_MCFM_qqVVqq(67,:) = (/ ADn_ , Up_ , ADn_ , Up_ /)
   Hash_MCFM_qqVVqq(68,:) = (/ AStr_ , Chm_ , AStr_ , Chm_ /)
   Hash_MCFM_qqVVqq(69,:) = (/ AUp_ , Chm_ , ADn_ , Str_ /)
   Hash_MCFM_qqVVqq(70,:) = (/ ADn_ , Str_ , AUp_ , Chm_ /)
   Hash_MCFM_qqVVqq(71,:) = (/ AUp_ , Up_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(72,:) = (/ AChm_ , Chm_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(73,:) = (/ ADn_ , Dn_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(74,:) = (/ AStr_ , Str_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(75,:) = (/ ABot_ , Bot_ , ABot_ , Bot_ /)

   Hash_MCFM_qqVVqq(76,:) = (/ AChm_ , Up_ , AChm_ , Up_ /)
   Hash_MCFM_qqVVqq(77,:) = (/ AStr_ , Dn_ , AStr_ , Dn_ /)
   Hash_MCFM_qqVVqq(78,:) = (/ ABot_ , Dn_ , ABot_ , Dn_ /)
   Hash_MCFM_qqVVqq(79,:) = (/ ABot_ , Str_ , ABot_ , Str_ /)
   Hash_MCFM_qqVVqq(80,:) = (/ AStr_ , Up_ , AStr_ , Up_ /)
   Hash_MCFM_qqVVqq(81,:) = (/ ABot_ , Up_ , ABot_ , Up_ /)
   Hash_MCFM_qqVVqq(82,:) = (/ ABot_ , Chm_ , ABot_ , Chm_ /)
   Hash_MCFM_qqVVqq(83,:) = (/ AChm_ , Dn_ , AChm_ , Dn_ /)
   Hash_MCFM_qqVVqq(84,:) = (/ AUp_ , Dn_ , AUp_ , Dn_ /)
   Hash_MCFM_qqVVqq(85,:) = (/ AChm_ , Str_ , AChm_ , Str_ /)
   Hash_MCFM_qqVVqq(86,:) = (/ AStr_ , Dn_ , AChm_ , Up_ /)
   Hash_MCFM_qqVVqq(87,:) = (/ AChm_ , Up_ , AStr_ , Dn_ /)

   Hash_MCFM_qqVVqq(88,:) = (/ Chm_ , AUp_ , AUp_ , Chm_ /)
   Hash_MCFM_qqVVqq(89,:) = (/ Str_ , ADn_ , ADn_ , Str_ /)
   Hash_MCFM_qqVVqq(90,:) = (/ Bot_ , ADn_ , ADn_ , Bot_ /)
   Hash_MCFM_qqVVqq(91,:) = (/ Bot_ , AStr_ , AStr_ , Bot_ /)
   Hash_MCFM_qqVVqq(92,:) = (/ Str_ , AUp_ , AUp_ , Str_ /)
   Hash_MCFM_qqVVqq(93,:) = (/ Bot_ , AUp_ , AUp_ , Bot_ /)
   Hash_MCFM_qqVVqq(94,:) = (/ Bot_ , AChm_ , AChm_ , Bot_ /)
   Hash_MCFM_qqVVqq(95,:) = (/ Chm_ , ADn_ , ADn_ , Chm_ /)
   Hash_MCFM_qqVVqq(96,:) = (/ Up_ , ADn_ , ADn_ , Up_ /)
   Hash_MCFM_qqVVqq(97,:) = (/ Chm_ , AStr_ , AStr_ , Chm_ /)
   Hash_MCFM_qqVVqq(98,:) = (/ Chm_ , AUp_ , ADn_ , Str_ /)
   Hash_MCFM_qqVVqq(99,:) = (/ Str_ , ADn_ , AUp_ , Chm_ /)
   Hash_MCFM_qqVVqq(100,:) = (/ Up_ , AUp_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(101,:) = (/ Chm_ , AChm_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(102,:) = (/ Dn_ , ADn_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(103,:) = (/ Str_ , AStr_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(104,:) = (/ Bot_ , ABot_ , ABot_ , Bot_ /)

   Hash_MCFM_qqVVqq(105,:) = (/ Up_ , AChm_ , AChm_ , Up_ /)
   Hash_MCFM_qqVVqq(106,:) = (/ Dn_ , AStr_ , AStr_ , Dn_ /)
   Hash_MCFM_qqVVqq(107,:) = (/ Dn_ , ABot_ , ABot_ , Dn_ /)
   Hash_MCFM_qqVVqq(108,:) = (/ Str_ , ABot_ , ABot_ , Str_ /)
   Hash_MCFM_qqVVqq(109,:) = (/ Up_ , AStr_ , AStr_ , Up_ /)
   Hash_MCFM_qqVVqq(110,:) = (/ Up_ , ABot_ , ABot_ , Up_ /)
   Hash_MCFM_qqVVqq(111,:) = (/ Chm_ , ABot_ , ABot_ , Chm_ /)
   Hash_MCFM_qqVVqq(112,:) = (/ Dn_ , AChm_ , AChm_ , Dn_ /)
   Hash_MCFM_qqVVqq(113,:) = (/ Dn_ , AUp_ , AUp_ , Dn_ /)
   Hash_MCFM_qqVVqq(114,:) = (/ Str_ , AChm_ , AChm_ , Str_ /)
   Hash_MCFM_qqVVqq(115,:) = (/ Dn_ , AStr_ , AChm_ , Up_ /)
   Hash_MCFM_qqVVqq(116,:) = (/ Up_ , AChm_ , AStr_ , Dn_ /)
   Hash_MCFM_qqVVqq(117,:) = (/ Up_ , AUp_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(118,:) = (/ Dn_ , ADn_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(119,:) = (/ Dn_ , ADn_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(120,:) = (/ Str_ , AStr_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(121,:) = (/ Up_ , AUp_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(122,:) = (/ Up_ , AUp_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(123,:) = (/ Chm_ , AChm_ , ABot_ , Bot_ /)

   Hash_MCFM_qqVVqq(124,:) = (/ Dn_ , ADn_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(125,:) = (/ Dn_ , ADn_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(126,:) = (/ Str_ , AStr_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(127,:) = (/ Dn_ , AUp_ , AChm_ , Str_ /)
   Hash_MCFM_qqVVqq(128,:) = (/ Up_ , ADn_ , AStr_ , Chm_ /)
   Hash_MCFM_qqVVqq(129,:) = (/ Chm_ , AChm_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(130,:) = (/ Str_ , AStr_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(131,:) = (/ Bot_ , ABot_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(132,:) = (/ Bot_ , ABot_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(133,:) = (/ Str_ , AStr_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(134,:) = (/ Bot_ , ABot_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(135,:) = (/ Bot_ , ABot_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(136,:) = (/ Chm_ , AChm_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(137,:) = (/ Up_ , AUp_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(138,:) = (/ Chm_ , AChm_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(139,:) = (/ Chm_ , AStr_ , ADn_ , Up_ /)
   Hash_MCFM_qqVVqq(140,:) = (/ Str_ , AChm_ , AUp_ , Dn_ /)

   Hash_MCFM_qqVVqq(141,:) = (/ AUp_ , Up_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(142,:) = (/ ADn_ , Dn_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(143,:) = (/ ADn_ , Dn_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(144,:) = (/ AStr_ , Str_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(145,:) = (/ AUp_ , Up_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(146,:) = (/ AUp_ , Up_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(147,:) = (/ AChm_ , Chm_ , ABot_ , Bot_ /)
   Hash_MCFM_qqVVqq(148,:) = (/ ADn_ , Dn_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(149,:) = (/ ADn_ , Dn_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(150,:) = (/ AStr_ , Str_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(151,:) = (/ AUp_ , Dn_ , AChm_ , Str_ /)
   Hash_MCFM_qqVVqq(152,:) = (/ ADn_ , Up_ , AStr_ , Chm_ /)
   Hash_MCFM_qqVVqq(153,:) = (/ AChm_ , Chm_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(154,:) = (/ AStr_ , Str_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(155,:) = (/ ABot_ , Bot_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(156,:) = (/ ABot_ , Bot_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(157,:) = (/ AStr_ , Str_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(158,:) = (/ ABot_ , Bot_ , AUp_ , Up_ /)
   Hash_MCFM_qqVVqq(159,:) = (/ ABot_ , Bot_ , AChm_ , Chm_ /)
   Hash_MCFM_qqVVqq(160,:) = (/ AChm_ , Chm_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(161,:) = (/ AUp_ , Up_ , ADn_ , Dn_ /)
   Hash_MCFM_qqVVqq(162,:) = (/ AChm_ , Chm_ , AStr_ , Str_ /)
   Hash_MCFM_qqVVqq(163,:) = (/ AStr_ , Chm_ , ADn_ , Up_ /)
   Hash_MCFM_qqVVqq(164,:) = (/ AChm_ , Str_ , AUp_ , Dn_ /)
end subroutine

subroutine Init_Hash_MCFM_qqVVqq_Generation()
implicit none
integer :: ch,i,j
   Hash_MCFM_qqVVqq_Gen=Not_a_particle_
   do i=1,5; do j=1,5
      Hash_MCFM_qqVVqq_Gen(20*i+4*j-23,:)=(/-i,-j,0,0/)
      Hash_MCFM_qqVVqq_Gen(20*i+4*j-22,:)=(/-i,j,0,0/)
      Hash_MCFM_qqVVqq_Gen(20*i+4*j-21,:)=(/i,-j,0,0/)
      Hash_MCFM_qqVVqq_Gen(20*i+4*j-20,:)=(/i,j,0,0/)
   enddo; enddo
   return
end subroutine

subroutine Init_Hash_OnshellVBF()
implicit none
integer,parameter :: zz=1, ww=0
   Hash_OnshellVBF(  1,1:3) = (/ 2, 1, zz/)
   Hash_OnshellVBF(  2,1:3) = (/ 2, 1, ww/)
   Hash_OnshellVBF(  3,1:3) = (/ 2,-2, zz/)
   Hash_OnshellVBF(  4,1:3) = (/ 2,-2, ww/)
   Hash_OnshellVBF(  5,1:3) = (/ 3, 2, zz/)
   Hash_OnshellVBF(  6,1:3) = (/ 3, 2, ww/)
   Hash_OnshellVBF(  7,1:3) = (/ 1,-1, zz/)
   Hash_OnshellVBF(  8,1:3) = (/ 1,-1, ww/)
   Hash_OnshellVBF(  9,1:3) = (/ 2,-4, zz/)
   Hash_OnshellVBF( 10,1:3) = (/ 2,-4, ww/)
   Hash_OnshellVBF( 11,1:3) = (/ 1,-3, zz/)
   Hash_OnshellVBF( 12,1:3) = (/ 1,-3, ww/)
   Hash_OnshellVBF( 13,1:3) = (/ 4, 1, zz/)
   Hash_OnshellVBF( 14,1:3) = (/ 4, 1, ww/)
   Hash_OnshellVBF( 15,1:3) = (/ 2, 2, zz/)
   Hash_OnshellVBF( 16,1:3) = (/-1,-2, zz/)
   Hash_OnshellVBF( 17,1:3) = (/-1,-2, ww/)
   Hash_OnshellVBF( 18,1:3) = (/ 2,-1, zz/)
   Hash_OnshellVBF( 19,1:3) = (/ 1, 1, zz/)
   Hash_OnshellVBF( 20,1:3) = (/ 3,-1, zz/)
   Hash_OnshellVBF( 21,1:3) = (/ 3,-1, ww/)
   Hash_OnshellVBF( 22,1:3) = (/ 2,-3, zz/)
   Hash_OnshellVBF( 23,1:3) = (/-2,-3, zz/)
   Hash_OnshellVBF( 24,1:3) = (/-2,-3, ww/)
   Hash_OnshellVBF( 25,1:3) = (/-1,-4, zz/)
   Hash_OnshellVBF( 26,1:3) = (/-1,-4, ww/)
   Hash_OnshellVBF( 27,1:3) = (/ 3,-3, zz/)
   Hash_OnshellVBF( 28,1:3) = (/ 3,-3, ww/)
   Hash_OnshellVBF( 29,1:3) = (/ 1,-2, zz/)
   Hash_OnshellVBF( 30,1:3) = (/ 3, 1, zz/)
   Hash_OnshellVBF( 31,1:3) = (/ 4,-2, zz/)
   Hash_OnshellVBF( 32,1:3) = (/ 4,-2, ww/)
   Hash_OnshellVBF( 33,1:3) = (/ 4, 2, zz/)
   Hash_OnshellVBF( 34,1:3) = (/ 5, 2, zz/)
   Hash_OnshellVBF( 35,1:3) = (/ 2,-5, zz/)
   Hash_OnshellVBF( 36,1:3) = (/-3,-4, zz/)
   Hash_OnshellVBF( 37,1:3) = (/-3,-4, ww/)
   Hash_OnshellVBF( 38,1:3) = (/ 4, 3, zz/)
   Hash_OnshellVBF( 39,1:3) = (/ 4, 3, ww/)
   Hash_OnshellVBF( 40,1:3) = (/ 1,-4, zz/)
   Hash_OnshellVBF( 41,1:3) = (/ 5, 1, zz/)
   Hash_OnshellVBF( 42,1:3) = (/ 1,-5, zz/)
   Hash_OnshellVBF( 43,1:3) = (/ 4,-4, zz/)
   Hash_OnshellVBF( 44,1:3) = (/ 4,-4, ww/)
   Hash_OnshellVBF( 45,1:3) = (/-1,-3, zz/)
   Hash_OnshellVBF( 46,1:3) = (/-1,-1, zz/)
   Hash_OnshellVBF( 47,1:3) = (/ 3,-2, zz/)
   Hash_OnshellVBF( 48,1:3) = (/ 4,-1, zz/)
   Hash_OnshellVBF( 49,1:3) = (/-1,-5, zz/)
   Hash_OnshellVBF( 50,1:3) = (/ 5,-1, zz/)
   Hash_OnshellVBF( 51,1:3) = (/-2,-2, zz/)
   Hash_OnshellVBF( 52,1:3) = (/-2,-4, zz/)
   Hash_OnshellVBF( 53,1:3) = (/ 3, 3, zz/)
   Hash_OnshellVBF( 54,1:3) = (/-3,-3, zz/)
   Hash_OnshellVBF( 55,1:3) = (/ 3,-4, zz/)
   Hash_OnshellVBF( 56,1:3) = (/ 4,-3, zz/)
   Hash_OnshellVBF( 57,1:3) = (/ 5,-2, zz/)
   Hash_OnshellVBF( 58,1:3) = (/-2,-5, zz/)
   Hash_OnshellVBF( 59,1:3) = (/ 5, 3, zz/)
   Hash_OnshellVBF( 60,1:3) = (/ 3,-5, zz/)
   Hash_OnshellVBF( 61,1:3) = (/-3,-5, zz/)
   Hash_OnshellVBF( 62,1:3) = (/ 5,-3, zz/)
   Hash_OnshellVBF( 63,1:3) = (/-4,-5, zz/)
   Hash_OnshellVBF( 64,1:3) = (/ 5, 4, zz/)
   Hash_OnshellVBF( 65,1:3) = (/ 4,-5, zz/)
   Hash_OnshellVBF( 66,1:3) = (/ 5,-4, zz/)
   Hash_OnshellVBF( 67,1:3) = (/ 5,-5, zz/)
   Hash_OnshellVBF( 68,1:3) = (/ 4, 4, zz/)
   Hash_OnshellVBF( 69,1:3) = (/-4,-4, zz/)
   Hash_OnshellVBF( 70,1:3) = (/-5,-5, zz/)
   Hash_OnshellVBF( 71,1:3) = (/ 5, 5, zz/)

   Hash_OnshellVBF( 72:,:)  = 0
   Hash_OnshellVBF( 72:,3)  = -1
end subroutine

subroutine Init_Hash_OnshellVBF_nosplit()
implicit none
integer, parameter :: nijchannels = 68
integer, parameter :: zzww=2, zz=1, ww=0
   ! u-d
   Hash_OnshellVBF_nosplit(  1,1:3) = (/ pdfUp_,pdfDn_,zzww/)
   Hash_OnshellVBF_nosplit(  2,1:3) = (/ pdfUp_,pdfStr_,zzww/)
   Hash_OnshellVBF_nosplit(  3,1:3) = (/ pdfUp_,pdfBot_,zzww/)
   Hash_OnshellVBF_nosplit(  4,1:3) = (/ pdfChm_,pdfDn_,zzww/)
   Hash_OnshellVBF_nosplit(  5,1:3) = (/ pdfChm_,pdfStr_,zzww/)
   Hash_OnshellVBF_nosplit(  6,1:3) = (/ pdfChm_,pdfBot_,zzww/)

   ! ub-db
   Hash_OnshellVBF_nosplit(  7,1:3) = (/ pdfAUp_,pdfADn_,zzww/)
   Hash_OnshellVBF_nosplit(  8,1:3) = (/ pdfAUp_,pdfAStr_,zzww/)
   Hash_OnshellVBF_nosplit(  9,1:3) = (/ pdfAUp_,pdfABot_,zzww/)
   Hash_OnshellVBF_nosplit( 10,1:3) = (/ pdfAChm_,pdfADn_,zzww/)
   Hash_OnshellVBF_nosplit( 11,1:3) = (/ pdfAChm_,pdfAStr_,zzww/)
   Hash_OnshellVBF_nosplit( 12,1:3) = (/ pdfAChm_,pdfABot_,zzww/)

   ! u-ub
   Hash_OnshellVBF_nosplit( 13,1:3) = (/ pdfUp_,pdfAUp_,ww/)
   Hash_OnshellVBF_nosplit( 14,1:3) = (/ pdfUp_,pdfAUp_,zz/)
   Hash_OnshellVBF_nosplit( 15,1:3) = (/ pdfUp_,pdfAChm_,ww/)
   Hash_OnshellVBF_nosplit( 16,1:3) = (/ pdfUp_,pdfAChm_,zz/)
   Hash_OnshellVBF_nosplit( 17,1:3) = (/ pdfChm_,pdfAUp_,ww/)
   Hash_OnshellVBF_nosplit( 18,1:3) = (/ pdfChm_,pdfAUp_,zz/)
   Hash_OnshellVBF_nosplit( 19,1:3) = (/ pdfChm_,pdfAChm_,ww/)
   Hash_OnshellVBF_nosplit( 20,1:3) = (/ pdfChm_,pdfAChm_,zz/)

   ! d-db
   Hash_OnshellVBF_nosplit( 21,1:3) = (/ pdfDn_,pdfADn_,ww/)
   Hash_OnshellVBF_nosplit( 22,1:3) = (/ pdfDn_,pdfADn_,zz/)
   Hash_OnshellVBF_nosplit( 23,1:3) = (/ pdfDn_,pdfAStr_,ww/)
   Hash_OnshellVBF_nosplit( 24,1:3) = (/ pdfDn_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 25,1:3) = (/ pdfDn_,pdfABot_,ww/)
   Hash_OnshellVBF_nosplit( 26,1:3) = (/ pdfDn_,pdfABot_,zz/)
   Hash_OnshellVBF_nosplit( 27,1:3) = (/ pdfStr_,pdfADn_,ww/)
   Hash_OnshellVBF_nosplit( 28,1:3) = (/ pdfStr_,pdfADn_,zz/)
   Hash_OnshellVBF_nosplit( 29,1:3) = (/ pdfStr_,pdfAStr_,ww/)
   Hash_OnshellVBF_nosplit( 30,1:3) = (/ pdfStr_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 31,1:3) = (/ pdfStr_,pdfABot_,ww/)
   Hash_OnshellVBF_nosplit( 32,1:3) = (/ pdfStr_,pdfABot_,zz/)

   ! b-db
   ! Some of these will be inefficient
   Hash_OnshellVBF_nosplit( 33,1:3) = (/ pdfBot_,pdfADn_,ww/)
   Hash_OnshellVBF_nosplit( 34,1:3) = (/ pdfBot_,pdfADn_,zz/)
   Hash_OnshellVBF_nosplit( 35,1:3) = (/ pdfBot_,pdfAStr_,ww/)
   Hash_OnshellVBF_nosplit( 36,1:3) = (/ pdfBot_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 37,1:3) = (/ pdfBot_,pdfABot_,ww/)
   Hash_OnshellVBF_nosplit( 38,1:3) = (/ pdfBot_,pdfABot_,zz/)

   ! u-u
   Hash_OnshellVBF_nosplit( 39,1:3) = (/ pdfUp_,pdfUp_,zz/)
   Hash_OnshellVBF_nosplit( 40,1:3) = (/ pdfUp_,pdfChm_,zz/)
   Hash_OnshellVBF_nosplit( 41,1:3) = (/ pdfChm_,pdfChm_,zz/)

   ! d-d
   Hash_OnshellVBF_nosplit( 42,1:3) = (/ pdfDn_,pdfDn_,zz/)
   Hash_OnshellVBF_nosplit( 43,1:3) = (/ pdfDn_,pdfStr_,zz/)
   Hash_OnshellVBF_nosplit( 44,1:3) = (/ pdfDn_,pdfBot_,zz/)
   Hash_OnshellVBF_nosplit( 45,1:3) = (/ pdfStr_,pdfStr_,zz/)
   Hash_OnshellVBF_nosplit( 46,1:3) = (/ pdfStr_,pdfBot_,zz/)
   Hash_OnshellVBF_nosplit( 47,1:3) = (/ pdfBot_,pdfBot_,zz/)

   ! u-db
   Hash_OnshellVBF_nosplit( 48,1:3) = (/ pdfUp_,pdfADn_,zz/)
   Hash_OnshellVBF_nosplit( 49,1:3) = (/ pdfUp_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 50,1:3) = (/ pdfUp_,pdfABot_,zz/)
   Hash_OnshellVBF_nosplit( 51,1:3) = (/ pdfChm_,pdfADn_,zz/)
   Hash_OnshellVBF_nosplit( 52,1:3) = (/ pdfChm_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 53,1:3) = (/ pdfChm_,pdfABot_,zz/)

   ! d-ub
   Hash_OnshellVBF_nosplit( 54,1:3) = (/ pdfDn_,pdfAUp_,zz/)
   Hash_OnshellVBF_nosplit( 55,1:3) = (/ pdfDn_,pdfAChm_,zz/)
   Hash_OnshellVBF_nosplit( 56,1:3) = (/ pdfStr_,pdfAUp_,zz/)
   Hash_OnshellVBF_nosplit( 57,1:3) = (/ pdfStr_,pdfAChm_,zz/)
   Hash_OnshellVBF_nosplit( 58,1:3) = (/ pdfBot_,pdfAUp_,zz/)
   Hash_OnshellVBF_nosplit( 59,1:3) = (/ pdfBot_,pdfAChm_,zz/)

   ! ub-ub
   Hash_OnshellVBF_nosplit( 60,1:3) = (/ pdfAUp_,pdfAUp_,zz/)
   Hash_OnshellVBF_nosplit( 61,1:3) = (/ pdfAUp_,pdfAChm_,zz/)
   Hash_OnshellVBF_nosplit( 62,1:3) = (/ pdfAChm_,pdfAChm_,zz/)

   ! db-db
   Hash_OnshellVBF_nosplit( 63,1:3) = (/ pdfADn_,pdfADn_,zz/)
   Hash_OnshellVBF_nosplit( 64,1:3) = (/ pdfADn_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 65,1:3) = (/ pdfADn_,pdfABot_,zz/)

   Hash_OnshellVBF_nosplit( 66,1:3) = (/ pdfAStr_,pdfAStr_,zz/)
   Hash_OnshellVBF_nosplit( 67,1:3) = (/ pdfAStr_,pdfABot_,zz/)
   Hash_OnshellVBF_nosplit( 68,1:3) = (/ pdfABot_,pdfABot_,zz/)

   !
   Hash_OnshellVBF_nosplit( nijchannels+1:,:)  = 0
   Hash_OnshellVBF_nosplit( nijchannels+1:,3)  = -1
return
end subroutine

subroutine Init_Hash_OnshellHJJ()
implicit none
   Hash_OnshellHJJ(  1,1:3) = (/ 0, 0, 1/)
   Hash_OnshellHJJ(  2,1:3) = (/ 0, 0, 2/)
   Hash_OnshellHJJ(  3,1:3) = (/ 2, 0, 1/)
   Hash_OnshellHJJ(  4,1:3) = (/ 1, 0, 1/)
   Hash_OnshellHJJ(  5,1:3) = (/ 0,-1, 1/)
   Hash_OnshellHJJ(  6,1:3) = (/ 0,-2, 1/)
   Hash_OnshellHJJ(  7,1:3) = (/ 3, 0, 1/)
   Hash_OnshellHJJ(  8,1:3) = (/ 0,-3, 1/)
   Hash_OnshellHJJ(  9,1:3) = (/ 4, 0, 1/)
   Hash_OnshellHJJ( 10,1:3) = (/ 0,-4, 1/)
   Hash_OnshellHJJ( 11,1:3) = (/ 5, 0, 1/)
   Hash_OnshellHJJ( 12,1:3) = (/ 0,-5, 1/)
   Hash_OnshellHJJ( 13,1:3) = (/ 2, 2, 1/)
   Hash_OnshellHJJ( 14,1:3) = (/ 2, 1, 1/)
   Hash_OnshellHJJ( 15,1:3) = (/ 2,-2, 1/)
   Hash_OnshellHJJ( 16,1:3) = (/ 2,-2, 2/)
   Hash_OnshellHJJ( 17,1:3) = (/ 2,-2, 3/)
   Hash_OnshellHJJ( 18,1:3) = (/ 1,-1, 1/)
   Hash_OnshellHJJ( 19,1:3) = (/ 1,-1, 2/)
   Hash_OnshellHJJ( 20,1:3) = (/ 1,-1, 3/)
   Hash_OnshellHJJ( 21,1:3) = (/ 2,-1, 1/)
   Hash_OnshellHJJ( 22,1:3) = (/ 1, 1, 1/)
   Hash_OnshellHJJ( 23,1:3) = (/ 3, 2, 1/)
   Hash_OnshellHJJ( 24,1:3) = (/ 2,-3, 1/)
   Hash_OnshellHJJ( 25,1:3) = (/ 4, 2, 1/)
   Hash_OnshellHJJ( 26,1:3) = (/ 1,-2, 1/)
   Hash_OnshellHJJ( 27,1:3) = (/ 2,-4, 1/)
   Hash_OnshellHJJ( 28,1:3) = (/ 1,-3, 1/)
   Hash_OnshellHJJ( 29,1:3) = (/ 2,-5, 1/)
   Hash_OnshellHJJ( 30,1:3) = (/ 3, 1, 1/)
   Hash_OnshellHJJ( 31,1:3) = (/ 5, 2, 1/)
   Hash_OnshellHJJ( 32,1:3) = (/ 4, 1, 1/)
   Hash_OnshellHJJ( 33,1:3) = (/ 1,-4, 1/)
   Hash_OnshellHJJ( 34,1:3) = (/ 5, 1, 1/)
   Hash_OnshellHJJ( 35,1:3) = (/ 1,-5, 1/)
   Hash_OnshellHJJ( 36,1:3) = (/-1,-2, 1/)
   Hash_OnshellHJJ( 37,1:3) = (/ 3,-1, 1/)
   Hash_OnshellHJJ( 38,1:3) = (/ 3,-3, 1/)
   Hash_OnshellHJJ( 39,1:3) = (/ 3,-3, 2/)
   Hash_OnshellHJJ( 40,1:3) = (/ 3,-3, 3/)
   Hash_OnshellHJJ( 41,1:3) = (/-2,-3, 1/)
   Hash_OnshellHJJ( 42,1:3) = (/-1,-3, 1/)
   Hash_OnshellHJJ( 43,1:3) = (/ 3,-2, 1/)
   Hash_OnshellHJJ( 44,1:3) = (/-1,-1, 1/)
   Hash_OnshellHJJ( 45,1:3) = (/ 4,-1, 1/)
   Hash_OnshellHJJ( 46,1:3) = (/-2,-2, 1/)
   Hash_OnshellHJJ( 47,1:3) = (/-1,-4, 1/)
   Hash_OnshellHJJ( 48,1:3) = (/ 4,-2, 1/)
   Hash_OnshellHJJ( 49,1:3) = (/-2,-4, 1/)
   Hash_OnshellHJJ( 50,1:3) = (/ 5,-1, 1/)
   Hash_OnshellHJJ( 51,1:3) = (/-1,-5, 1/)
   Hash_OnshellHJJ( 52,1:3) = (/ 4,-4, 1/)
   Hash_OnshellHJJ( 53,1:3) = (/ 4,-4, 2/)
   Hash_OnshellHJJ( 54,1:3) = (/ 4,-4, 3/)
   Hash_OnshellHJJ( 55,1:3) = (/ 4,-3, 1/)
   Hash_OnshellHJJ( 56,1:3) = (/ 4, 3, 1/)
   Hash_OnshellHJJ( 57,1:3) = (/ 3,-4, 1/)
   Hash_OnshellHJJ( 58,1:3) = (/-2,-5, 1/)
   Hash_OnshellHJJ( 59,1:3) = (/-3,-4, 1/)
   Hash_OnshellHJJ( 60,1:3) = (/ 5,-2, 1/)
   Hash_OnshellHJJ( 61,1:3) = (/ 3, 3, 1/)
   Hash_OnshellHJJ( 62,1:3) = (/-3,-3, 1/)
   Hash_OnshellHJJ( 63,1:3) = (/ 5, 3, 1/)
   Hash_OnshellHJJ( 64,1:3) = (/ 3,-5, 1/)
   Hash_OnshellHJJ( 65,1:3) = (/ 5,-3, 1/)
   Hash_OnshellHJJ( 66,1:3) = (/-3,-5, 1/)
   Hash_OnshellHJJ( 67,1:3) = (/ 5, 4, 1/)
   Hash_OnshellHJJ( 68,1:3) = (/ 5,-4, 1/)
   Hash_OnshellHJJ( 69,1:3) = (/-4,-5, 1/)
   Hash_OnshellHJJ( 70,1:3) = (/ 4,-5, 1/)
   Hash_OnshellHJJ( 71,1:3) = (/ 5,-5, 1/)
   Hash_OnshellHJJ( 72,1:3) = (/ 5,-5, 2/)
   Hash_OnshellHJJ( 73,1:3) = (/ 5,-5, 3/)
   Hash_OnshellHJJ( 74,1:3) = (/-4,-4, 1/)
   Hash_OnshellHJJ( 75,1:3) = (/ 4, 4, 1/)
   Hash_OnshellHJJ( 76,1:3) = (/-5,-5, 1/)
   Hash_OnshellHJJ( 77,1:3) = (/ 5, 5, 1/)

   Hash_OnshellHJJ( 78:,:) = 0
   Hash_OnshellHJJ( 78:,3) = -1
return
end subroutine

subroutine Init_Hash_OnshellHJJ_nosplit()
implicit none
integer, parameter :: nijchannels = 77
   Hash_OnshellHJJ_nosplit(  1,1:3) = (/ 0, 0, 1/)
   Hash_OnshellHJJ_nosplit(  2,1:3) = (/ 0, 0, 2/)
   Hash_OnshellHJJ_nosplit(  3,1:3) = (/ 2, 0, 1/)
   Hash_OnshellHJJ_nosplit(  4,1:3) = (/ 1, 0, 1/)
   Hash_OnshellHJJ_nosplit(  5,1:3) = (/ 0,-1, 1/)
   Hash_OnshellHJJ_nosplit(  6,1:3) = (/ 0,-2, 1/)
   Hash_OnshellHJJ_nosplit(  7,1:3) = (/ 3, 0, 1/)
   Hash_OnshellHJJ_nosplit(  8,1:3) = (/ 0,-3, 1/)
   Hash_OnshellHJJ_nosplit(  9,1:3) = (/ 4, 0, 1/)
   Hash_OnshellHJJ_nosplit( 10,1:3) = (/ 0,-4, 1/)
   Hash_OnshellHJJ_nosplit( 11,1:3) = (/ 5, 0, 1/)
   Hash_OnshellHJJ_nosplit( 12,1:3) = (/ 0,-5, 1/)
   Hash_OnshellHJJ_nosplit( 13,1:3) = (/ 2, 2, 1/)
   Hash_OnshellHJJ_nosplit( 14,1:3) = (/ 2, 1, 1/)
   Hash_OnshellHJJ_nosplit( 15,1:3) = (/ 2,-2, 1/)
   Hash_OnshellHJJ_nosplit( 16,1:3) = (/ 2,-2, 2/)
   Hash_OnshellHJJ_nosplit( 17,1:3) = (/ 2,-2, 3/)
   Hash_OnshellHJJ_nosplit( 18,1:3) = (/ 1,-1, 1/)
   Hash_OnshellHJJ_nosplit( 19,1:3) = (/ 1,-1, 2/)
   Hash_OnshellHJJ_nosplit( 20,1:3) = (/ 1,-1, 3/)
   Hash_OnshellHJJ_nosplit( 21,1:3) = (/ 2,-1, 1/)
   Hash_OnshellHJJ_nosplit( 22,1:3) = (/ 1, 1, 1/)
   Hash_OnshellHJJ_nosplit( 23,1:3) = (/ 3, 2, 1/)
   Hash_OnshellHJJ_nosplit( 24,1:3) = (/ 2,-3, 1/)
   Hash_OnshellHJJ_nosplit( 25,1:3) = (/ 4, 2, 1/)
   Hash_OnshellHJJ_nosplit( 26,1:3) = (/ 1,-2, 1/)
   Hash_OnshellHJJ_nosplit( 27,1:3) = (/ 2,-4, 1/)
   Hash_OnshellHJJ_nosplit( 28,1:3) = (/ 1,-3, 1/)
   Hash_OnshellHJJ_nosplit( 29,1:3) = (/ 2,-5, 1/)
   Hash_OnshellHJJ_nosplit( 30,1:3) = (/ 3, 1, 1/)
   Hash_OnshellHJJ_nosplit( 31,1:3) = (/ 5, 2, 1/)
   Hash_OnshellHJJ_nosplit( 32,1:3) = (/ 4, 1, 1/)
   Hash_OnshellHJJ_nosplit( 33,1:3) = (/ 1,-4, 1/)
   Hash_OnshellHJJ_nosplit( 34,1:3) = (/ 5, 1, 1/)
   Hash_OnshellHJJ_nosplit( 35,1:3) = (/ 1,-5, 1/)
   Hash_OnshellHJJ_nosplit( 36,1:3) = (/-1,-2, 1/)
   Hash_OnshellHJJ_nosplit( 37,1:3) = (/ 3,-1, 1/)
   Hash_OnshellHJJ_nosplit( 38,1:3) = (/ 3,-3, 1/)
   Hash_OnshellHJJ_nosplit( 39,1:3) = (/ 3,-3, 2/)
   Hash_OnshellHJJ_nosplit( 40,1:3) = (/ 3,-3, 3/)
   Hash_OnshellHJJ_nosplit( 41,1:3) = (/-2,-3, 1/)
   Hash_OnshellHJJ_nosplit( 42,1:3) = (/-1,-3, 1/)
   Hash_OnshellHJJ_nosplit( 43,1:3) = (/ 3,-2, 1/)
   Hash_OnshellHJJ_nosplit( 44,1:3) = (/-1,-1, 1/)
   Hash_OnshellHJJ_nosplit( 45,1:3) = (/ 4,-1, 1/)
   Hash_OnshellHJJ_nosplit( 46,1:3) = (/-2,-2, 1/)
   Hash_OnshellHJJ_nosplit( 47,1:3) = (/-1,-4, 1/)
   Hash_OnshellHJJ_nosplit( 48,1:3) = (/ 4,-2, 1/)
   Hash_OnshellHJJ_nosplit( 49,1:3) = (/-2,-4, 1/)
   Hash_OnshellHJJ_nosplit( 50,1:3) = (/ 5,-1, 1/)
   Hash_OnshellHJJ_nosplit( 51,1:3) = (/-1,-5, 1/)
   Hash_OnshellHJJ_nosplit( 52,1:3) = (/ 4,-4, 1/)
   Hash_OnshellHJJ_nosplit( 53,1:3) = (/ 4,-4, 2/)
   Hash_OnshellHJJ_nosplit( 54,1:3) = (/ 4,-4, 3/)
   Hash_OnshellHJJ_nosplit( 55,1:3) = (/ 4,-3, 1/)
   Hash_OnshellHJJ_nosplit( 56,1:3) = (/ 4, 3, 1/)
   Hash_OnshellHJJ_nosplit( 57,1:3) = (/ 3,-4, 1/)
   Hash_OnshellHJJ_nosplit( 58,1:3) = (/-2,-5, 1/)
   Hash_OnshellHJJ_nosplit( 59,1:3) = (/-3,-4, 1/)
   Hash_OnshellHJJ_nosplit( 60,1:3) = (/ 5,-2, 1/)
   Hash_OnshellHJJ_nosplit( 61,1:3) = (/ 3, 3, 1/)
   Hash_OnshellHJJ_nosplit( 62,1:3) = (/-3,-3, 1/)
   Hash_OnshellHJJ_nosplit( 63,1:3) = (/ 5, 3, 1/)
   Hash_OnshellHJJ_nosplit( 64,1:3) = (/ 3,-5, 1/)
   Hash_OnshellHJJ_nosplit( 65,1:3) = (/ 5,-3, 1/)
   Hash_OnshellHJJ_nosplit( 66,1:3) = (/-3,-5, 1/)
   Hash_OnshellHJJ_nosplit( 67,1:3) = (/ 5, 4, 1/)
   Hash_OnshellHJJ_nosplit( 68,1:3) = (/ 5,-4, 1/)
   Hash_OnshellHJJ_nosplit( 69,1:3) = (/-4,-5, 1/)
   Hash_OnshellHJJ_nosplit( 70,1:3) = (/ 4,-5, 1/)
   Hash_OnshellHJJ_nosplit( 71,1:3) = (/ 5,-5, 1/)
   Hash_OnshellHJJ_nosplit( 72,1:3) = (/ 5,-5, 2/)
   Hash_OnshellHJJ_nosplit( 73,1:3) = (/ 5,-5, 3/)
   Hash_OnshellHJJ_nosplit( 74,1:3) = (/-4,-4, 1/)
   Hash_OnshellHJJ_nosplit( 75,1:3) = (/ 4, 4, 1/)
   Hash_OnshellHJJ_nosplit( 76,1:3) = (/-5,-5, 1/)
   Hash_OnshellHJJ_nosplit( 77,1:3) = (/ 5, 5, 1/)

   Hash_OnshellHJJ_nosplit( nijchannels+1:,:) = 0
   Hash_OnshellHJJ_nosplit( nijchannels+1:,3) = -1
return
end subroutine

subroutine Init_Hash_GENchannel()
implicit none
   Hash_GENchannel(  1,1:3) = (/-5,-5, 1/)
   Hash_GENchannel(  2,1:3) = (/-5,-4, 1/)
   Hash_GENchannel(  3,1:3) = (/-5,-3, 1/)
   Hash_GENchannel(  4,1:3) = (/-5,-2, 1/)
   Hash_GENchannel(  5,1:3) = (/-5,-1, 1/)
   Hash_GENchannel(  6,1:3) = (/-5, 0, 1/)
   Hash_GENchannel(  7,1:3) = (/-5, 1, 1/)
   Hash_GENchannel(  8,1:3) = (/-5, 2, 1/)
   Hash_GENchannel(  9,1:3) = (/-5, 3, 1/)
   Hash_GENchannel( 10,1:3) = (/-5, 4, 1/)
   Hash_GENchannel( 11,1:3) = (/-5, 5, 1/)
   Hash_GENchannel( 12,1:3) = (/-4,-5, 1/)
   Hash_GENchannel( 13,1:3) = (/-4,-4, 1/)
   Hash_GENchannel( 14,1:3) = (/-4,-3, 1/)
   Hash_GENchannel( 15,1:3) = (/-4,-2, 1/)
   Hash_GENchannel( 16,1:3) = (/-4,-1, 1/)
   Hash_GENchannel( 17,1:3) = (/-4, 0, 1/)
   Hash_GENchannel( 18,1:3) = (/-4, 1, 1/)
   Hash_GENchannel( 19,1:3) = (/-4, 2, 1/)
   Hash_GENchannel( 20,1:3) = (/-4, 3, 1/)
   Hash_GENchannel( 21,1:3) = (/-4, 4, 1/)
   Hash_GENchannel( 22,1:3) = (/-4, 5, 1/)
   Hash_GENchannel( 23,1:3) = (/-3,-5, 1/)
   Hash_GENchannel( 24,1:3) = (/-3,-4, 1/)
   Hash_GENchannel( 25,1:3) = (/-3,-3, 1/)
   Hash_GENchannel( 26,1:3) = (/-3,-2, 1/)
   Hash_GENchannel( 27,1:3) = (/-3,-1, 1/)
   Hash_GENchannel( 28,1:3) = (/-3, 0, 1/)
   Hash_GENchannel( 29,1:3) = (/-3, 1, 1/)
   Hash_GENchannel( 30,1:3) = (/-3, 2, 1/)
   Hash_GENchannel( 31,1:3) = (/-3, 3, 1/)
   Hash_GENchannel( 32,1:3) = (/-3, 4, 1/)
   Hash_GENchannel( 33,1:3) = (/-3, 5, 1/)
   Hash_GENchannel( 34,1:3) = (/-2,-5, 1/)
   Hash_GENchannel( 35,1:3) = (/-2,-4, 1/)
   Hash_GENchannel( 36,1:3) = (/-2,-3, 1/)
   Hash_GENchannel( 37,1:3) = (/-2,-2, 1/)
   Hash_GENchannel( 38,1:3) = (/-2,-1, 1/)
   Hash_GENchannel( 39,1:3) = (/-2, 0, 1/)
   Hash_GENchannel( 40,1:3) = (/-2, 1, 1/)
   Hash_GENchannel( 41,1:3) = (/-2, 2, 1/)
   Hash_GENchannel( 42,1:3) = (/-2, 3, 1/)
   Hash_GENchannel( 43,1:3) = (/-2, 4, 1/)
   Hash_GENchannel( 44,1:3) = (/-2, 5, 1/)
   Hash_GENchannel( 45,1:3) = (/-1,-5, 1/)
   Hash_GENchannel( 46,1:3) = (/-1,-4, 1/)
   Hash_GENchannel( 47,1:3) = (/-1,-3, 1/)
   Hash_GENchannel( 48,1:3) = (/-1,-2, 1/)
   Hash_GENchannel( 49,1:3) = (/-1,-1, 1/)
   Hash_GENchannel( 50,1:3) = (/-1, 0, 1/)
   Hash_GENchannel( 51,1:3) = (/-1, 1, 1/)
   Hash_GENchannel( 52,1:3) = (/-1, 2, 1/)
   Hash_GENchannel( 53,1:3) = (/-1, 3, 1/)
   Hash_GENchannel( 54,1:3) = (/-1, 4, 1/)
   Hash_GENchannel( 55,1:3) = (/-1, 5, 1/)
   Hash_GENchannel( 56,1:3) = (/ 0,-5, 1/)
   Hash_GENchannel( 57,1:3) = (/ 0,-4, 1/)
   Hash_GENchannel( 58,1:3) = (/ 0,-3, 1/)
   Hash_GENchannel( 59,1:3) = (/ 0,-2, 1/)
   Hash_GENchannel( 60,1:3) = (/ 0,-1, 1/)
   Hash_GENchannel( 61,1:3) = (/ 0, 0, 1/)
   Hash_GENchannel( 62,1:3) = (/ 0, 1, 1/)
   Hash_GENchannel( 63,1:3) = (/ 0, 2, 1/)
   Hash_GENchannel( 64,1:3) = (/ 0, 3, 1/)
   Hash_GENchannel( 65,1:3) = (/ 0, 4, 1/)
   Hash_GENchannel( 66,1:3) = (/ 0, 5, 1/)
   Hash_GENchannel( 67,1:3) = (/ 1,-5, 1/)
   Hash_GENchannel( 68,1:3) = (/ 1,-4, 1/)
   Hash_GENchannel( 69,1:3) = (/ 1,-3, 1/)
   Hash_GENchannel( 70,1:3) = (/ 1,-2, 1/)
   Hash_GENchannel( 71,1:3) = (/ 1,-1, 1/)
   Hash_GENchannel( 72,1:3) = (/ 1, 0, 1/)
   Hash_GENchannel( 73,1:3) = (/ 1, 1, 1/)
   Hash_GENchannel( 74,1:3) = (/ 1, 2, 1/)
   Hash_GENchannel( 75,1:3) = (/ 1, 3, 1/)
   Hash_GENchannel( 76,1:3) = (/ 1, 4, 1/)
   Hash_GENchannel( 77,1:3) = (/ 1, 5, 1/)
   Hash_GENchannel( 78,1:3) = (/ 2,-5, 1/)
   Hash_GENchannel( 79,1:3) = (/ 2,-4, 1/)
   Hash_GENchannel( 80,1:3) = (/ 2,-3, 1/)
   Hash_GENchannel( 81,1:3) = (/ 2,-2, 1/)
   Hash_GENchannel( 82,1:3) = (/ 2,-1, 1/)
   Hash_GENchannel( 83,1:3) = (/ 2, 0, 1/)
   Hash_GENchannel( 84,1:3) = (/ 2, 1, 1/)
   Hash_GENchannel( 85,1:3) = (/ 2, 2, 1/)
   Hash_GENchannel( 86,1:3) = (/ 2, 3, 1/)
   Hash_GENchannel( 87,1:3) = (/ 2, 4, 1/)
   Hash_GENchannel( 88,1:3) = (/ 2, 5, 1/)
   Hash_GENchannel( 89,1:3) = (/ 3,-5, 1/)
   Hash_GENchannel( 90,1:3) = (/ 3,-4, 1/)
   Hash_GENchannel( 91,1:3) = (/ 3,-3, 1/)
   Hash_GENchannel( 92,1:3) = (/ 3,-2, 1/)
   Hash_GENchannel( 93,1:3) = (/ 3,-1, 1/)
   Hash_GENchannel( 94,1:3) = (/ 3, 0, 1/)
   Hash_GENchannel( 95,1:3) = (/ 3, 1, 1/)
   Hash_GENchannel( 96,1:3) = (/ 3, 2, 1/)
   Hash_GENchannel( 97,1:3) = (/ 3, 3, 1/)
   Hash_GENchannel( 98,1:3) = (/ 3, 4, 1/)
   Hash_GENchannel( 99,1:3) = (/ 3, 5, 1/)
   Hash_GENchannel(100,1:3) = (/ 4,-5, 1/)
   Hash_GENchannel(101,1:3) = (/ 4,-4, 1/)
   Hash_GENchannel(102,1:3) = (/ 4,-3, 1/)
   Hash_GENchannel(103,1:3) = (/ 4,-2, 1/)
   Hash_GENchannel(104,1:3) = (/ 4,-1, 1/)
   Hash_GENchannel(105,1:3) = (/ 4, 0, 1/)
   Hash_GENchannel(106,1:3) = (/ 4, 1, 1/)
   Hash_GENchannel(107,1:3) = (/ 4, 2, 1/)
   Hash_GENchannel(108,1:3) = (/ 4, 3, 1/)
   Hash_GENchannel(109,1:3) = (/ 4, 4, 1/)
   Hash_GENchannel(110,1:3) = (/ 4, 5, 1/)
   Hash_GENchannel(111,1:3) = (/ 5,-5, 1/)
   Hash_GENchannel(112,1:3) = (/ 5,-4, 1/)
   Hash_GENchannel(113,1:3) = (/ 5,-3, 1/)
   Hash_GENchannel(114,1:3) = (/ 5,-2, 1/)
   Hash_GENchannel(115,1:3) = (/ 5,-1, 1/)
   Hash_GENchannel(116,1:3) = (/ 5, 0, 1/)
   Hash_GENchannel(117,1:3) = (/ 5, 1, 1/)
   Hash_GENchannel(118,1:3) = (/ 5, 2, 1/)
   Hash_GENchannel(119,1:3) = (/ 5, 3, 1/)
   Hash_GENchannel(120,1:3) = (/ 5, 4, 1/)
   Hash_GENchannel(121,1:3) = (/ 5, 5, 1/)
return
end subroutine

subroutine Init_Hash_THchannel()
implicit none
   Hash_THchannel(  1,1:3) = (/-5,-4,111/)
   Hash_THchannel(  2,1:3) = (/-5,-2,111/)
   Hash_THchannel(  3,1:3) = (/-5, 1,111/)
   Hash_THchannel(  4,1:3) = (/-5, 3,111/)
   Hash_THchannel(  5,1:3) = (/-4,-5,111/)
   Hash_THchannel(  6,1:3) = (/-4, 3,113/)
   Hash_THchannel(  7,1:3) = (/-3, 4,112/)
   Hash_THchannel(  8,1:3) = (/-3, 5,110/)
   Hash_THchannel(  9,1:3) = (/-2,-5,111/)
   Hash_THchannel( 10,1:3) = (/-2, 1,113/)
   Hash_THchannel( 11,1:3) = (/-1, 2,112/)
   Hash_THchannel( 12,1:3) = (/-1, 5,110/)
   Hash_THchannel( 13,1:3) = (/ 1,-5,111/)
   Hash_THchannel( 14,1:3) = (/ 1,-2,113/)
   Hash_THchannel( 15,1:3) = (/ 2,-1,112/)
   Hash_THchannel( 16,1:3) = (/ 2, 5,110/)
   Hash_THchannel( 17,1:3) = (/ 3,-5,111/)
   Hash_THchannel( 18,1:3) = (/ 3,-4,113/)
   Hash_THchannel( 19,1:3) = (/ 4,-3,112/)
   Hash_THchannel( 20,1:3) = (/ 4, 5,110/)
   Hash_THchannel( 21,1:3) = (/ 5,-3,110/)
   Hash_THchannel( 22,1:3) = (/ 5,-1,110/)
   Hash_THchannel( 23,1:3) = (/ 5, 2,110/)
   Hash_THchannel( 24,1:3) = (/ 5, 4,110/)

   Hash_THchannel( 25:,:) = 0
   Hash_THchannel( 25:,3) = -1
return
end subroutine

subroutine Init_Hash_PPXchannel()
implicit none
integer :: i
   Hash_PPXchannel = 0
   Hash_PPXchannel(:,3) = -1

   Hash_PPXchannel(1,1:3) = (/0,0,1/)
   do i=1,5
      Hash_PPXchannel(2*i  ,1:3) = (/-i,i,1/)
      Hash_PPXchannel(2*i+1,1:3) = (/i,-i,1/)
   enddo
return
end subroutine


subroutine SetupHashes()
implicit none
   if (.not. hashcoll_hashes_initialized) then
      call Init_Hash_MCFM_qqVVqq()
      call Init_Hash_MCFM_qqVVqq_Generation()
      call Init_Hash_OnshellVBF()
      call Init_Hash_OnshellVBF_nosplit()
      call Init_Hash_OnshellHJJ()
      call Init_Hash_OnshellHJJ_nosplit()
      call Init_Hash_GENchannel()
      call Init_Hash_THchannel()
      call Init_Hash_PPXchannel()

      hashcoll_hashes_initialized = .true.
   endif
end subroutine

! Copy functions
subroutine get_MCFM_qqVVqq_Hash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_MCFM_qqVVqq
return
end subroutine

subroutine get_MCFM_qqVVqq_GenHash(ijSel)
implicit none
integer, intent(out) :: ijSel(1:Hash_MCFM_qqVVqq_Gen_Size,1:4)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_MCFM_qqVVqq_Gen
return
end subroutine

subroutine get_VBFchannelHash(ijSel)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_OnshellVBF
return
end subroutine

subroutine get_VBFchannelHash_nosplit(ijSel,nijchannels)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
integer, intent(out) :: nijchannels
   nijchannels = 68
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_OnshellVBF_nosplit
return
end subroutine

subroutine get_HJJchannelHash(ijSel)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_OnshellHJJ
return
end subroutine

subroutine get_HJJchannelHash_nosplit(ijSel,nijchannels)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
integer, intent(out) :: nijchannels
   nijchannels=77
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_OnshellHJJ_nosplit
return
end subroutine

subroutine get_GENchannelHash(ijSel)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_GENchannel
return
end subroutine

subroutine get_THchannelHash(ijSel)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_THchannel
return
end subroutine

SUBROUTINE get_PPXchannelHash(ijSel)
implicit none
integer, intent(out) :: ijSel(1:121,1:3)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel=Hash_PPXchannel
RETURN
END SUBROUTINE


! Reference functions
subroutine getRef_MCFM_qqVVqq_Hash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_MCFM_qqVVqq
return
end subroutine

subroutine getRef_MCFM_qqVVqq_GenHash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_MCFM_qqVVqq_Gen
return
end subroutine

subroutine getRef_VBFchannelHash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_OnshellVBF
return
end subroutine

subroutine getRef_VBFchannelHash_nosplit(ijSel,nijchannels)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
integer, intent(out) :: nijchannels
   nijchannels = 68
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_OnshellVBF_nosplit
return
end subroutine

subroutine getRef_HJJchannelHash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_OnshellHJJ
return
end subroutine

subroutine getRef_HJJchannelHash_nosplit(ijSel,nijchannels)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
integer, intent(out) :: nijchannels
   nijchannels=77
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_OnshellHJJ_nosplit
return
end subroutine

subroutine getRef_GENchannelHash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_GENchannel
return
end subroutine

subroutine getRef_THchannelHash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_THchannel
return
end subroutine

SUBROUTINE getRef_PPXchannelHash(ijSel)
implicit none
integer, pointer, intent(out) :: ijSel(:,:)
   if (.not. hashcoll_hashes_initialized) then
      call SetupHashes()
   endif
   ijSel => Hash_PPXchannel
RETURN
END SUBROUTINE


subroutine removeOffshellChannelFromHash(ijSel,iremove,imax,jmax)
implicit none
integer, intent(in) :: imax,jmax,iremove
integer :: ijSel(:,:),k
   if (iremove.gt.0 .and. iremove .lt. imax) then
      do k=iremove+1,imax
         ijSel(k-1,1:jmax) = ijSel(k,1:jmax)
      enddo
      ijSel(imax,:)=Not_a_particle_
   elseif (iremove .eq. imax) then
      ijSel(imax,:)=Not_a_particle_
   endif
end subroutine

subroutine removeOffshellChannelFromHashRef(ijSel,iremove,imax,jmax)
implicit none
integer, intent(in) :: imax,jmax,iremove
integer, pointer :: ijSel(:,:)
integer :: k
   if (iremove.gt.0 .and. iremove .lt. imax) then
      do k=iremove+1,imax
         ijSel(k-1,1:jmax) = ijSel(k,1:jmax)
      enddo
      ijSel(imax,:)=Not_a_particle_
   elseif (iremove .eq. imax) then
      ijSel(imax,:)=Not_a_particle_
   endif
end subroutine


END MODULE
