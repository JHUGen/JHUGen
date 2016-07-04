PROGRAM TestProgram
use modHiggs
use modZprime
use modGraviton
implicit none
real(8) :: p(4,6),MatElSq,M_Reso,Ga_Reso
integer :: MY_IDUP(6:9)
real(8),  parameter :: GeV=1d0/100d0
real(8) :: Hggcoupl(1:3),Hzzcoupl(1:4),Zqqcoupl(1:2),Zzzcoupl(1:2),Gggcoupl(1:5),Gqqcoupl(1:2),Gzzcoupl(1:10)


! input unit = GeV/100 such that 125GeV is 1.25 in the code
  M_Reso  = 125d0 * GeV
  Ga_Reso = 0.1d0 * GeV
  Hggcoupl(1:3) = (/1d0,0d0,0d0/)
  Hzzcoupl(1:4) = (/1d0,0d0,0d0,0d0/)
  Zqqcoupl(1:2) = (/1d0,1d0/)
  Zzzcoupl(1:2) = (/0d0,1d0/)
  Gggcoupl(1:5) = (/1d0,0d0,0d0,0d0,0d0/)
  Gqqcoupl(1:2) = (/1d0,1d0/)
  Gzzcoupl(1:10)= (/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,0d0/)
 
 
! particle ID: +7=e+,  -7=e-,  +8=mu+,  -8=mu-
  MY_IDUP(6:9) = (/+7,-7,+7,-7/)

! p(1:4,i) = (E(i),px(i),py(i),pz(i))
! i=1,2: glu1,glu2 (outgoing convention)
! i=3,4: correspond to MY_IDUP(7),MY_IDUP(6)
! i=5,6: correspond to MY_IDUP(9),MY_IDUP(8)
  p(1:4,1) = (/         -0.0645264033200954d0,          0.0000000000000000d0,          0.0000000000000000d0,         -0.0645264033200954d0   /)
  p(1:4,2) = (/         -6.0537234356954572d0,          0.0000000000000000d0,          0.0000000000000000d0,          6.0537234356954572d0   /)
  p(1:4,3) = (/          4.7598878302889442d0,          0.4544925597087586d0,          0.0917597774970785d0,         -4.7372511874858247d0   /)
  p(1:4,4) = (/          0.8177159853499207d0,         -0.2768889802512220d0,         -0.0202805643015353d0,         -0.7691427851991084d0   /)
  p(1:4,5) = (/          0.3191706236713395d0,         -0.0470102494218651d0,         -0.0079466854602927d0,         -0.3155895651859249d0   /)
  p(1:4,6) = (/          0.2214753997053482d0,         -0.1305933300356715d0,         -0.0635325277352505d0,         -0.1672134945045036d0   /)


! note: numbers correspond to includeInterference=.false.
! spin-0 corresponds to SM Higgs; spin-2 has minimal couplings
  print *, "PS point 1"

  call EvalAmp_gg_H_VV(p(1:4,1:6),M_Reso,Ga_Reso,Hggcoupl,Hzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (spin-0)",MatElSq
  print *, "result should be (spin-0)",0.0045682366425370D0
  print *, "ratio",MatElSq/0.0045682366425370D0
  print *, ""
  call EvalAmp_qqb_Zprime_VV(p(1:4,1:6),M_Reso,Ga_Reso,Zqqcoupl,Zzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (spin-1)",MatElSq
  print *, "result should be (spin-1)",0.0020357978978982D0
  print *, "ratio",MatElSq/0.0020357978978982D0
  print *, ""
  call EvalAmp_gg_G_VV(p(1:4,1:6),M_Reso,Ga_Reso,Gggcoupl,Gzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (gg spin-2)",MatElSq
  print *, "result should be (gg spin-2)",0.0307096869320374d0
  print *, "ratio",MatElSq/0.0307096869320374d0
  print *, ""
  call EvalAmp_qqb_G_VV(p(1:4,1:6),M_Reso,Ga_Reso,Gqqcoupl,Gzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (qq spin-2)",MatElSq
  print *, "result should be (qq spin-2)",0.0004838377647021d0
  print *, "ratio",MatElSq/0.0004838377647021d0




!   MY_IDUP(6:9) = (/+7,-7,+8,-8/)
!   p(1:4,1) = (/ -0.0370152249548727d0,          0.0000000000000000d0,          0.0000000000000000d0,         -0.0370152249548727d0  /)
!   p(1:4,2) = (/ -10.5530899913815492d0,          0.0000000000000000d0,          0.0000000000000000d0,         10.5530899913815492d0  /)
!   p(1:4,3) = (/ 8.5411802617610935d0,          0.1767561504734450d0,          0.3952687972845690d0,         -8.5301981281245940d0 /)
!   p(1:4,4) = (/ 0.7050609986225407d0,          0.0155515571468995d0,         -0.2232185028886611d0,         -0.6686124892769094d0 /)
!   p(1:4,5) = (/ 1.2474949073683825d0,         -0.1650712873082496d0,         -0.1508369698905747d0,         -1.2272910097163903d0/)
!   p(1:4,6) = (/ 0.0963690485844057d0,         -0.0272364203120949d0,         -0.0212133245053332d0,         -0.0899731393087841d0 /)
!    
! 
!   print *, ""
!   print *, "PS point 2"
! 
!   call EvalAmp_gg_H_VV(p(1:4,1:6),M_Reso,Ga_Reso,MY_IDUP(6:9),MatElSq)
!   print *, "Matr.el. squared (spin-0)",MatElSq
!   print *, "result should be (spin-0)",0.0339676512275350d0
!   print *, "ratio",MatElSq/0.0339676512275350d0
!   print *, ""
!   call EvalAmp_gg_G_VV(p(1:4,1:6),M_Reso,Ga_Reso,MY_IDUP(6:9),MatElSq)
!   print *, "Matr.el. squared (gg spin-2)",MatElSq
!   print *, "result should be (gg spin-2)",0.0112810555699413d0
!   print *, "ratio",MatElSq/0.0112810555699413d0
!   print *, ""
!   call EvalAmp_qqb_G_VV(p(1:4,1:6),M_Reso,Ga_Reso,MY_IDUP(6:9),MatElSq)
!   print *, "Matr.el. squared (qq spin-2)",MatElSq
!   print *, "result should be (qq spin-2)",0.0112810555699413d0
!   print *, "ratio",MatElSq/0.0112810555699413d0

  


END PROGRAM
