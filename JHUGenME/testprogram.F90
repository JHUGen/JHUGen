PROGRAM TestProgram
use modHiggs
use modZprime
use modGraviton
implicit none
real(8) :: p(4,6),MatElSq,M_Reso,Ga_Reso
integer :: MY_IDUP(6:9)
real(8),  parameter :: GeV=1d0/100d0
complex(8) :: Hggcoupl(1:3),Hzzcoupl(1:4)
complex(8) :: Zqqcoupl(1:2),Zzzcoupl(1:2)
complex(8) :: Gggcoupl(1:5),Gqqcoupl(1:2),Gzzcoupl(1:10)

! input unit = GeV/100 such that 125GeV is 1.25 in the code
  M_Reso  = 125d0 * GeV
  Ga_Reso = 0.1d0 * GeV
  Hggcoupl(1:3) = (/ (1d0,0d0), (0d0,0d0), (0d0,0d0) /)
  Hzzcoupl(1:4) = (/ (1d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0) /)
  Zqqcoupl(1:2) = (/ (1d0,0d0), (1d0,0d0) /)
  Zzzcoupl(1:2) = (/ (0d0,0d0), (1d0,0d0) /)
  Gggcoupl(1:5) = (/ (1d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0) /)
  Gqqcoupl(1:2) = (/ (1d0,0d0), (1d0,0d0) /)
  Gzzcoupl(1:10)= (/ (1d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0), (1d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0), (0d0,0d0) /)
 
 
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
   print *, ""






! this checks the decay rates without production dynamics
  print *, "PS point 3 (no production dynamics)"
  p(1:4,1) = (/   -1.2500000000000000d0,    0.0000000000000000d0,    0.0000000000000000d0,    0.0000000000000000d0   /)
  p(1:4,2) = (/    1.0000000000000000d9,    1.0000000000000000d9,    1.0000000000000000d9,    1.0000000000000000d9   /)  ! DUMMY
  p(1:4,3) = (/    0.5723124900092045d0,    0.2083223910863685d0,    0.1441561361468674d0,   -0.5131884410270751d0   /)
  p(1:4,4) = (/    0.3270364795109065d0,   -0.1343151455361828d0,    0.0368904244639839d0,    0.2958908535141779d0   /)
  p(1:4,5) = (/    0.0584259632074581d0,   -0.0206108978979021d0,   -0.0424960639382021d0,   -0.0343928570247043d0   /)
  p(1:4,6) = (/    0.2922250672724307d0,   -0.0533963476522836d0,   -0.1385504966726492d0,    0.2516904445376015d0   /)

  call EvalAmp_Zprime_VV(p(1:4,1:6),M_Reso,Ga_Reso,Zzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (spin-1)",MatElSq
  print *, "result should be (spin-1)",5.25259158248293146d-010
  print *, "ratio",MatElSq/5.25259158248293146d-010
  print *, ""
  call EvalAmp_G_VV(p(1:4,1:6),M_Reso,Ga_Reso,Gzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (spin-2)",MatElSq
  print *, "result should be (spin-2)",3.50330723427981412d-009
  print *, "ratio",MatElSq/3.50330723427981412d-009
  print *, ""


  print *, "PS point 4 (no production dynamics)"
  p(1:4,1) = (/   -1.2500000000000000d0,    0.0000000000000000d0,    0.0000000000000000d0,    0.0000000000000000d0   /)
  p(1:4,2) = (/    1.0000000000000000d9,    1.0000000000000000d9,    1.0000000000000000d9,    1.0000000000000000d9   /)  ! DUMMY
  p(1:4,3) = (/    0.5501955504289553d0,    0.1590232342830207d0,    0.1201521110389673d0,   -0.5128257256445587d0   /)
  p(1:4,4) = (/    0.1226104406534881d0,   -0.0150215063102395d0,   -0.1053277769957477d0,    0.0609404126877065d0   /)
  p(1:4,5) = (/    0.1283008077137973d0,   -0.0450661397083963d0,    0.1152711293482556d0,    0.0338039502214440d0   /)
  p(1:4,6) = (/    0.4488932012037591d0,   -0.0989355882643850d0,   -0.1300954633914751d0,    0.4180813627354082d0   /)

  call EvalAmp_Zprime_VV(p(1:4,1:6),M_Reso,Ga_Reso,Zzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (spin-1)",MatElSq
  print *, "result should be (spin-1)",2.41036103083747314d-011
  print *, "ratio",MatElSq/2.41036103083747314d-011
  print *, ""
  call EvalAmp_G_VV(p(1:4,1:6),M_Reso,Ga_Reso,Gzzcoupl,MY_IDUP(6:9),MatElSq)
  print *, "Matr.el. squared (spin-2)",MatElSq
  print *, "result should be (spin-2)",1.68408821989468668d-010
  print *, "ratio",MatElSq/1.68408821989468668d-010



END PROGRAM
