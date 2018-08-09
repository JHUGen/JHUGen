!--YaofuZhou-----------------------------------------
module ModHH
  use ModParameters
!  use ModKinematics
  use ModMisc
  use ModVHaux
#if useCollier==1
  use Collier
  use ModggboxHH1mm
  use ModggboxHH1mp
  use ModggboxHH1pm
  use ModggboxHH1pp
  use ModggboxHH2mm
  use ModggboxHH2mp
  use ModggboxHH2pm
  use ModggboxHH2pp
  use ModggboxHH3mm
  use ModggboxHH3mp
  use ModggboxHH3pm
  use ModggboxHH3pp
  use ModggboxHH4mm
  use ModggboxHH4mp
  use ModggboxHH4pm
  use ModggboxHH4pp
  use ModggboxHH5mm
  use ModggboxHH5mp
  use ModggboxHH5pm
  use ModggboxHH5pp
  use ModggboxHH6mm
  use ModggboxHH6mp
  use ModggboxHH6pm
  use ModggboxHH6pp

  use Modggbox5HH1mm
  use Modggbox5HH1mp
  use Modggbox5HH1pm
  use Modggbox5HH1pp
  use Modggbox5HH2mm
  use Modggbox5HH2mp
  use Modggbox5HH2pm
  use Modggbox5HH2pp
  use Modggbox5HH3mm
  use Modggbox5HH3mp
  use Modggbox5HH3pm
  use Modggbox5HH3pp
  use Modggbox5HH4mm
  use Modggbox5HH4mp
  use Modggbox5HH4pm
  use Modggbox5HH4pp
  use Modggbox5HH5mm
  use Modggbox5HH5mp
  use Modggbox5HH5pm
  use Modggbox5HH5pp
  use Modggbox5HH6mm
  use Modggbox5HH6mp
  use Modggbox5HH6pm
  use Modggbox5HH6pp
#endif
  implicit none


  !----- notation for subroutines
  public :: Amp_HH

contains

subroutine Amp_HH(Mom,mass,helicity,id,amp)
  implicit none
  integer, intent(in) :: id(9)
  real(8), intent(in) :: helicity(9)
  real(8), intent(in) :: Mom(1:4,1:9)
  complex(8), intent(out) :: amp
  real(8) :: mass(3:5,1:2)
  complex(8) :: PROP3, PROP4, PROP5
  real(8) :: q3_q3, q4_q4, q5_q5
  complex(8) :: gHHH
  complex(8) :: TriH, BoxHH, Box5HH
  complex(8) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8) :: sprod(1:4,1:4)
  integer :: i

#if useCollier==1
  amp=czero

  q3_q3 = scr(Mom(:,3),Mom(:,3))
  q4_q4 = scr(Mom(:,4),Mom(:,4))
  q5_q5 = scr(Mom(:,5),Mom(:,5))

  PROP3 = -PROPAGATOR(q3_q3,mass(3,1),mass(3,2))
  PROP4 = -PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
  PROP5 = -PROPAGATOR(q5_q5,mass(5,1),mass(5,2))

  gHHH = -ci*3d0*(M_reso**2)/vev

  call spinoru2(4,(/Mom(1:4,1),Mom(1:4,2),Mom(1:4,6),Mom(1:4,7)/),Spaa,Spbb,sprod)
  call ggTriH(Mom,Spaa,Spbb,sprod,helicity,TriH)
  call ggBoxHH(Mom,Spaa,Spbb,sprod,helicity,BoxHH)
  call ggBox5HH(Mom,Spaa,Spbb,sprod,helicity,Box5HH)

  amp = kappa*BoxHH + ci*kappa_tilde*Box5HH + TriH*PROP3*gHHH
!print*, gHHH
!pause
  amp = amp * (-1d0) * gs**2 !-1 = i^2 from g_s each
!print*, gs, amp
  amp = amp * PROP4                                                                  &
  *(    kappa      *FFS(id(6), Mom(:,6), helicity(6), id(7), Mom(:,7), helicity(7))  &
    +ci*kappa_tilde*FFP(id(6), Mom(:,6), helicity(6), id(7), Mom(:,7), helicity(7))) &
  *(-ci*massfrun(getMass(convertLHEreverse(id(6))),get_minv(Mom(:,4)))/vev)

  amp = amp * PROP5                                                                  &
  *(    kappa      *FFS(id(8), Mom(:,8), helicity(8), id(9), Mom(:,9), helicity(9))  &
    +ci*kappa_tilde*FFP(id(8), Mom(:,8), helicity(8), id(9), Mom(:,9), helicity(9))) &
  *(-ci*massfrun(getMass(convertLHEreverse(id(8))),get_minv(Mom(:,5)))/vev)
!print*,PROP3
!print*,PROP4
!print*,PROP5
!print*,kappa*      FFS(id(6), Mom(:,6), helicity(6), id(7), Mom(:,7), helicity(7))
!print*,kappa*      FFS(id(8), Mom(:,8), helicity(8), id(9), Mom(:,9), helicity(9))
!print*,(-ci/vev*massfrun(getMass(convertLHEreverse(id(6))),get_minv(Mom(:,4))))
!print*,(-ci/vev*massfrun(getMass(convertLHEreverse(id(8))),get_minv(Mom(:,5))))
!print*,massfrun(getMass(convertLHEreverse(id(6))),get_minv(Mom(:,4)))
!print*,"-------"
  return
#else
print *, "Need to link COLLIER for this process."
print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
print *, "specified in the makefile."
stop 1
#endif
end subroutine Amp_HH





subroutine ggTriH(Mom,Spaa,Spbb,sprod,helicity,TriH)
  implicit none
  real(8) , intent(in) :: Mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  real(8) :: m_top_run
  complex(8) C0
  complex(8) TriH

#if useCollier==1
  call C0_cll(C0,dcmplx(0d0),dcmplx(sprod(1,2)),dcmplx(0d0),dcmplx(M_Top**2),dcmplx(M_Top**2),dcmplx(M_Top**2))
  C0 = C0 * ci * pisq / (2d0*pi)**4

  !m_top_run = massfrun(m_top,Mu_Ren)
  m_top_run = m_top

  if(    helicity(1).gt.0d0 .and. helicity(2).gt.0d0)then
    TriH= -2d0 * m_top_run * Spbb(1,2) * (2d0 + C0 * (4d0 * m_top_run**2 + Spaa(1,2) * Spbb(1,2))) / Spaa(1,2)
  elseif(helicity(1).lt.0d0 .and. helicity(2).lt.0d0)then
    TriH= -2d0 * m_top_run * Spaa(1,2) * (2d0 + C0 * (4d0 * m_top_run**2 + Spaa(1,2) * Spbb(1,2))) / Spbb(1,2)
  else
    TriH= 0d0
  endif

  TriH = 2d0 * TriH * (-ci) * kappa * m_top_run/vev !2 for 2 digrams with opposite loop momentum

  return
#else
print *, "Need to link COLLIER for this process."
print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
print *, "specified in the makefile."
stop 1
#endif
end subroutine ggTriH



subroutine ggBoxHH(Mom,Spaa,Spbb,sprod,helicity,BoxHH)
  implicit none
  real(8) , intent(in) :: Mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  real(8) :: m_top_run
  complex(8) BoxHH
  complex(8) ggboxHH1,ggboxHH2,ggboxHH3,ggboxHH4,ggboxHH5,ggboxHH6  

#if useCollier==1
  if(    helicity(1).gt.0d0.and.helicity(2).lt.0d0)then!+-
    call ggboxHH1pm(Mom,Spaa,Spbb,sprod,ggboxHH1)
    !call ggboxHH2pm(Mom,Spaa,Spbb,sprod,ggboxHH2)
    call ggboxHH3pm(Mom,Spaa,Spbb,sprod,ggboxHH3)
    !call ggboxHH4pm(Mom,Spaa,Spbb,sprod,ggboxHH4)
    call ggboxHH5pm(Mom,Spaa,Spbb,sprod,ggboxHH5)
    call ggboxHH6pm(Mom,Spaa,Spbb,sprod,ggboxHH6)
  elseif(helicity(1).lt.0d0.and.helicity(2).gt.0d0)then!-+
    call ggboxHH1mp(Mom,Spaa,Spbb,sprod,ggboxHH1)
    !call ggboxHH2mp(Mom,Spaa,Spbb,sprod,ggboxHH2)
    call ggboxHH3mp(Mom,Spaa,Spbb,sprod,ggboxHH3)
    !call ggboxHH4mp(Mom,Spaa,Spbb,sprod,ggboxHH4)
    call ggboxHH5mp(Mom,Spaa,Spbb,sprod,ggboxHH5)
    call ggboxHH6mp(Mom,Spaa,Spbb,sprod,ggboxHH6)
  elseif(helicity(1).gt.0d0.and.helicity(2).gt.0d0)then!++
    call ggboxHH1pp(Mom,Spaa,Spbb,sprod,ggboxHH1)
    !call ggboxHH2pp(Mom,Spaa,Spbb,sprod,ggboxHH2)
    call ggboxHH3pp(Mom,Spaa,Spbb,sprod,ggboxHH3)
    !call ggboxHH4pp(Mom,Spaa,Spbb,sprod,ggboxHH4)
    call ggboxHH5pp(Mom,Spaa,Spbb,sprod,ggboxHH5)
    call ggboxHH6pp(Mom,Spaa,Spbb,sprod,ggboxHH6)
  else                                                 !--
    call ggboxHH1mm(Mom,Spaa,Spbb,sprod,ggboxHH1)
    !call ggboxHH2mm(Mom,Spaa,Spbb,sprod,ggboxHH2)
    call ggboxHH3mm(Mom,Spaa,Spbb,sprod,ggboxHH3)
    !call ggboxHH4mm(Mom,Spaa,Spbb,sprod,ggboxHH4)
    call ggboxHH5mm(Mom,Spaa,Spbb,sprod,ggboxHH5)
    call ggboxHH6mm(Mom,Spaa,Spbb,sprod,ggboxHH6)
  endif

!print*,helicity(1:2)
!print*,ggboxHH1
!print*,ggboxHH2
!print*,ggboxHH3
!print*,ggboxHH4
!if(helicity(1).eq.helicity(2))then
!!print*,helicity(1:2)
!call ggboxHH5pp(Mom,Spaa,Spbb,sprod,ggboxHH5)
!call ggboxHH6pp(Mom,Spaa,Spbb,sprod,ggboxHH6)
!print*,ggboxHH5
!print*,ggboxHH6
!call ggboxHH5mm(Mom,Spaa,Spbb,sprod,ggboxHH5)
!call ggboxHH6mm(Mom,Spaa,Spbb,sprod,ggboxHH6)
!print*,ggboxHH5
!print*,ggboxHH6
!print*,"-----------"

!endif
  BoxHH = ggboxHH1*2d0 + ggboxHH3*2d0 + ggboxHH5 + ggboxHH6
  BoxHH = -BoxHH * ci * pisq / (2d0*pi)**4 ! -1 = -1 (fermion loop) i^4 (4 fermion propagators)

  return
#else
print *, "Need to link COLLIER for this process."
print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
print *, "specified in the makefile."
stop 1
#endif
end subroutine ggBoxHH





subroutine ggBox5HH(Mom,Spaa,Spbb,sprod,helicity,Box5HH)
  implicit none
  real(8) , intent(in) :: Mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  real(8) :: m_top_run
  complex(8) Box5HH
  complex(8) ggbox5HH1,ggbox5HH2,ggbox5HH3,ggbox5HH4,ggbox5HH5,ggbox5HH6  

#if useCollier==1
  if(    helicity(1).gt.0d0.and.helicity(2).lt.0d0)then!+-
    call ggbox5HH1pm(Mom,Spaa,Spbb,sprod,ggbox5HH1)
    !call ggbox5HH2pm(Mom,Spaa,Spbb,sprod,ggbox5HH2)
    call ggbox5HH3pm(Mom,Spaa,Spbb,sprod,ggbox5HH3)
    !call ggbox5HH4pm(Mom,Spaa,Spbb,sprod,ggbox5HH4)
    call ggbox5HH5pm(Mom,Spaa,Spbb,sprod,ggbox5HH5)
    call ggbox5HH6pm(Mom,Spaa,Spbb,sprod,ggbox5HH6)
  elseif(helicity(1).lt.0d0.and.helicity(2).gt.0d0)then!-+
    call ggbox5HH1mp(Mom,Spaa,Spbb,sprod,ggbox5HH1)
    !call ggbox5HH2mp(Mom,Spaa,Spbb,sprod,ggbox5HH2)
    call ggbox5HH3mp(Mom,Spaa,Spbb,sprod,ggbox5HH3)
    !call ggbox5HH4mp(Mom,Spaa,Spbb,sprod,ggbox5HH4)
    call ggbox5HH5mp(Mom,Spaa,Spbb,sprod,ggbox5HH5)
    call ggbox5HH6mp(Mom,Spaa,Spbb,sprod,ggbox5HH6)
  elseif(helicity(1).gt.0d0.and.helicity(2).gt.0d0)then!++
    call ggbox5HH1pp(Mom,Spaa,Spbb,sprod,ggbox5HH1)
    !call ggbox5HH2pp(Mom,Spaa,Spbb,sprod,ggbox5HH2)
    call ggbox5HH3pp(Mom,Spaa,Spbb,sprod,ggbox5HH3)
    !call ggbox5HH4pp(Mom,Spaa,Spbb,sprod,ggbox5HH4)
    call ggbox5HH5pp(Mom,Spaa,Spbb,sprod,ggbox5HH5)
    call ggbox5HH6pp(Mom,Spaa,Spbb,sprod,ggbox5HH6)
  else                                                 !--
    call ggbox5HH1mm(Mom,Spaa,Spbb,sprod,ggbox5HH1)
    !call ggbox5HH2mm(Mom,Spaa,Spbb,sprod,ggbox5HH2)
    call ggbox5HH3mm(Mom,Spaa,Spbb,sprod,ggbox5HH3)
    !call ggbox5HH4mm(Mom,Spaa,Spbb,sprod,ggbox5HH4)
    call ggbox5HH5mm(Mom,Spaa,Spbb,sprod,ggbox5HH5)
    call ggbox5HH6mm(Mom,Spaa,Spbb,sprod,ggbox5HH6)
  endif

!print*,helicity(1:2)
!print*,ggbox5HH1
!print*,ggbox5HH2
!print*,ggbox5HH3
!print*,ggbox5HH4
!if(helicity(1).eq.helicity(2))then
!!print*,helicity(1:2)
!call ggboxHH5pp(Mom,Spaa,Spbb,sprod,ggboxHH5)
!call ggboxHH6pp(Mom,Spaa,Spbb,sprod,ggboxHH6)
!print*,ggbox5HH5
!print*,ggbox5HH6
!call ggboxHH5mm(Mom,Spaa,Spbb,sprod,ggboxHH5)
!call ggboxHH6mm(Mom,Spaa,Spbb,sprod,ggboxHH6)
!print*,ggboxHH5
!print*,ggboxHH6
!print*,"==========="
!pause
!endif
  Box5HH = ggbox5HH1*2d0 + ggbox5HH3*2d0 + ggbox5HH5 + ggbox5HH6
  Box5HH = -Box5HH * ci * pisq / (2d0*pi)**4 ! -1 = -1 (fermion loop) i^4 (4 fermion propagators)

  return
#else
print *, "Need to link COLLIER for this process."
print *, "Please set either linkMELA or linkCollierLib to Yes in the makefile and recompile"
print *, "You will have to have a compiled JHUGenMELA or a compiled COLLIER in the directories"
print *, "specified in the makefile."
stop 1
#endif
end subroutine ggBox5HH



end module ModHH
!!--YaofuZhou-----------------------------------------
