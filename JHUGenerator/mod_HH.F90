!--YaofuZhou-----------------------------------------
module ModHH
  use ModParameters
  use ModKinematics
  use ModMisc
  use ModVHaux
  use Collier
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
  complex(8) :: TriH, BoxHH
  complex(8) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8) :: sprod(1:4,1:4)
  integer :: i

  amp=czero

  q3_q3 = scr(Mom(:,3),Mom(:,3))
  q4_q4 = scr(Mom(:,4),Mom(:,4))
  q5_q5 = scr(Mom(:,5),Mom(:,5))

  PROP3 = -PROPAGATOR(q3_q3,mass(3,1),mass(3,2))
  PROP4 = -PROPAGATOR(q4_q4,mass(4,1),mass(4,2))
  PROP5 = -PROPAGATOR(q5_q5,mass(5,1),mass(5,2))

  gHHH = -ci*3d0/2d0*dsqrt(gwsq)*(M_reso**2)/(M_W**2)
!gHHH = -ci*3d0*M_W/vev*M_reso**2/M_W**2

  call spinoru2(4,(/Mom(1:4,1),Mom(1:4,2),Mom(1:4,6),Mom(1:4,7)/),Spaa,Spbb,sprod)
  call ggTriH(Mom,Spaa,Spbb,sprod,helicity,TriH)
  call ggBoxHH(Mom,Spaa,Spbb,sprod,helicity,BoxHH)

  amp = BoxHH + TriH*PROP3*gHHH
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
end subroutine ggTriH



subroutine ggBoxHH(Mom,Spaa,Spbb,sprod,helicity,BoxHH)
  implicit none
  real(8) , intent(in) :: Mom(1:4,1:9)
  complex(8), intent(in) :: Spaa(1:4,1:4),Spbb(1:4,1:4)
  real(8), intent(in) :: sprod(1:4,1:4)
  real(8), intent(in) :: helicity(9)
  real(8) :: m_top_run
  complex(8) C0
  complex(8) BoxHH

  BoxHH=0d0

  return
end subroutine ggBoxHH






end module ModHH
!!--YaofuZhou-----------------------------------------
