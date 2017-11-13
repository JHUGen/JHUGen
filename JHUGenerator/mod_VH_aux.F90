!--YaofuZhou-----------------------------------------
module ModVHaux
  use ModParameters
  use ModMisc
  implicit none

  public :: PROPAGATOR
  public :: spinoru2
  public :: FFP
  public :: FFS

contains




!FFP.A
!VERSION 20130522

!returns i.Psi~(p1,s1).gamma5.Psi(p2,s2) for massless states
complex(8) function FFP(pdg_code1, p1, h1, pdg_code2, p2, h2)

  implicit none
  real(8), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

  real(8) p1(4), p2(4), h1, h2
  integer pdg_code1, pdg_code2
  real(8) sqrt_pp1Dpp2

  if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).gt.0d0)then
    FFP=0d0

  else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. (   dabs( p2(1)+p2(4) ).lt.epsilon ) )then
    FFP=0d0
  else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then
    FFP=-dsqrt(2d0*p1(1)*(p2(1)+p2(4)))
  else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then
    FFP= dsqrt(2d0*p2(1)*(p1(1)+p1(4)))
  else
    sqrt_pp1Dpp2 = dsqrt((p1(1)+p1(4))/(p2(1)+p2(4)))
    FFP=(p2(2)-(0d0,1d0)*p2(3))*sqrt_pp1Dpp2- (p1(2)-(0d0,1d0)*p1(3))/sqrt_pp1Dpp2
  endif

  FFP=FFP*(0d0,-1d0)

  if( (dble(pdg_code1)*h1) .lt. 0d0)then
    FFP=-dconjg(FFP)
  endif

  return
END function FFP



!FFS.F
!VERSION 20130522

!returns Psi~(p1,s1).Psi(p2,s2) for massless states
complex(8) function FFS(pdg_code1, p1, h1, pdg_code2, p2, h2)

  implicit none
  real(8), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision
  real(8) p1(4), p2(4), h1, h2
  integer pdg_code1, pdg_code2
  real(8) sqrt_pp1Dpp2

  FFS = (0d0,0d0)

  if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).gt.0d0)then
    FFS=0d0
  else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. (   dabs( p2(1)+p2(4) ).lt.epsilon ) )then
    FFS=0d0
  else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then
    FFS=-dsqrt(2d0*p1(1)*(p2(1)+p2(4)))
  else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then
    FFS= dsqrt(2d0*p2(1)*(p1(1)+p1(4)))
  else
    sqrt_pp1Dpp2 = dsqrt((p1(1)+p1(4))/(p2(1)+p2(4)))
    FFS=(p2(2)-(0d0,1d0)*p2(3))*sqrt_pp1Dpp2- (p1(2)-(0d0,1d0)*p1(3))/sqrt_pp1Dpp2
  endif

  if( (dble(pdg_code1)*h1) .lt. 0d0)then
    FFS=-dconjg(FFS)
  endif

  return
END function FFS

!PROPAGATOR.F
!VERSION 20130522
!
!PROPAGATOR() returns the generi!complex-valued propagator
!without tensor structure (numerator), given mass, invariant mass
!and width.
complex(8) function PROPAGATOR(inv_mass2, mass, width)
  implicit none

  real(8) inv_mass2, mass, width

!assuming auto-conversion. works with gfortran
  PROPAGATOR = -ci / ( inv_mass2 - mass**2 + ci*mass*width )
!     print *, PROPAGATOR

  return
END function PROPAGATOR



  !-- generic functions below
  !- MCFM spinors in non-MCFM momentum convention
subroutine spinoru2(n,p,za,zb,s)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: p(4,n)
  complex(8), intent(out) :: za(n,n), zb(n,n)
  real(8), intent(out) :: s(n,n)
  integer :: i,j
  complex(8) :: c23(n), f(n)
  real(8) :: rt(n)

  !---if one of the vectors happens to be zero this routine fails.
  do j=1,N
    za(j,j)=czero
    zb(j,j)=za(j,j)
    !-----positive energy case
    if (p(1,j) .gt. zero) then
        rt(j)=sqrt(abs(p(2,j)+p(1,j)))
        c23(j)=dcmplx(p(4,j),-p(3,j))
        f(j)=(one,zero)
    else
    !-----negative energy case
        rt(j)=sqrt(abs(-p(1,j)-p(2,j)))
        c23(j)=dcmplx(-p(4,j),p(3,j))
        f(j)=ci
    endif
  enddo

  do i=2,N
  do j=1,i-1
        s(i,j)=two*scr(p(:,i),p(:,j))
        za(i,j)=f(i)*f(j)  * ( c23(i)*dcmplx(rt(j)/(rt(i)+1d-16))-c23(j)*dcmplx(rt(i)/(rt(j)+1d-16)) )
        if (abs(s(i,j)).lt.1d-5) then
          zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
        else
          zb(i,j)=-dcmplx(s(i,j))/(za(i,j)+1d-16)
        endif
        za(j,i)=-za(i,j)
        zb(j,i)=-zb(i,j)
        s(j,i)=s(i,j)
    enddo
  enddo

  return

end subroutine spinoru2




end module ModVHaux
!!--YaofuZhou-----------------------------------------
