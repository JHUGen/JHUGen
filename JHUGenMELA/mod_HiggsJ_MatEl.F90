! matrix element squared of p(1) p(2) > h(3) j(4)
module ModHiggsJ
  implicit none
  private

  public :: EvalAmp_HJ

!-- general definitions, to be merged with Markus final structure
  integer, parameter  :: dp = selected_real_kind(15)

  include './variables.F90'

  real(8),  parameter :: pi =3.141592653589793238462643383279502884197d0  
  
  real(dp), parameter :: gs = sqrt(alphas*4.0_dp*pi)

  real(dp), parameter :: xn = 3.0_dp
  real(dp), parameter :: Ca = 3.0_dp
  real(dp), parameter :: Cf = 4.0_dp/3.0_dp

  real(dp), parameter :: avegg = 1.0_dp/4.0_dp/64.0_dp
  real(dp), parameter :: aveqg = 1.0_dp/4.0_dp/24.0_dp
  real(dp), parameter :: aveqq = 1.0_dp/4.0_dp/9.0_dp

contains

subroutine EvalAmp_HJ(p,res)
    real(dp), intent(in) :: p(4,4)
    complex(dp) :: heftcoupl
    real(dp), intent(out) :: res(-5:5,-5:5)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: restmp!, restmpid
    integer :: i,j

    heftcoupl = gs*alphas/(6.0_dp * pi * vev)
    res=0d0

    call spinoru(3,(/-p(:,1),-p(:,2),p(:,4)/),za,zb,sprod)

    call me2_ggg_tree(1,2,3,sprod,res(0,0))
    res(0,0) = res(0,0)*avegg
    
    do i=1,5
      call me2_qbqg_tree(2,3,1,sprod,res(0,i))
      res(0,i) = -res(0,i)*aveqg
      !print *,res(0,i)
    enddo

    do i=1,5
      call me2_qbqg_tree(1,3,2,sprod,res(i,0))
      res(i,0) = -res(i,0)*aveqg
      !print *,res(i,0)
    enddo

    do i=1,5
      call me2_qbqg_tree(3,2,1,sprod,res(0,-i))
      res(0,-i) = -res(0,-i)*aveqg
      !print *,res(0,-i)
    enddo

    do i=1,5
      call me2_qbqg_tree(1,3,2,sprod,res(-i,0))
      res(-i,0) = -res(-i,0)*aveqg
    !if(res(-i,0).le.0d0) then
    !print *,res(-i,0)
    !endif
    enddo

    do i=1,5
      call me2_qbqg_tree(1,2,3,sprod,res(-i,i))
      res(-i,i) = res(-i,i) * aveqq
!      res(i,-i) = res(i,-i) + restmp * aveqq
!      res(0, i) = res(0, i) + restmp * aveqg
!      res(0,-i) = res(0,-i) + restmp * aveqg
!      res(i, 0) = res(i, 0) + restmp * aveqg
!      res(-i,0) = res(-i,0) + restmp * aveqg
    enddo

    do i=1,5
      call me2_qbqg_tree(2,1,3,sprod,res(i,-i))
      res(i,-i) = res(i,-i)*aveqq
    enddo

    res = res * (heftcoupl**2)
!if(res(0,0).gt.10000d0)then
    !do i = -5,5
    !do j = -5,5
    !print *, res(i,j),i,j
    !enddo
    !enddo
!    print *, dsqrt(scr(p(:,1),p(:,1))), dsqrt(scr(p(:,2),p(:,2))), dsqrt(scr(p(:,3),p(:,3))), dsqrt(scr(p(:,4),p(:,4)))
!    print*, p
!endif

end subroutine EvalAmp_HJ





subroutine me2_ggg_tree(j1,j2,j3,sprod,me2)
   integer, intent(in) :: j1,j2,j3
   real(dp), intent(in) :: sprod(3,3)
   real(dp), intent(out) :: me2
   real(dp) :: s12, s13, s23, qsq
   real(dp), parameter :: col_gg = 4.0d0 * CA**2 * CF

   me2 = 0d0

   s12 = sprod(j1,j2)
   s13 = sprod(j1,j3)
   s23 = sprod(j2,j3)

   qsq = s12 + s13 + s23

   me2 = (s12**4 + s13**4 + s23**4 + qsq**4)/s12/s13/s23

   me2 = me2 * 2d0 * col_gg

   return

end subroutine me2_ggg_tree


!-- 0 -> qb(p1) q(p2) g(p3) H
subroutine me2_qbqg_tree(j1,j2,j3,sprod,me2)
   integer, intent(in) :: j1, j2, j3
   real(dp), intent(in) :: sprod(3,3)
   real(dp), intent(out) :: me2
   real(dp) :: s12, s13,s23
   real(dp), parameter :: col_qg = 2d0 * xn * CF

   me2 = 0d0

   s12 = sprod(j1,j2)
   s13 = sprod(j1,j3)
   s23 = sprod(j2,j3)

   me2 = 2d0 * (s13**2 + s23**2)/s12 * col_qg
!if(me2.le.0d0)then
!print *, s12,s13,s23
!endif
   return

end subroutine me2_qbqg_tree


  !-------------------------------------------------------------------------
  !-- generic functions below
  !- MCFM spinors
  subroutine spinoru(n,p,za,zb,s)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,n)
    complex(dp), intent(out) :: za(n,n), zb(n,n)
    real(dp), intent(out) :: s(n,n)
    integer :: i,j
    complex(dp) :: c23(n), f(n)
    real(dp) :: rt(n)
      
    !---if one of the vectors happens to be zero this routine fails.
    do j=1,N
       za(j,j)=(0d0,0d0)
       zb(j,j)=za(j,j)

       !-----positive energy case
       if (p(1,j) .gt. 0d0) then
          rt(j)=sqrt(p(2,j)+p(1,j))
          c23(j)=cmplx(p(4,j),-p(3,j),kind=dp)
          f(j)=(1d0,0d0)
       else
       !-----negative energy case
          rt(j)=sqrt(-p(1,j)-p(2,j))
          c23(j)=cmplx(-p(4,j),p(3,j),kind=dp)
          f(j)=(0d0,1d0)
       endif
    enddo

    do i=2,N
  
     do j=1,i-1
          s(i,j)=2d0*scr(p(:,i),p(:,j))
          za(i,j)=f(i)*f(j) &
               *(c23(i)*cmplx(rt(j)/rt(i),kind=dp)-c23(j)*cmplx(rt(i)/rt(j),kind=dp))
          
          if (abs(s(i,j)).lt.1d-5) then
             zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
          else
             zb(i,j)=-cmplx(s(i,j),kind=dp)/za(i,j)
          endif
          
          za(j,i)=-za(i,j)
          zb(j,i)=-zb(i,j)
          s(j,i)=s(i,j)
          
       enddo

    enddo

    return
    
  end subroutine spinoru


  function scr(p1,p2) !coincide with MinkowskyProduct() in mod_Misc
    real(dp), intent(in) :: p1(4), p2(4)
    real(dp) :: scr
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
  end function scr




end module ModHiggsJ