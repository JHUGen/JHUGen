      function lnrat(x,y)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: lnrat
************************************************************************
*     Author: R.K. Ellis                                               *
*     August, 1998.                                                    *
c     lnrat(x,y)=log(x-i*ep)-log(y-i*ep)                               *
c     this function is har.e-_dpwired for sign of epsilon we must adjust   *
c     sign of x and y to get the right sign for epsilon                *
************************************************************************
      include 'constants.f'
       real(dp):: x,y,htheta
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+half*sign(one,x)
!---  end statement function
      lnrat=cplx1(log(abs(x/y)))
     & -impi*(htheta(-x)-htheta(-y))
      return
      end

