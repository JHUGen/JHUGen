      double complex function Ftriangle(x)
      implicit none 
      include 'constants.f'
      double precision x,y
      double complex arg
      y=1d0-4d0*x
      if (y .gt. 0d0) then
      arg=dcmplx((1d0+sqrt(y))/(1d0-sqrt(y)))
      Ftriangle=+dcmplx(0.5d0)*(log(arg)-impi)**2
      elseif (y .le. 0d0) then
      Ftriangle=-dcmplx(2d0)*dcmplx(asin(half/sqrt(x)))**2
      endif
      return 
      end 

