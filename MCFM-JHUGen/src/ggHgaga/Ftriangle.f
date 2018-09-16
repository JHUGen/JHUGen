      function Ftriangle(x)
      implicit none
      include 'types.f'
      complex(dp):: Ftriangle
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: x,y
      complex(dp):: arg
      y=1._dp-4._dp*x
      if (y > 0._dp) then
      arg=cplx1((1._dp+sqrt(y))/(1._dp-sqrt(y)))
      Ftriangle=+chalf*(log(arg)-impi)**2
      elseif (y <= 0._dp) then
      Ftriangle=-ctwo*cplx1(asin(half/sqrt(x)))**2
      endif
      return 
      end 

