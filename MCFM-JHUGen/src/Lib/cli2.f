      function cli2(x)
      implicit none
      include 'types.f'
      complex(dp):: cli2
c--complex dilogarithm (spence-function)
      include 'constants.f'
      complex(dp):: x,y,li2taylor
      real(dp):: zeta2,zeta3,tiny,rr,r2,xr,xi
      common/const/zeta2,zeta3
      logical,save:: first=.true.
      include 'cplx.h'

      if (first) then
      first=.false.
      call bernini
      endif

      tiny=1.e-8_dp
      xr=real(x)
      xi=aimag(x)
      r2=xr*xr+xi*xi
      cli2=cplx2(zip,zip)
      if(r2<=tiny)then
        cli2=x+x**2/4._dp
        return
      endif
      rr=xr/r2
      if ((r2==one) .and. (xi==zip)) then
        if (xr==one) then
          cli2=cplx1(zeta2)
        else
          cli2=-cplx1(zeta2/two)
        endif
        return
      elseif ((r2>one) .and. (rr>half)) then
        y=(x-one)/x
        cli2=li2taylor(y)+zeta2-log(x)*log(one-x)+half*log(x)**2
        return
      elseif ((r2>one) .and. (rr<=half))then
        y=one/x
        cli2=-li2taylor(y)-zeta2-half*log(-x)**2
        return
      elseif ((r2<=one) .and. (xr>half)) then
        y=one-x
        cli2=-li2taylor(y)+zeta2-log(x)*log(one-x)
       return
      elseif ((r2<=one) .and. (xr<=half)) then
        y=x
        cli2=li2taylor(y)
        return
      endif
      end

      function li2taylor(x)

      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: li2taylor
c--taylor-expansion for complex dilogarithm (spence-function)
      integer:: nber,i,n
      parameter(nber=18)
      real(dp):: b2(nber)
      complex(dp):: x,z
      common/bernoulli/b2

      n=nber-1
      z=-log(one-x)
      li2taylor=b2(nber)
      do 111 i=n,1,-1
      li2taylor=z*li2taylor+b2(i)
111   continue
      li2taylor=z**2*li2taylor+z
      return
      end

      function facult(n)
c--real(dp):: version of faculty
      implicit none
      include 'types.f'
      real(dp):: facult
      integer:: i,n
      include 'constants.f'
      facult=one
      if(n==0)return
      do 999 i=1,n
        facult=facult*real(i,dp)
999   continue
      return
      end

      subroutine bernini

c--initialization of coefficients for polylogarithms
      implicit none
      include 'types.f'
      include 'constants.f'
      integer:: nber,i
      parameter(nber=18)
      real(dp):: b(nber),b2(nber),zeta2,zeta3,facult
      common/bernoulli/b2
      common/const/zeta2,zeta3


      b(1)=-1._dp/2._dp
      b(2)=1._dp/6._dp
      b(3)=0._dp
      b(4)=-1._dp/30._dp
      b(5)=0._dp
      b(6)=1._dp/42._dp
      b(7)=0._dp
      b(8)=-1._dp/30._dp
      b(9)=0._dp
      b(10)=5._dp/66._dp
      b(11)=0._dp
      b(12)=-691._dp/2730._dp
      b(13)=0._dp
      b(14)=7._dp/6._dp
      b(15)=0._dp
      b(16)=-3617._dp/510._dp
      b(17)=0._dp
      b(18)=43867._dp/798._dp
      zeta2=pi**2/6._dp
      zeta3=1.202056903159594_dp

      do 995 i=1,nber
        b2(i)=b(i)/facult(i+1)
995   continue

      return
      end

