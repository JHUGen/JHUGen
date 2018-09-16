      function li2(x)
      implicit none
      include 'types.f'
      complex(dp):: li2
c--complex dilogarithm (spence-function)
      implicit real(dp) (a-h,o-z)
      include 'first.f'
      complex(dp):: x,y,cli2
      real(dp):: zeta2,zeta3
      common/const/zeta2,zeta3

      if (first) then
      first=.false.
      call bernini
      endif

      tiny=1.e-8_dp
      xr=real(x)
      xi=aimag(x)
      r2=xr*xr+xi*xi
      li2=cplx2(0._dp,0._dp)
      if(r2<=tiny)then
        li2=x+x**2/4._dp
        return
      endif
      rr=xr/r2
      if ((r2==1._dp) .and. (xi==0._dp)) then
        if (xr==1._dp) then
          li2=cplx2(zeta2)
        else
          li2=-cplx2(zeta2/2._dp)
        endif
        return
      elseif ((r2>1._dp) .and. (rr>0.5_dp)) then
        y=(x-1._dp)/x
        li2=cli2(y)+zeta2-log(x)*log(1._dp-x)+0.5_dp*log(x)**2
        return
      elseif ((r2>1._dp) .and. (rr<=0.5_dp))then
        y=1._dp/x
        li2=-cli2(y)-zeta2-0.5_dp*log(-x)**2
        return
      elseif ((r2<=1._dp) .and. (xr>0.5_dp)) then
        y=1._dp-x
        li2=-cli2(y)+zeta2-log(x)*log(1._dp-x)
       return
      elseif ((r2<=1._dp) .and. (xr<=0.5_dp)) then
        y=x
        li2=cli2(y)
        return
      endif
      end

      function cli2(x)
      implicit none
      include 'types.f'
      complex(dp):: cli2
c--taylor-expansion for complex dilogarithm (spence-function)
      implicit real(dp) (a-h,o-z)
      parameter(nber=18)
      real(dp):: b2(nber)
      complex(dp):: x,z
      common/bernoulli/b2

      n=nber-1
      z=-log(1._dp-x)
      cli2=b2(nber)
      do 111 i=n,1,-1
        cli2=z*cli2+b2(i)
111   continue
      cli2=z**2*cli2+z
      return
      end

      function facult(n)
      implicit none
      include 'types.f'
      real(dp):: facult
c--real(dp):: version of faculty
      implicit real(dp) (a-h,o-z)
      facult=1._dp
      if(n==0)return
      do 999 i=1,n
        facult=facult*real(i,dp)
999   continue
      return
      end

      subroutine bernini
      implicit none
      include 'types.f'
c--initialization of coefficients for polylogarithms

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
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

