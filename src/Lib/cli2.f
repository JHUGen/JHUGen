      double complex function cli2(x)
c--complex dilogarithm (spence-function)
      implicit double precision (a-h,o-z)
      double complex x,y,li2taylor
      double precision zeta2,zeta3
      common/const/zeta2,zeta3
      logical first
      data first/.true./
      save first

      if (first) then
      first=.false.
      call bernini
      endif

      zero=1.d-8
      xr=dble(x)
      xi=dimag(x)
      r2=xr*xr+xi*xi
      cli2=dcmplx(0d0,0d0)
      if(r2.le.zero)then
        cli2=x+x**2/4.d0
        return
      endif
      rr=xr/r2
      if ((r2.eq.1.d0) .and. (xi.eq.0.d0)) then
        if (xr.eq.1.d0) then
          cli2=dcmplx(zeta2)
        else
          cli2=-dcmplx(zeta2/2.d0)
        endif
        return
      elseif ((r2.gt.1.d0) .and. (rr.gt.0.5d0)) then
        y=(x-1.d0)/x
        cli2=li2taylor(y)+zeta2-cdlog(x)*cdlog(1.d0-x)+0.5d0*cdlog(x)**2
        return
      elseif ((r2.gt.1.d0) .and. (rr.le.0.5d0))then
        y=1.d0/x
        cli2=-li2taylor(y)-zeta2-0.5d0*cdlog(-x)**2
        return
      elseif ((r2.le.1.d0) .and. (xr.gt.0.5d0)) then
        y=1.d0-x
        cli2=-li2taylor(y)+zeta2-cdlog(x)*cdlog(1.d0-x)
       return
      elseif ((r2.le.1.d0) .and. (xr.le.0.5d0)) then
        y=x
        cli2=li2taylor(y)
        return
      endif
      end
 
      double complex function li2taylor(x)
c--taylor-expansion for complex dilogarithm (spence-function)
      implicit double precision (a-h,o-z)
      parameter(nber=18)
      double precision b2(nber) 
      double complex x,z
      common/bernoulli/b2

      n=nber-1
      z=-cdlog(1.d0-x)
      li2taylor=b2(nber)
      do 111 i=n,1,-1
        li2taylor=z*li2taylor+b2(i)
111   continue
      li2taylor=z**2*li2taylor+z
      return
      end
 
      double precision function facult(n)
c--double precision version of faculty
      implicit double precision (a-h,o-z)
      facult=1.d0
      if(n.eq.0)return
      do 999 i=1,n
        facult=facult*dfloat(i)
999   continue
      return
      end
 
      subroutine bernini
c--initialization of coefficients for polylogarithms
      implicit none
      include 'constants.f'
      integer nber,i
      parameter(nber=18)
      double precision b(nber),b2(nber),zeta2,zeta3,facult
      common/bernoulli/b2
      common/const/zeta2,zeta3
 
 
      b(1)=-1.d0/2.d0
      b(2)=1.d0/6.d0
      b(3)=0.d0
      b(4)=-1.d0/30.d0
      b(5)=0.d0
      b(6)=1.d0/42.d0
      b(7)=0.d0
      b(8)=-1.d0/30.d0
      b(9)=0.d0
      b(10)=5.d0/66.d0
      b(11)=0.d0
      b(12)=-691.d0/2730.d0
      b(13)=0.d0
      b(14)=7.d0/6.d0
      b(15)=0.d0
      b(16)=-3617.d0/510.d0
      b(17)=0.d0
      b(18)=43867.d0/798.d0
      zeta2=pi**2/6.d0
      zeta3=1.202056903159594d0
 
      do 995 i=1,nber
        b2(i)=b(i)/facult(i+1)
995   continue
 
      return
      end
 
