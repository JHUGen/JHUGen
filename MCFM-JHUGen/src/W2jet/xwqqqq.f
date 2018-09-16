      function xwqqqq(i1,i2,i3,i4,n1,n2)
      implicit none
      include 'types.f'
      real(dp):: xwqqqq

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'basic.f'
      include 'ckm.f'
      include 'ckm1.f'

      real(dp):: fac,xma,xmb,xmc,xmd
      integer:: n1,n2,nk,i1,i2,i3,i4

      complex(dp):: mlll1,mlll2,mlrl1,mlrl2
      complex(dp):: mrll1,mrll2


      if (n1 .ne. n2) then
Case 3 Distinct quarks- charged coupling
C*************************************************************

C++++ subcase a
c     qi(-p1)+qj(-p2)-->qk(p3)+qj(p4)+W+/-  (k .ne. i,j)

      fac=0._dp
      do nk=-nf,nf
      if ((nk == n1) .or. (nk == -n2)) then
      continue
      else
      fac=fac+(fl*gl(n1,nk)*VV(n1,nk))**2
      endif
      enddo

c      left-left matrix element
      mlll1=+Lla(i1,i2,i3,i4)
c      left-right matrix element
      mlrl1=+Lla(i1,i4,i3,i2)

      xma=fac*aveqq*Von4*(abs(mlll1)**2+abs(mlrl1)**2)

C++++ subcase b
c     qi(-p1)+qj(-p2)-->qi(p3)+qk(p4)+W+/-  (k .ne. i,j)

      fac=0._dp
      do nk=-nf,nf
      if ((nk == -n1) .or. (nk == n2)) then
      continue
      else
      fac=fac+(fl*gl(n2,nk)*VV(n2,nk))**2
      endif
      enddo

c      left-left matrix element
      mlll1=+Lla(i2,i1,i4,i3)
c      right-left matrix element
      mrll1=+Lla(i2,i3,i4,i1)

      xmb=fac*aveqq*Von4*(abs(mlll1)**2+abs(mrll1)**2)


c++++ subcase c

c     qi(-p1)+qj(-p2)-->qj(p3)+qj(p4)+W+/-  (k == j)


      fac=(fl*gl(n1,-n2)*VV(n1,-n2))**2

c      left-left matrix element
      mlll1=+Lla(i1,i2,i3,i4)
      mlll2=-Lla(i1,i2,i4,i3)
c      left-right matrix element
      mlrl1=+Lla(i1,i4,i3,i2)
      mlrl2=-Lla(i1,i3,i4,i2)

      xmc=fac*half*aveqq*Von4*(
     & +abs(mlll1)**2+abs(mlll2)**2-two/XN*real(mlll1*conjg(mlll2))
     & +abs(mlrl1)**2+abs(mlrl2)**2)

c++++ subcase d
c     qi(-p1)+qj(-p2)-->qi(p3)+qi(p4)+W+/-  (k == i)

      fac=(fl*gl(n2,-n1)*VV(n2,-n1))**2

c      left-left matrix element
      mlll1=+Lla(i2,i1,i4,i3)
      mlll2=-Lla(i2,i1,i3,i4)
c      right-left matrix element
      mrll1=+Lla(i2,i3,i4,i1)
      mrll2=-Lla(i2,i4,i3,i1)

      xmd=fac*half*aveqq*Von4*(
     & +abs(mlll1)**2+abs(mlll2)**2-two/XN*real(mlll1*conjg(mlll2))
     & +abs(mrll1)**2+abs(mrll2)**2)

      xwqqqq=xma+xmb+xmc+xmd

      return

C************************************************************
      elseif (n1 == n2) then
Case 4 Identical quarks - charged coupling
c     qi(-p1)+qi(-p2)-->qk(p3)+qi(p4)+W+/-  (k .ne. i)

      fac=0._dp
      do nk=-nf,nf
      fac=fac+(fl*gl(n1,nk)*VV(n1,nk))**2
      enddo
c      write(6,*) 'fl',fl
c      write(6,*) 'n1,nk,VV(n1,nk)',n1,nk,VV(n1,nk)
c      write(6,*) 'gl(n1,nk)',n1,nk,gl(n1,nk)
c      write(6,*) 'n1',n1
c      write(6,*) 'n2',n2
c      write(6,*) 'fac',fac

c      write(6,*) 'entering here'

c----left-left matrix element
      mlll1=+Lla(i1,i2,i3,i4)
      mlll2=-Lla(i2,i1,i3,i4)

c--- new expression, testing
c      mlll2=-Lla(i1,i2,i4,i3)
c--- new expression, testing

c----left-right matrix element
      mlrl1=+Lla(i1,i4,i3,i2)

c----right-left matrix element
      mrll2=-Lla(i2,i4,i3,i1)

c--- new expression, testing
c      mrll2=-Lla(i1,i3,i4,i2)
c--- new expression, testing

      xwqqqq=fac*aveqq*Von4*(
     & +abs(mlll1)**2+abs(mlll2)**2
     & -two/XN*real(mlll1*conjg(mlll2))
     & +abs(mlrl1)**2
     & +abs(mrll2)**2)

      endif
      return
      end






