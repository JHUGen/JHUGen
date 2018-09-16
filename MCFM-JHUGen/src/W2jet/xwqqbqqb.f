      function xwqqbqqb(i1,i2,i3,i4,n1,n2)
      implicit none
      include 'types.f'
      real(dp):: xwqqbqqb

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'basic.f'
      include 'ckm.f'
      include 'ckm1.f'

      real(dp):: fac,xma,xmb,xmc,xmd,xme
      integer:: n1,n2,nj,nk,i1,i2,i3,i4


      complex(dp):: mlll1,mlll2,mlrl1,mlrl2
      complex(dp):: mrll1,mrll2



      if (n1 .ne. -n2) then
Case 3 Distinct quarks- charged coupling
C*************************************************************

C++++ subcase a
c     qi(-p1)+qbarj(-p2)-->qk(p3)+qbarj(p4)+W+/-  (k .ne. i,j)

      fac=0._dp
      do nk=-nf,nf
      fac=fac+(fl*gl(n1,nk)*VV(n1,nk))**2
c      write(6,*) 'n1,nk,fac',n1,nk,fac
c      write(6,*) 'gl(n1,nk)',gl(n1,nk)
c      write(6,*) 'VV(n1,nk)',VV(n1,nk)
      enddo
      fac=fac-(fl*gl(n1,-n2)*VV(n1,n2))**2
c      write(6,*) 'n1',n1
c      write(6,*) 'n2',n2

c      left-left matrix element
      mlll1=+Lla(i1,i2,i3,i4)
c      left-right matrix element
      mlrl1=+Lla(i1,i4,i3,i2)

      xma=fac*aveqq*Von4*(abs(mlll1)**2+abs(mlrl1)**2)

C++++ subcase b
c     qi(-p1)+qbarj(-p2)-->qi(p3)+qbark(p4)+W+/-  (k .ne. i,j)

      fac=0._dp
      do nk=-nf,nf
      fac=fac+(fl*gl(n2,nk)*VV(n2,nk))**2
      enddo
      fac=fac-(fl*gl(n2,-n1)*VV(n2,-n1))**2
c      left-left matrix element
      mlll1=+Lla(i2,i1,i4,i3)
c      right-left matrix element
      mrll1=+Lla(i2,i3,i4,i1)

      xmb=fac*aveqq*Von4*(abs(mlll1)**2+abs(mrll1)**2)


c++++ subcase c

c     qi(-p1)+qbarj(-p2)-->qj(p3)+qbarj(p4)+W+/-  (k .ne. j)


      fac=(fl*gl(n1,n2)*VV(n1,n2))**2

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
c     qi(-p1)+qbarj(-p2)-->qi(p3)+qbari(p4)+W+/-  (k == i)

      fac=(fl*gl(n2,n1)*VV(n2,n1))**2

c      left-left matrix element
      mlll1=+Lla(i2,i1,i4,i3)
      mlll2=-Lla(i1,i2,i4,i3)
c      right-left matrix element
      mrll1=+Lla(i2,i3,i4,i1)
      mrll2=-Lla(i1,i3,i4,i2)

      xmd=fac*half*aveqq*Von4*(
     & +abs(mlll1)**2+abs(mlll2)**2-two/XN*real(mlll1*conjg(mlll2))
     & +abs(mrll1)**2+abs(mrll2)**2)

c++++ subcase e

c     qi(-p1)+qbarj(-p2)-->qk(p3)+qbark(p4)+W+/-  (k .ne. i,j)

c      left-left matrix element
      fac=(fl*gl(n1,n2)*VV(n1,n2))**2*(nf-2)
      mlll2=-Lla(i1,i2,i4,i3)
c      right-left matrix element
      mlrl2=-Lla(i1,i3,i4,i2)
      xme=fac*half*aveqq*Von4*(
     & +abs(mlll2)**2+abs(mlrl2)**2)

      xwqqbqqb=xma+xmb+xmc+xm.e+_dpxme

      return

C************************************************************
      elseif (n1 == -n2) then

Case 4 Identical qi-qbari - charged coupling
c++++ subcase a
c     qi(-p1)+qbari(-p2)-->qj(p3)+qbari(p4)+W+/-  (k .ne. i)

      fac=0._dp
      do nk=-nf,nf
      fac=fac+(fl*gl(n1,nk)*VV(n1,nk))**2
      enddo

c----left-left matrix element
      mlll1=+Lla(i1,i2,i3,i4)
      mlll2=-Lla(i2,i1,i3,i4)

c----left-right matrix element
      mlrl1=+Lla(i1,i4,i3,i2)
      mlrl2=-Lla(i2,i4,i3,i1)

      xma=fac*aveqq*Von4*(
     & +abs(mlll1)**2+abs(mlll2)**2
     & -two/XN*real(mlll1*conjg(mlll2))
     & +abs(mlrl1)**2+abs(mlrl2)**2)

c++++ subcase b
c     qi(-p1)+qbari(-p2)-->qi(p3)+qbark(p4)+W+/-  (k .ne. i)

      fac=0._dp
      do nk=-nf,nf
      fac=fac+(fl*gl(n2,nk)*VV(n2,nk))**2
      enddo

c----left-left matrix element
      mlll1=+Lla(i2,i1,i4,i3)
      mlll2=-Lla(i2,i1,i3,i4)

c----right-left matrix element
      mrll1=+Lla(i2,i3,i4,i1)
      mrll2=-Lla(i2,i4,i3,i1)


      xmb=fac*aveqq*Von4*(
     & +abs(mlll1)**2+abs(mlll2)**2
     & -two/XN*real(mlll1*conjg(mlll2))
     & +abs(mrll1)**2+abs(mrll2)**2)


c++++ subcase c
c     qi(-p1)+qbari(-p2)-->qj(p3)+qbark(p4)+W+/-  (j,k .ne. i)

      fac=0._dp

      do nj=-nf,nf
      do nk=-nf,nf
      if ((nj == n1) .or. (nk == n2)) then
      continue
      else
      fac=fac+(fl*gl(nj,nk)*VV(nj,nk))**2
      endif
      enddo
      enddo

c----left-left matrix element
      mlll2=+Lla(i2,i4,i3,i1)

c----left-right matrix element
      mlrl2=+Lla(i2,i1,i3,i4)

      xmc=fac*aveqq*Von4*(abs(mlll2)**2+abs(mlrl2)**2)

      xwqqbqqb=xma+xmb+xmc

      endif

      return
      end

