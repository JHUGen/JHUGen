      subroutine xwqqggg(j1,j2,j3,j4,j5,j6,j7,mqqb)
      implicit none
      include 'types.f'

************************************************************************
*     Author J.M.Campbell, February 2001                               *
*     Returns the amplitudes squared for the process                   *
*     0 ---> q(p1)+g(p2)+g(p3)+g(p4)+qbar(p5)+l(p6)+a(p7)              *
*     mqqb is the matrix element for left-handed quark and lepton      *
*     Averaging over 2 initial state gluons is assumed, and            *
*     no final state average is included                               *
*                                                                      *
*     Specifying colourchoice = 1 --> leading colour only              *
*     Specifying colourchoice = 2 --> sub-leading colour only          *
*     Specifying colourchoice = 3 --> sub-sub-leading colour only      *
*     Specifying colourchoice = 0 --> TOTAL                            *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'lc.f'
      integer:: i2(6),i3(6),i4(6),j,lh,h2,h3,h4,hq,h(7)
      integer:: j1,j2,j3,j4,j5,j6,j7
      real(dp):: mqqb,m1,m2,m0,fac
      complex(dp):: tempm0,m(6),amp_qqggg
      mqqb=0._dp

      fac=avegg*gsq**3*esq**2*xn**3*cf*8._dp
c--- extra factor of 8 due to colour matrix normalization (rt2**6)

      i2(1)=j2
      i3(1)=j3
      i4(1)=j4
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3
      i2(3)=j4
      i3(3)=j2
      i4(3)=j3
      i2(4)=j3
      i3(4)=j4
      i4(4)=j2
      i2(5)=j3
      i3(5)=j2
      i4(5)=j4
      i2(6)=j4
      i3(6)=j3
      i4(6)=j2

c--- left-left helicity configuration only for fermions
      hq=1
      lh=1

      do h2=1,2
      do h3=1,2
      do h4=1,2

        h(j2)=h2
        h(j3)=h3
        h(j4)=h4
        tempm0=czip
        m2=zip

        do j=1,6
          m(j)=amp_qqggg(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                                       i4(j),h(i4(j)),j5,lh,j6,j7)
          tempm0=tempm0+m(j)
          m2=m2+abs(m(j))**2
        enddo

        if ((colourchoice == 1) .or. (colourchoice == 0)) then
          mqqb=mqqb+fac*m2
        endif
        if ((colourchoice == 2) .or. (colourchoice == 0)) then
c--- here we have (2,3,4)+(2,4,3)+(4,2,3) [4 is photon-like]
c---         plus (3,4,2)+(3,2,4)+(2,3,4) [2 is photon-like]
c---         plus (4,2,3)+(4,3,2)+(3,4,2) [3 is photon-like]
c--- (plus perms)
          m1=abs(m(1)+m(2)+m(3))**2
     &      +abs(m(4)+m(5)+m(1))**2
     &      +abs(m(3)+m(6)+m(4))**2
     &      +abs(m(5)+m(4)+m(6))**2
     &      +abs(m(2)+m(1)+m(5))**2
     &      +abs(m(6)+m(3)+m(2))**2
          mqqb=mqqb+fac*(-m1/xnsq)
        endif
        if ((colourchoice == 3) .or. (colourchoice == 0)) then
          m0=abs(tempm0)**2
          mqqb=mqqb+fac*m0*(xnsq+1._dp)/xnsq**2
        endif

      enddo
      enddo
      enddo

      return
      end
