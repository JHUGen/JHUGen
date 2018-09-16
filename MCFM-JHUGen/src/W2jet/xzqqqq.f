      function xzqqqq(i,n1,n2)
      implicit none
      include 'types.f'
      real(dp):: xzqqqq

c*******************************************************************
c     the matrix elements of the
C     helicity amplitudes for the QCD process
c     q(-q1)+q(-q2)+a(-q6) --> q(q3)+q(q4)+l(q5)
c     all squared
c     multiplied by (((a+l)^2-M**2)^2+M^2*Gam^2)/(a+l)^4/g^4
c     averaged over initial colour or spin
c*******************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'ckm1.f'
      complex(dp):: flgl1,flgl2,frgl1,frgl2,flgr1,flgr2,frgr1,frgr2


      integer:: n1,n2

      complex(dp):: mlll1,mlll2,mlrl1,mlrl2,mllr1,mllr2,mlrr1,mlrr2
      complex(dp):: mrll1,mrll2,mrrl1,mrrl2,mrlr1,mrlr2,mrrr1,mrrr2
      complex(dp):: bwf
      common/bwf/bwf
!$omp threadprivate(/bwf/)
      complex(dp):: Lla(4,4,4,4),Lal(4,4,4,4)
      complex(dp):: Rla(4,4,4,4),Ral(4,4,4,4)
      common/basic/Lla,Lal,Rla,Ral
!$omp threadprivate(/basic/)

      integer:: i(4)
      xzqqqq=0._dp

      flgl1=fl*gl(n1,n1)-bwf*e(n1,n1)
      flgl2=fl*gl(n2,n2)-bwf*e(n2,n2)
      frgl1=fr*gl(n1,n1)-bwf*e(n1,n1)
      frgl2=fr*gl(n2,n2)-bwf*e(n2,n2)

      flgr1=fl*gr(n1,n1)-bwf*e(n1,n1)
      flgr2=fl*gr(n2,n2)-bwf*e(n2,n2)
      frgr1=fr*gr(n1,n1)-bwf*e(n1,n1)
      frgr2=fr*gr(n2,n2)-bwf*e(n2,n2)


c----left-left matrix element
      mlll1=
     & +flgl1*Lla(i(1),i(2),i(3),i(4))+flgl2*Lla(i(2),i(1),i(4),i(3))
      mllr1=
     & +frgl1*Lal(i(1),i(2),i(3),i(4))+frgl2*Lal(i(2),i(1),i(4),i(3))

c----right-right matrix element
      mrrr1=
     & +frgr1*Rla(i(1),i(2),i(3),i(4))+frgr2*Rla(i(2),i(1),i(4),i(3))
      mrrl1=
     & +flgr1*Ral(i(1),i(2),i(3),i(4))+flgr2*Ral(i(2),i(1),i(4),i(3))

c----left-right matrix element
      mlrl1=
     & +flgl1*Lla(i(1),i(4),i(3),i(2))+flgr2*Ral(i(2),i(3),i(4),i(1))
      mlrr1=
     & +frgl1*Lal(i(1),i(4),i(3),i(2))+frgr2*Rla(i(2),i(3),i(4),i(1))

c----right-left matrix element
      mrll1=
     & +flgr1*Ral(i(1),i(4),i(3),i(2))+flgl2*Lla(i(2),i(3),i(4),i(1))
      mrlr1=
     & +frgr1*Rla(i(1),i(4),i(3),i(2))+frgl2*Lal(i(2),i(3),i(4),i(1))


C*************************************************************
Case 1 Non-identical quarks diagonal coupling                *
C     qi(-p1)+qj(-p2) --> qi(p3)+qj(p4) + Z0,Gamma           *
c*************************************************************

      if (n1 .ne. n2) then
      xzqqqq=aveqq*von4*(
     & +abs(mlll1)**2+abs(mllr1)**2
     & +abs(mrrl1)**2+abs(mrrr1)**2
     & +abs(mlrl1)**2+abs(mlrr1)**2
     & +abs(mrll1)**2+abs(mrlr1)**2)

      return

      elseif (n1 == n2) then

C************************************************************
Case 2 Identical quarks diagonal coupling                   *
C     qi(-p1)+qi(-p2) --> qi(p3)+qi(p4) + Z0,Gamma          *
C************************************************************


c----left-left matrix element
      mlll2=-flgl1*(Lla(i(1),i(2),i(4),i(3))+Lla(i(2),i(1),i(3),i(4)))
      mllr2=-frgl1*(Lal(i(1),i(2),i(4),i(3))+Lal(i(2),i(1),i(3),i(4)))

c----right-right matrix element
      mrrr2=-frgr1*(Rla(i(1),i(2),i(4),i(3))+Rla(i(2),i(1),i(3),i(4)))
      mrrl2=-flgr1*(Ral(i(1),i(2),i(4),i(3))+Ral(i(2),i(1),i(3),i(4)))

c----left-right matrix element
      mlrl2=
     & -flgl1*Lla(i(1),i(3),i(4),i(2))-flgr2*Ral(i(2),i(4),i(3),i(1))
      mlrr2=
     & -frgl1*Lal(i(1),i(3),i(4),i(2))-frgr2*Rla(i(2),i(4),i(3),i(1))

c----right-left matrix element
      mrll2=
     & -flgr1*Ral(i(1),i(3),i(4),i(2))-flgl2*Lla(i(2),i(4),i(3),i(1))
      mrlr2=
     & -frgr1*Rla(i(1),i(3),i(4),i(2))-frgl2*Lal(i(2),i(4),i(3),i(1))


      xzqqqq=half*aveqq*Von4*(
     & +abs(mrrr1)**2+abs(mrrr2)**2
     & -two/XN*real(mrrr1*conjg(mrrr2),dp)
     & +abs(mrlr1)**2+abs(mrlr2)**2
     & +abs(mrrl1)**2+abs(mrrl2)**2
     & -two/XN*real(mrrl1*conjg(mrrl2),dp)
     & +abs(mrll1)**2+abs(mrll2)**2
     & +abs(mlrr1)**2+abs(mlrr2)**2
     & +abs(mllr1)**2+abs(mllr2)**2
     & -two/XN*real(mllr1*conjg(mllr2),dp)
     & +abs(mlrl1)**2+abs(mlrl2)**2
     & +abs(mlll1)**2+abs(mlll2)**2
     & -two/XN*real(mlll1*conjg(mlll2),dp))

      return
      endif
      return
      end
