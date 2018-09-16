      subroutine qqb_zh_gaga(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1998.                                                  *
*  Matrix element squared averaged over initial colors and spins       *
*     q(-p1)+qbar(-p2) -->  H  + Z                                     *
*                           |    |                                     *
*                           |     ->fermion(p3)+antifermion(p4)        *
*                           |                                          *
*                            ---> gamma(p5)+gamma(p6)                  *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer:: j,k
      real(dp):: p(mxpart,4)
      real(dp):: s,prop,fac,q1423,q2413,s56,v2(2)
      real(dp):: msq(-nf:nf,-nf:nf),hdecay
      real(dp):: msqhgamgam

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      v2(1)=l1
      v2(2)=r1

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo
      s56=s(5,6)

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop

      hdecay=msqhgamgam(s56)/((s56-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

c-- Old form of this matrix element (modified to facilitate extension
c--- to H->WW decay)
c      spinave only (color factors cancel)
c       fac=two*spinave*gw**8*(xw/(1._dp-xw))**2*mbsq*(s56-4._dp*mb**2)/prop


      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)

      do j=-nf,nf
      if (j == 0) go to 40
      k=-j
      if ((j > 0) .and. (k < 0)) then
      msq(j,k)=
     &  +((l(j)*v2(1))**2+(r(j)*v2(2))**2)*q1423
     &  +((l(j)*v2(2))**2+(r(j)*v2(1))**2)*q2413
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=
     &  +((l(k)*v2(1))**2+(r(k)*v2(2))**2)*q2413
     &  +((l(k)*v2(2))**2+(r(k)*v2(1))**2)*q1423
      endif
 40   continue
      enddo
      return
      end

c + L1(j)*L1(j)*L2(j)*L2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s14*s23*cos^-4 )
 
c + R1(j)*R1(j)*R2(j)*R2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s14*s23*cos^-4 )

c + L1(j)*L1(j)*R2(j)*R2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s13*s24*cos^-4 )
 
c + R1(j)*R1(j)*L2(j)*L2(j)*gw^8*[1/4/XN]*mb^2*[s56-4*mb^2]
c  * ( 2*zprop1^-1*zprop2^-1*hprop^-1*sinw^4*s13*s24*cos^-4 )
 

