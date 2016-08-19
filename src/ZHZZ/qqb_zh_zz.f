      subroutine qqb_zh_zz(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1998.                                                  *
*  Matrix element squared averaged over initial colors and spins       *
*     q(-p1)+qbar(-p2) -->  H  + Z                                     *
*                           |    |                                     *
*                           |     ->fermion(p3)+antifermion(p4)        *
*                           |                                          *
*                           -> Z(e^-(p5),e^+(p6)) Z(mu^-(p7),mu^+(p8)) *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer j,k
      double precision p(mxpart,4)
      double precision s,prop,fac,q1423,q2413,s5678,v2(2)
      double precision msq(-nf:nf,-nf:nf),hdecay

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      v2(1)=l1
      v2(2)=r1

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      s5678=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      fac=xn*4d0*(xw/(1d0-xw))**2*gwsq**3*wmass**2/prop
      hdecay=gwsq**3*zmass**2*4d0*xw**2/(one-xw)*
     . ( ((l1*l2)**2+(r1*r2)**2)*s(5,7)*s(6,8)
     .  +((r1*l2)**2+(r2*l1)**2)*s(5,8)*s(6,7))
      hdecay=hdecay/(((s5678-hmass**2)**2+(hmass*hwidth)**2)
     .   *((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
     .   *((s(7,8)-zmass**2)**2+(zmass*zwidth)**2))
      fac=fac*hdecay

      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)

      do j=-nf,nf
      if (j .eq. 0) go to 40
      k=-j
      if ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=
     .  +((l(j)*v2(1))**2+(r(j)*v2(2))**2)*q1423
     .  +((l(j)*v2(2))**2+(r(j)*v2(1))**2)*q2413
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=
     .  +((l(k)*v2(1))**2+(r(k)*v2(2))**2)*q2413
     .  +((l(k)*v2(2))**2+(r(k)*v2(1))**2)*q1423
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
 
