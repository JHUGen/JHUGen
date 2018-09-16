      subroutine qqb_QQb_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process with gluon line in contracted with the vector n  *
*     This is the four dimensional result for                          *
*     Heavy quark production in order alfa_s^2                         *
*     f(P1) + f(P2) --> Q(-P3) + Qbar(-P4)                             *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'msqv_cs.f'
      include 'breit.f'
      include 'first.f'
      integer:: j,k,in
C--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4),msqn(0:2)

      if (first) then
      first=.false.
      write(6,*) 'qqb_QQb_gvec:mass2',mass2
      endif


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(4,p,s)

      call checkndotp(p,n,in)

      if     (in == 1) then
        call qqb_QQbn(1,2,mass2,p,n,msqn)
      elseif (in == 2) then
        call qqb_QQbn(2,1,mass2,p,n,msqn)
      endif

      do j=0,2
      msqv_cs(j,0,0)=avegg*gsq**2*msqn(j)
      enddo
      msq(0,0)=msqv_cs(0,0,0)+msqv_cs(1,0,0)+msqv_cs(2,0,0)

      return
      end




      subroutine qqb_QQbn(i1,i2,mass,p,n,msqn)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      integer:: i1,i2
C--in is the label of the parton dotted with n
      real(dp):: n(4),nDn,nDt,nDtb,nDp2,t1,t2,ro,mass,p(mxpart,4),
     & msqn(0:2)

      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2
      nDt=n(4)*p(3,4)-n(1)*p(3,1)-n(2)*p(3,2)-n(3)*p(3,3)
      nDp2=n(4)*p(i2,4)-n(1)*p(i2,1)-n(2)*p(i2,2)-n(3)*p(i2,3)
      nDtb=-nDp2-nDt
      t1=-s(i1,3)/s(i1,i2)
      t2=-s(i2,3)/s(i1,i2)
      ro=4._dp*mass**2/s(i1,i2)

      msqn(0)=V/xn*(nDn*(1._dp/(t1*t2)-2._dp)
     & -2._dp*ro/s(1,2)*((nDt+t1*nDp2)/(t1*t2))**2)
      msqn(i1)=V*xn*(nDn*(2._dp*t1-1._dp/t2+t1**2+t2**2)
     & +2._dp*ro/s(1,2)*(nDtb/t2+nDp2)**2)
      msqn(i2)=V*xn*(nDn*(2._dp*t2-1._dp/t1+t1**2+t2**2)
     & +2._dp*ro/s(1,2)*(nDt/t1+nDp2)**2)

c      qqb_QQbx=2._dp*V*(V/(xn*t1*t2)-2._dp*xn)
c     & *(nDn*(t1*t2-0.5_dp)+ro*(nDt+t1*nDp2)**2/(t1*t2*s(i1,i2)))
      return
      end
