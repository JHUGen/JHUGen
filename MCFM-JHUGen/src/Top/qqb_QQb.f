      subroutine qqb_QQb(p,msq)
      implicit none
      include 'types.f'


************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C      This is the four dimensional result for                         *
C      Quark antiquark annihilation in order alfa_s^2                  *
C      q(P1) + qbar(P2) --> Q(-P3) + Qbar(-P4)                         *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'msq_cs.f'
      include 'breit.f'
      include 'first.f'

      integer:: j,k,cs
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: wtqqb,wtgg,t1,t2,ro

      if (first) then
      first=.false.
      write(6,*) 'Heavy Quark mass:',mass2
      endif

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      do cs=0,2
      msq_cs(cs,j,k)=0._dp
      enddo
      enddo
      enddo
      call dotem(4,p,s)

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)
      ro=4._dp*mass2**2/s(1,2)

      wtqqb=gsq**2*4._dp/9._dp*(t1**2+t2**2+ro/2._dp)

c      wtgg=gsq**2*(1._dp/6._dp/t1/t2-3._dp/8._dp)
c     & *(t1**2+t2**2+ro-0.25_dp*ro**2/(t1*t2))
      msq_cs(1,0,0)=avegg*V*xn*gsq**2
     & *(-2._dp*(1._dp+ro)*(1._dp-1._dp/t2)-4._dp*t1**2-0.5_dp*(ro/t2)**2)
      msq_cs(2,0,0)=avegg*V*xn*gsq**2
     & *(-2._dp*(1._dp+ro)*(1._dp-1._dp/t1)-4._dp*t2**2-0.5_dp*(ro/t1)**2)
      msq_cs(0,0,0)=-avegg*V/xn*2._dp*gsq**2
     & *(-2._dp+(1._dp+ro*(1._dp-0.5_dp*ro))/t1/t2
     &         -0.25_dp*(ro/t1)**2-0.25_dp*(ro/t2)**2)

      wtgg=msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)

C---fill qb-q, gg and q-qb elements
c--- the msq_cs entries for qqb and qbq are arbitrary
c--- divisions that are needed for summing up the _z contribution
c--- in virtint
      do j=-nf,nf
      k=-j
      if ((j == 0) .and. (k==0)) then
          msq(j,k)=wtgg
      elseif ((j > 0) .and. (k<0)) then
          msq(j,k)=wtqqb
          msq_cs(0,j,k)=wtqqb/3._dp
          msq_cs(1,j,k)=wtqqb/3._dp
          msq_cs(2,j,k)=wtqqb/3._dp
      elseif ((j < 0) .and. (k>0)) then
          msq(j,k)=wtqqb
          msq_cs(0,j,k)=wtqqb/3._dp
          msq_cs(1,j,k)=wtqqb/3._dp
          msq_cs(2,j,k)=wtqqb/3._dp
      endif
      enddo

      return
      end
