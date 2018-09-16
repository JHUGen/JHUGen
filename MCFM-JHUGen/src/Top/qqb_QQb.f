      subroutine qqb_QQb(p,msq) 
      implicit none

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
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'msq_cs.f'
      include 'breit.f'
      include 'first.f'
      
      integer j,k,cs
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision wtqqb,wtgg,t1,t2,ro

      if (first) then
      first=.false.
      write(6,*) 'Heavy Quark mass:',mass2
      endif 

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      do cs=0,2
      msq_cs(cs,j,k)=0d0
      enddo
      enddo
      enddo
      call dotem(4,p,s)

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)
      ro=4d0*mass2**2/s(1,2)

      wtqqb=gsq**2*4d0/9d0*(t1**2+t2**2+ro/2d0)

c      wtgg=gsq**2*(1d0/6d0/t1/t2-3d0/8d0)
c     . *(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))
      msq_cs(1,0,0)=avegg*V*xn*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t2)-4d0*t1**2-0.5d0*(ro/t2)**2)
      msq_cs(2,0,0)=avegg*V*xn*gsq**2
     . *(-2d0*(1d0+ro)*(1d0-1d0/t1)-4d0*t2**2-0.5d0*(ro/t1)**2)
      msq_cs(0,0,0)=-avegg*V/xn*2d0*gsq**2
     . *(-2d0+(1d0+ro*(1d0-0.5d0*ro))/t1/t2
     .         -0.25d0*(ro/t1)**2-0.25d0*(ro/t2)**2)

      wtgg=msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)

C---fill qb-q, gg and q-qb elements
c--- the msq_cs entries for qqb and qbq are arbitrary
c--- divisions that are needed for summing up the _z contribution
c--- in virtint
      do j=-nf,nf
      k=-j
      if ((j .eq. 0) .and. (k.eq.0)) then
          msq(j,k)=wtgg
      elseif ((j .gt. 0) .and. (k.lt.0)) then
          msq(j,k)=wtqqb
          msq_cs(0,j,k)=wtqqb/3d0
          msq_cs(1,j,k)=wtqqb/3d0
          msq_cs(2,j,k)=wtqqb/3d0
      elseif ((j .lt. 0) .and. (k.gt.0)) then
          msq(j,k)=wtqqb
          msq_cs(0,j,k)=wtqqb/3d0
          msq_cs(1,j,k)=wtqqb/3d0
          msq_cs(2,j,k)=wtqqb/3d0
      endif
      enddo

      return
      end
