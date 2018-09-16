      subroutine qqb_QQbdku_gvec(p,n,in,msqv)
      implicit none
      include 'types.f'


c----Matrix element for tt production
C----averaged over initial colours and spins
c    line in contracted with the vector n(mu)
c     g(-p1)+g(-p2)--> t(p3,p4,p5)+tb(p6,p7,p8)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'msqv_cs.f'

C in is the label of the contracted line
      integer:: j,k,in,icol,nu
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),q(mxpart,4),s
      real(dp):: n(4),nDp1
      real(dp):: pt1Dpe1,pb1Dpn1,pt2Dpe2,pb2Dpn2,decay
      complex(dp)::  prop
      s(j,k)=2._dp
     & *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      prop=cplx2(s(3,4)-wmass**2,wmass*wwidth)
     &    *cplx2(s(7,8)-wmass**2,wmass*wwidth)
     &    *cplx2(zip,mt*twidth)**2


C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      do icol=0,2
        msqv_cs(icol,j,k)=0._dp
      enddo
      enddo
      enddo

      nDp1=n(4)*p(in,4)-n(1)*p(in,1)-n(2)*p(in,2)-n(3)*p(in,3)
      if (abs(nDp1)>1.e-3_dp*abs(p(1,4))) then
         write(*,*) 'Error for in=:',in
         write(*,*) 'cutoff',1.e-3_dp*abs(p(1,4))
         write(6,*) 'qqb_QQbdk_gvec:nDp1',nDp1
         call flush(6)
         stop
      endif

      do nu=1,4
      q(1,nu)=p(1,nu)
      q(2,nu)=p(2,nu)
      q(3,nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      q(4,nu)=+p(6,nu)+p(7,nu)+p(8,nu)
      do j=5,mxpart
      q(j,nu)=0._dp
      enddo
      enddo

      pt1Dpe1=q(3,4)*p(4,4)-q(3,1)*p(4,1)-q(3,2)*p(4,2)-q(3,3)*p(4,3)
      pt2Dpe2=q(4,4)*p(7,4)-q(4,1)*p(7,1)-q(4,2)*p(7,2)-q(4,3)*p(7,3)
      pb1Dpn1=p(3,4)*p(5,4)-p(3,1)*p(5,1)-p(3,2)*p(5,2)-p(3,3)*p(5,3)
      pb2Dpn2=p(6,4)*p(8,4)-p(6,1)*p(8,1)-p(6,2)*p(8,2)-p(6,3)*p(8,3)
      decay=4._dp*gwsq**4*pt1Dpe1*pb1Dpn1*pt2Dpe2*pb2Dpn2/abs(prop)**2

      call qqb_QQb_gvec(q,n,in,msqv)

      do icol=0,2
      msqv_cs(icol,0,0)=msqv_cs(icol,0,0)*decay
      enddo

      msqv(0,0)=msqv_cs(0,0,0)+msqv_cs(1,0,0)+msqv_cs(2,0,0)


      return
      end

