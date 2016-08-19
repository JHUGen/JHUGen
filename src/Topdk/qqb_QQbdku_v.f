      subroutine qqb_QQbdku_v(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     June 2010.                                                       *
*     Virtual matrix elements squared for the process                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)  *
*     Uncorrelated decays                                              * 
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),q(mxpart,4),p(mxpart,4),s
      double precision pt1Dpe1,pb1Dpn1,pt2Dpe2,pb2Dpn2,decay
      double complex  prop
      s(j,k)=2d0
     . *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      prop=dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     .    *dcmplx(s(7,8)-wmass**2,wmass*wwidth)
     .    *dcmplx(zip,mt*twidth)**2

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      q(1,nu)=p(1,nu)
      q(2,nu)=p(2,nu)
      q(3,nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      q(4,nu)=+p(6,nu)+p(7,nu)+p(8,nu)
      do j=5,mxpart
      q(j,nu)=0d0
      enddo
      enddo

      pt1Dpe1=q(3,4)*p(4,4)-q(3,1)*p(4,1)-q(3,2)*p(4,2)-q(3,3)*p(4,3)
      pt2Dpe2=q(4,4)*p(7,4)-q(4,1)*p(7,1)-q(4,2)*p(7,2)-q(4,3)*p(7,3)
      pb1Dpn1=p(3,4)*p(5,4)-p(3,1)*p(5,1)-p(3,2)*p(5,2)-p(3,3)*p(5,3)
      pb2Dpn2=p(6,4)*p(8,4)-p(6,1)*p(8,1)-p(6,2)*p(8,2)-p(6,3)*p(8,3)
      decay=4d0*gwsq**4*pt1Dpe1*pb1Dpn1*pt2Dpe2*pb2Dpn2/abs(prop)**2

      call qqb_QQb_v(q,msq)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=msq(j,k)*decay
      enddo
      enddo
      return
      end
