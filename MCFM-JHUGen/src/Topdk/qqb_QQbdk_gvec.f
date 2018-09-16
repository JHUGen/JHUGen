      subroutine qqb_QQbdk_gvec(p,n,in,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 2011.                                                  *
*----Matrix element for tt production                                  *
*----averaged over initial colours and spins                           *
*    line in contracted with the vector n(mu)                          *
*     g(-p1)+g(-p2)--> t(p3,p4,p5)+tb(p6,p7,p8)                        *
*                                                                      *
*     g(-p1) +g(-p2)= nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)    *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'plabel.f'
      include 'masses.f'
      include 'msqv_cs.f'
      include 'zprods_com.f'
C in is the label of the contracted line
      integer:: j,k,in,icol,nu
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),res(0:2),s,mt2
      real(dp):: n(4),fac,q(mxpart,4),s1t,s2t,s12
      real(dp):: c1,c2,c6,c8,c9,nDn,nDp4
      complex(dp):: prop
      complex(dp):: zanb(mxpart,mxpart)
      common/zprodsgvec/zanb
!$omp threadprivate(/zprodsgvec/)
      s(j,k)=2._dp
     & *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      call checkndotp(p,n,in)
C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      do icol=0,2
        msqv_cs(icol,j,k)=0._dp
      enddo
      enddo
      enddo


      s1t=s(1,3)+s(1,4)+s(1,5)
      s2t=s(2,3)+s(2,4)+s(2,5)
      s12=s(1,2)
      prop=cplx2(s(3,4)-wmass**2,wmass*wwidth)
     &    *cplx2(s(7,8)-wmass**2,wmass*wwidth)
     &    *cplx2(zip,mt*twidth)**2
      fac=V*xn/4._dp*gwsq**4*gsq**2/abs(prop)**2*s(5,3)*s(6,8)
      fac=fac*2._dp ! extra factor from ggn.frm

      mt2=mt**2
      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2
      nDp4=n(4)*p(4,4)-n(1)*p(4,1)-n(2)*p(4,2)-n(3)*p(4,3)
      c1=mt2/(s(3,4)+s(4,5))
      c2=mt2/(s(6,7)+s(7,8))
      c6=mt2/s1t+1._dp
      c8=mt2/s2t+1._dp
      c9=0.5_dp*nDn/nDp4

c--- pt=p3+p4+p5, ptb=p6+p7+p8
      do nu=1,4
      q(1,nu)=p(1,nu)
      q(2,nu)=p(2,nu)
      q(3,nu)=+(p(3,nu)+p(4,nu)+p(5,nu)-c1*p(4,nu))
      q(4,nu)=p(4,nu)
      q(5,nu)=p(6,nu)+(1._dp-c2)*p(7,nu)+p(8,nu)
      q(6,nu)=p(3,nu)+p(4,nu)+p(5,nu)+(1._dp-c6)*p(1,nu)
      q(7,nu)=p(7,nu)
      q(8,nu)=p(3,nu)+p(4,nu)+p(5,nu)+(1._dp-c8)*p(2,nu)
      q(9,nu)=n(nu)-c9*p(4,nu)
      enddo

      call spinoru(9,q,za,zb)
      do j=1,8
      do k=1,8
      zanb(j,k)=za(j,9)*zb(9,k)+c9*za(j,4)*zb(4,k)
      enddo
      enddo

      call ggttww1n(in,s1t,s2t,s12,c6,c8,res)

C--include factor for hadronic decays
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

      do icol=0,2
      msqv_cs(icol,0,0)=avegg*fac*res(icol)
      enddo
      msqv(0,0)=msqv_cs(0,0,0)+msqv_cs(1,0,0)+msqv_cs(2,0,0)

      return
      end
