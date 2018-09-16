      subroutine qqb_gmgmjt_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
C***********************************************************************
*    Author: J.M. Campbell                                             *
*    March, 2013.                                                      *
c     Matrix element for gamma+gamma+jet production                    *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p5)              *
c     q(-p1)+qbar(-p2) --> gamma(p3) + gamma(p4) + g(p5)               *
C***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'nflav.f'
      include 'zprods_com.f'
      integer:: j,k,in
C--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),cfac
      real(dp):: gmgmjetn,fac,qqb,qbq,qg,gq,qbg,gqb,n(4)
      complex(dp):: zanb(mxpart,mxpart),zbna(mxpart,mxpart)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(5,p,s)

      fac=4._dp*cf*xn*esq**2*gsq
      call spinoru(5,p,za,zb)
      call spinork(5,p,zanb,zbna,n)

      if (in == 1) then
        call checkndotp(p,n,1)
        gqb=aveqg*fac*gmgmjetn(zanb,5,2,1,3,4)
        gq=aveqg*fac*gmgmjetn(zanb,2,5,1,3,4)
      elseif (in == 2) then
        call checkndotp(p,n,2)
        qg=aveqg*fac*gmgmjetn(zanb,1,5,2,3,4)
        qbg=aveqg*fac*gmgmjetn(zanb,5,1,2,3,4)
      elseif (in == 5) then
        call checkndotp(p,n,5)
        qbq=+aveqq*fac*gmgmjetn(zanb,2,1,5,3,4)
        qqb=+aveqq*fac*gmgmjetn(zanb,1,2,5,3,4)
      endif

      do j=1,nf
        cfac=Q(j)**4
        msq(j,-j)=cfac*qqb
        msq(-j,j)=cfac*qbq
        msq(j,0)=cfac*qg
        msq(0,j)=cfac*gq
        msq(-j,0)=cfac*qbg
        msq(0,-j)=cfac*gqb
      enddo


      return
      end

