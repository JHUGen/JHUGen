      subroutine qqb_z(p,msq)
      implicit none
      include 'types.f'

C-----Author John Campbell
C-----June 2000
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s,fac,s34
      complex(dp):: prop,qqb,qbq

c---statement function
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo
      s34=s(3,4)
c      if (s34 < 4._dp*mbsq) return

      fac=4._dp*esq**2*xn

c--   calculate propagators
      fac=aveqq*fac/s34**2
      prop=s34/cplx2((s34-zmass**2),zmass*zwidth)

      call spinoru(4,p,za,zb)

c---case qbar-q or q-qbar
c      qqb=fac*s(1,4)**2
c      qbq=fac*s(2,4)**2
      qqb=za(2,3)*zb(4,1)
      qbq=za(1,3)*zb(4,2)

      do j=-nf,nf
      k=-j
          if ((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
          elseif ((j > 0) .and. (k < 0)) then
            msq(j,k)=+abs((Q(j)*q1+L(j)*l1*prop)*qqb)**2
     &               +abs((Q(j)*q1+R(j)*r1*prop)*qqb)**2
     &               +abs((Q(j)*q1+L(j)*r1*prop)*qbq)**2
     &               +abs((Q(j)*q1+R(j)*l1*prop)*qbq)**2
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=+abs((Q(k)*q1+L(k)*l1*prop)*qbq)**2
     &               +abs((Q(k)*q1+R(k)*r1*prop)*qbq)**2
     &               +abs((Q(k)*q1+L(k)*r1*prop)*qqb)**2
     &               +abs((Q(k)*q1+R(k)*l1*prop)*qqb)**2
          endif
          msq(j,k)=msq(j,k)*fac
      enddo

      return
      end
