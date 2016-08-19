      subroutine qqb_w(p,msq)
      implicit none
c----Matrix element for W production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c---
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qqb,qbq,s
     
c--statement function 
      s(j,k)=2d0*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     .           -p(j,2)*p(k,2)-p(j,3)*p(k,3))


      fac=gw**4*xn
c--   calculate propagator
      fac=aveqq*fac/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
c---case dbar-u or ubar-d
      qqb=fac*s(1,4)**2
      qbq=fac*s(2,4)**2

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
          if ((j .gt. 0) .and. (k .lt. 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo
      return
      end
