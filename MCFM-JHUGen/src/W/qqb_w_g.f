      subroutine qqb_w_g(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     u(-p1)+dbar(-p2)--> W^+(n(p3)+e^+(p4))   + g(p5)
C For nwz=-1
c     d(-p1)+ubar(-p2)--> W^-(e^-(p3)+nbar(p4))+ g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq,w1jet


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zip
      enddo
      enddo

      call dotem(5,p,s)
c---calculate the propagator
      fac=gwsq**2*gsq*V

      qqbWg= +aveqq*fac*w1jet(1,2,3,4,5)
      gqbWqb=-aveqg*fac*w1jet(5,2,3,4,1)
      qgWq=  -aveqg*fac*w1jet(1,5,3,4,2)

      qbqWg= +aveqq*fac*w1jet(2,1,3,4,5)
      qbgWqb=-aveqg*fac*w1jet(5,1,3,4,2)
      gqWq=  -aveqg*fac*w1jet(2,5,3,4,1)

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
      endif

      enddo
      enddo
      return
      end

      function w1jet(j1,j2,j3,j4,j5)
      implicit none
      include 'types.f'
      real(dp):: w1jet
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5
      real(dp):: prop

      prop=((s(j3,j4)-wmass**2)**2+(wmass*wwidth)**2)

      w1jet=(s(j1,j4)**2+s(j2,j3)**2)*s(j3,j4)/(s(j1,j5)*s(j2,j5)*prop)

      return
      end

