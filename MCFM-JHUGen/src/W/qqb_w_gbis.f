      subroutine qqb_w_gbis(p,msq)
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
      include 'mxpart.f'
      include 'nf.f'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'hpls.f'
      integer:: j,k
      integer,parameter:: region2=2,region3=3,region4=4
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,propsq
      real(dp):: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq,Wampqqbgsq

      msq(:,:)=zip

      call spinoru(5,p,za,zb)
c---calculate the propagator
      propsq=(s(3,4)-wmass**2)**2+(wmass*wwidth)**2
      fac=gwsq**2*gsq*V/propsq


c      qqbWg= +aveqq*fac*Wampqqbgsq(1,2,5,3,4,region2,za,zb)
c      gqbWqb=+aveqg*fac*Wampqqbgsq(5,2,1,3,4,region4,za,zb)
c      qgWq = +aveqg*fac*Wampqqbgsq(1,5,2,3,4,region3,za,zb)

c      qbqWg= +aveqq*fac*Wampqqbgsq(2,1,5,3,4,region2,za,zb)
c      qbgWqb=+aveqg*fac*Wampqqbgsq(5,1,2,3,4,region3,za,zb)
c      gqWq  =+aveqg*fac*Wampqqbgsq(2,5,1,3,4,region4,za,zb)

      qqbWg= aveqq*fac*Wampqqbgsq(2,1,5,3,4,region2,za,zb)
      gqbWqb=aveqg*fac*Wampqqbgsq(2,5,1,3,4,region3,za,zb)
      qgWq = aveqg*fac*Wampqqbgsq(5,1,2,3,4,region4,za,zb)

      qbqWg= aveqq*fac*Wampqqbgsq(1,2,5,3,4,region2,za,zb)
      qbgWqb=aveqg*fac*Wampqqbgsq(1,5,2,3,4,region3,za,zb)
      gqWq  =aveqg*fac*Wampqqbgsq(5,2,1,3,4,region4,za,zb)

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

