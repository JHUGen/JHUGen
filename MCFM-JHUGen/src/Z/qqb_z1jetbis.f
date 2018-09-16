      subroutine qqb_z1jetbis(p,msq)
      implicit none
      include 'types.f'

C-----Authors: John Campbell, Keith Ellis
C-----June 2000 and December 2001
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j,k,hq,hl
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      complex(dp):: prop,al,be,ga,de
      real(dp):: AqqbZg2(2,2),AqbqZg2(2,2),AqgZq2(2,2),
     &               AqbgZqb2(2,2),AgqbZqb2(2,2),AgqZq2(2,2)
      integer,parameter:: swap(2)=(/2,1/)
      integer,parameter:: region2=2,region3=3,region4=4

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(5,p,za,zb)

      al=cone
      be=cone
      ga=czip
      de=half*s(1,2)*(al-be-ga)/(s(1,2)+s(2,3)+s(3,1))

c---protect from soft and collinear singularities
c      if  ((-s(1,5) < cutoff) .or. (-s(2,5) < cutoff)) return

C-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) < 4._dp*mbsq) return

      prop=s(3,4)/cmplx((s(3,4)-zmass**2),zmass*zwidth,dp)
      fac=4._dp*V*esq**2*gsq


      call zgampsbis(0,2,1,5,3,4,region2,za,zb,AqqbZg2)
      call zgampsbis(0,2,5,1,3,4,region3,za,zb,AgqbZqb2)
      call zgampsbis(0,5,1,2,3,4,region4,za,zb,AqgZq2)

      call zgampsbis(0,1,2,5,3,4,region2,za,zb,AqbqZg2)
      call zgampsbis(0,1,5,2,3,4,region3,za,zb,AqbgZqb2)
      call zgampsbis(0,5,2,1,3,4,region4,za,zb,AgqZq2)

      do j=-nflav,nflav
      do k=-nflav,nflav

      if( (j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 20

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*AqqbZg2(1,1)
     &            +abs(Q(j)*q1+L(j)*r1*prop)**2*AqqbZg2(1,2)
     &            +abs(Q(j)*q1+R(j)*l1*prop)**2*AqqbZg2(2,1)
     &            +abs(Q(j)*q1+R(j)*r1*prop)**2*AqqbZg2(2,2)
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*AqbqZg2(1,1)
     &            +abs(Q(k)*q1+L(k)*r1*prop)**2*AqbqZg2(1,2)
     &            +abs(Q(k)*q1+R(k)*l1*prop)**2*AqbqZg2(2,1)
     &            +abs(Q(k)*q1+R(k)*r1*prop)**2*AqbqZg2(2,2)
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*AqgZq2(1,1)
     &            +abs(Q(j)*q1+L(j)*r1*prop)**2*AqgZq2(1,2)
     &            +abs(Q(j)*q1+R(j)*l1*prop)**2*AqgZq2(2,1)
     &            +abs(Q(j)*q1+R(j)*r1*prop)**2*AqgZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(-j)*q1+L(-j)*l1*prop)**2*AqbgZqb2(1,1)
     &            +abs(Q(-j)*q1+L(-j)*r1*prop)**2*AqbgZqb2(1,2)
     &            +abs(Q(-j)*q1+R(-j)*l1*prop)**2*AqbgZqb2(2,1)
     &            +abs(Q(-j)*q1+R(-j)*r1*prop)**2*AqbgZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*AgqZq2(1,1)
     &            +abs(Q(k)*q1+L(k)*r1*prop)**2*AgqZq2(1,2)
     &            +abs(Q(k)*q1+R(k)*l1*prop)**2*AgqZq2(2,1)
     &            +abs(Q(k)*q1+R(k)*r1*prop)**2*AgqZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(-k)*q1+L(-k)*l1*prop)**2*AgqbZqb2(1,1)
     &            +abs(Q(-k)*q1+L(-k)*r1*prop)**2*AgqbZqb2(1,2)
     &            +abs(Q(-k)*q1+R(-k)*l1*prop)**2*AgqbZqb2(2,1)
     &            +abs(Q(-k)*q1+R(-k)*r1*prop)**2*AgqbZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      endif

 20   continue
      enddo
      enddo
      return
      end


