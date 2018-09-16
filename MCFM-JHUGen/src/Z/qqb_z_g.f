      subroutine qqb_z_g(p,msq)
      implicit none
      include 'types.f'

C-----Author John Campbell
C-----June 2000
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      complex(dp):: AqqbZg(2,2,2),AqbqZg(2,2,2),AqgZq(2,2,2),
     &               AqbgZqb(2,2,2),AgqbZqb(2,2,2),AgqZq(2,2,2),prop

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(5,p,s)
      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities
c      if  ((-s(1,5) < cutoff) .or. (-s(2,5) < cutoff)) return

C-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) < 4._dp*mbsq) return

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      fac=4._dp*V*esq**2*gsq

c      qqbZg= +aveqq*s(3,4)**2*fac*z1jet(1,2,3,4,5)
c      gqbZqb=-aveqg*s(3,4)**2*fac*z1jet(5,2,3,4,1)
c      qgZq=  -aveqg*s(3,4)**2*fac*z1jet(1,5,3,4,2)
c      qbqZg= +aveqq*s(3,4)**2*fac*z1jet(2,1,3,4,5)
c      qbgZqb=-aveqg*s(3,4)**2*fac*z1jet(5,1,3,4,2)
c      gqZq=  -aveqg*s(3,4)**2*fac*z1jet(2,5,3,4,1)

      call zgamps(1,2,3,4,5,za,zb,AqqbZg)
      call zgamps(5,2,3,4,1,za,zb,AgqbZqb)
      call zgamps(1,5,3,4,2,za,zb,AqgZq)
      call zgamps(2,1,3,4,5,za,zb,AqbqZg)
      call zgamps(5,1,3,4,2,za,zb,AqbgZqb)
      call zgamps(2,5,3,4,1,za,zb,AgqZq)

      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=abs((Q(j)*q1+L(j)*l1*prop)*AqqbZg(1,1,1))**2
     &            +abs((Q(j)*q1+L(j)*l1*prop)*AqqbZg(1,1,2))**2
     &            +abs((Q(j)*q1+L(j)*r1*prop)*AqqbZg(1,2,1))**2
     &            +abs((Q(j)*q1+L(j)*r1*prop)*AqqbZg(1,2,2))**2
     &            +abs((Q(j)*q1+R(j)*l1*prop)*AqqbZg(2,1,1))**2
     &            +abs((Q(j)*q1+R(j)*l1*prop)*AqqbZg(2,1,2))**2
     &            +abs((Q(j)*q1+R(j)*r1*prop)*AqqbZg(2,2,1))**2
     &            +abs((Q(j)*q1+R(j)*r1*prop)*AqqbZg(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=abs((Q(k)*q1+L(k)*l1*prop)*AqbqZg(1,1,1))**2
     &            +abs((Q(k)*q1+L(k)*l1*prop)*AqbqZg(1,1,2))**2
     &            +abs((Q(k)*q1+L(k)*r1*prop)*AqbqZg(1,2,1))**2
     &            +abs((Q(k)*q1+L(k)*r1*prop)*AqbqZg(1,2,2))**2
     &            +abs((Q(k)*q1+R(k)*l1*prop)*AqbqZg(2,1,1))**2
     &            +abs((Q(k)*q1+R(k)*l1*prop)*AqbqZg(2,1,2))**2
     &            +abs((Q(k)*q1+R(k)*r1*prop)*AqbqZg(2,2,1))**2
     &            +abs((Q(k)*q1+R(k)*r1*prop)*AqbqZg(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=abs((Q(j)*q1+L(j)*l1*prop)*AqgZq(1,1,1))**2
     &            +abs((Q(j)*q1+L(j)*l1*prop)*AqgZq(1,1,2))**2
     &            +abs((Q(j)*q1+L(j)*r1*prop)*AqgZq(1,2,1))**2
     &            +abs((Q(j)*q1+L(j)*r1*prop)*AqgZq(1,2,2))**2
     &            +abs((Q(j)*q1+R(j)*l1*prop)*AqgZq(2,1,1))**2
     &            +abs((Q(j)*q1+R(j)*l1*prop)*AqgZq(2,1,2))**2
     &            +abs((Q(j)*q1+R(j)*r1*prop)*AqgZq(2,2,1))**2
     &            +abs((Q(j)*q1+R(j)*r1*prop)*AqgZq(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=abs((Q(-j)*q1+L(-j)*l1*prop)*AqbgZqb(1,1,1))**2
     &            +abs((Q(-j)*q1+L(-j)*l1*prop)*AqbgZqb(1,1,2))**2
     &            +abs((Q(-j)*q1+L(-j)*r1*prop)*AqbgZqb(1,2,1))**2
     &            +abs((Q(-j)*q1+L(-j)*r1*prop)*AqbgZqb(1,2,2))**2
     &            +abs((Q(-j)*q1+R(-j)*l1*prop)*AqbgZqb(2,1,1))**2
     &            +abs((Q(-j)*q1+R(-j)*l1*prop)*AqbgZqb(2,1,2))**2
     &            +abs((Q(-j)*q1+R(-j)*r1*prop)*AqbgZqb(2,2,1))**2
     &            +abs((Q(-j)*q1+R(-j)*r1*prop)*AqbgZqb(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=abs((Q(k)*q1+L(k)*l1*prop)*AgqZq(1,1,1))**2
     &            +abs((Q(k)*q1+L(k)*l1*prop)*AgqZq(1,1,2))**2
     &            +abs((Q(k)*q1+L(k)*r1*prop)*AgqZq(1,2,1))**2
     &            +abs((Q(k)*q1+L(k)*r1*prop)*AgqZq(1,2,2))**2
     &            +abs((Q(k)*q1+R(k)*l1*prop)*AgqZq(2,1,1))**2
     &            +abs((Q(k)*q1+R(k)*l1*prop)*AgqZq(2,1,2))**2
     &            +abs((Q(k)*q1+R(k)*r1*prop)*AgqZq(2,2,1))**2
     &            +abs((Q(k)*q1+R(k)*r1*prop)*AgqZq(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=abs((Q(-k)*q1+L(-k)*l1*prop)*AgqbZqb(1,1,1))**2
     &            +abs((Q(-k)*q1+L(-k)*l1*prop)*AgqbZqb(1,1,2))**2
     &            +abs((Q(-k)*q1+L(-k)*r1*prop)*AgqbZqb(1,2,1))**2
     &            +abs((Q(-k)*q1+L(-k)*r1*prop)*AgqbZqb(1,2,2))**2
     &            +abs((Q(-k)*q1+R(-k)*l1*prop)*AgqbZqb(2,1,1))**2
     &            +abs((Q(-k)*q1+R(-k)*l1*prop)*AgqbZqb(2,1,2))**2
     &            +abs((Q(-k)*q1+R(-k)*r1*prop)*AgqbZqb(2,2,1))**2
     &            +abs((Q(-k)*q1+R(-k)*r1*prop)*AgqbZqb(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      endif

   19 continue
      enddo
      enddo
      return
      end


c      function z1jet(j1,j2,j3,j4,j5)
c      implicit none
      include 'types.f'
      real(dp):: z1jet
c
c      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'sprods_com.f'
c      integer:: j1,j2,j3,j4,j5
c      real(dp):: s12,s15,s25

c      s12=s(j1,j2)
c      s15=s(j1,j5)
c      s25=s(j2,j5)
c---calculate the propagator
c      z1jet=(s(j1,j4)**2+s(j2,j3)**2)/(s25*s15*s(j3,j4))
cc      z1jet=
cc     & (s12*(2._dp*s(j1,j4)*s(j2,j3)+s(j1,j4)*s(j3,j5)+s(j2,j3)*s(j4,j5))
cc     & +s15*(s(j1,j4)*s(j2,j3)+s(j1,j4)*s(j3,j5)-s(j2,j3)*s(j2,j4))
cc     & +s25*(s(j1,j4)*s(j2,j3)+s(j2,j3)*s(j4,j5)-s(j1,j3)*s(j1,j4)))
cc     & /(s15*s25*s(j3,j4)**2)
c      return
c      end


c      subroutine zgamps(j1,j2,j3,j4,j5,za,zb,amps)
      implicit none
      include 'types.f'
c
c      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'zprods_decl.f'
c      complex(dp):: amps(2,2,2)
c      integer:: h1,h2,j1,j2,j3,j4,j5
c-- amplitude helicities are amps(quark,lepton,gluon)

c      amps(1,1,1)=za(j2,j3)/za(j1,j5)/za(j2,j5)
c     &             *(za(j2,j1)*zb(j4,j1)+za(j2,j5)*zb(j4,j5))
c
c      amps(1,1,2)=zb(j4,j1)/zb(j1,j5)/zb(j2,j5)
c     &             *(za(j2,j3)*zb(j2,j1)+za(j3,j5)*zb(j1,j5))
c
c      amps(1,2,1)=za(j2,j4)/za(j1,j5)/za(j2,j5)
c     &             *(za(j2,j1)*zb(j3,j1)+za(j2,j5)*zb(j3,j5))
c
c      amps(1,2,2)=zb(j3,j1)/zb(j1,j5)/zb(j2,j5)
c     &             *(za(j2,j4)*zb(j2,j1)+za(j4,j5)*zb(j1,j5))
c
c      do h1=1,2
c      do h2=1,2
c        amps(2,h1,h2)=-conjg(amps(1,3-h1,3-h2))
c      enddo
c      enddo
c
c
c      return
c      end

