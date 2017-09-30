      subroutine qqb_z1jet(p,msq)
      implicit none
C-----Authors: John Campbell, Keith Ellis
C-----June 2000 and December 2001
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j,k,hq,hl
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double complex prop
      double precision AqqbZg2(2,2),AqbqZg2(2,2),AqgZq2(2,2),
     .               AqbgZqb2(2,2),AgqbZqb2(2,2),AgqZq2(2,2)
      integer,parameter:: swap(2)=(/2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities
c      if  ((-s(1,5) .lt. cutoff) .or. (-s(2,5) .lt. cutoff)) return

C-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) .lt. 4d0*mbsq) return

      prop=s(3,4)/Dcmplx((s(3,4)-zmass**2),zmass*zwidth)
      fac=4d0*V*esq**2*gsq

c      qqbZg= +aveqq*s(3,4)**2*fac*z1jet(1,2,3,4,5)
c      gqbZqb=-aveqg*s(3,4)**2*fac*z1jet(5,2,3,4,1)
c      qgZq=  -aveqg*s(3,4)**2*fac*z1jet(1,5,3,4,2)
c      qbqZg= +aveqq*s(3,4)**2*fac*z1jet(2,1,3,4,5)
c      qbgZqb=-aveqg*s(3,4)**2*fac*z1jet(5,1,3,4,2)
c      gqZq=  -aveqg*s(3,4)**2*fac*z1jet(2,5,3,4,1)

      call zgamps2(1,2,3,4,5,za,zb,AqqbZg2)
      call zgamps2(1,5,3,4,2,za,zb,AqgZq2)
      call zgamps2(2,5,3,4,1,za,zb,AgqZq2)
      do hq=1,2
      do hl=1,2
      AqbqZg2(hq,hl)=AqqbZg2(hq,swap(hl))
      AqbgZqb2(hq,hl)=AqgZq2(hq,swap(hl))
      AgqbZqb2(hq,hl)=AgqZq2(hq,swap(hl))
      enddo
      enddo

c      call zgamps2(2,1,3,4,5,za,zb,AqbqZg2)
c      call zgamps2(5,1,3,4,2,za,zb,AqbgZqb2)
c      call zgamps2(5,2,3,4,1,za,zb,AgqbZqb2)

      do j=-nflav,nflav
      do k=-nflav,nflav

      if( (j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 20

      if     ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=0d0
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=cdabs(Q(j)*q1+L(j)*l1*prop)**2*AqqbZg2(1,1)
     .            +cdabs(Q(j)*q1+L(j)*r1*prop)**2*AqqbZg2(1,2)
     .            +cdabs(Q(j)*q1+R(j)*l1*prop)**2*AqqbZg2(2,1)
     .            +cdabs(Q(j)*q1+R(j)*r1*prop)**2*AqqbZg2(2,2)
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=cdabs(Q(k)*q1+L(k)*l1*prop)**2*AqbqZg2(1,1)
     .            +cdabs(Q(k)*q1+L(k)*r1*prop)**2*AqbqZg2(1,2)
     .            +cdabs(Q(k)*q1+R(k)*l1*prop)**2*AqbqZg2(2,1)
     .            +cdabs(Q(k)*q1+R(k)*r1*prop)**2*AqbqZg2(2,2)
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=cdabs(Q(j)*q1+L(j)*l1*prop)**2*AqgZq2(1,1)
     .            +cdabs(Q(j)*q1+L(j)*r1*prop)**2*AqgZq2(1,2)
     .            +cdabs(Q(j)*q1+R(j)*l1*prop)**2*AqgZq2(2,1)
     .            +cdabs(Q(j)*q1+R(j)*r1*prop)**2*AqgZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=cdabs(Q(-j)*q1+L(-j)*l1*prop)**2*AqbgZqb2(1,1)
     .            +cdabs(Q(-j)*q1+L(-j)*r1*prop)**2*AqbgZqb2(1,2)
     .            +cdabs(Q(-j)*q1+R(-j)*l1*prop)**2*AqbgZqb2(2,1)
     .            +cdabs(Q(-j)*q1+R(-j)*r1*prop)**2*AqbgZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=cdabs(Q(k)*q1+L(k)*l1*prop)**2*AgqZq2(1,1)
     .            +cdabs(Q(k)*q1+L(k)*r1*prop)**2*AgqZq2(1,2)
     .            +cdabs(Q(k)*q1+R(k)*l1*prop)**2*AgqZq2(2,1)
     .            +cdabs(Q(k)*q1+R(k)*r1*prop)**2*AgqZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=cdabs(Q(-k)*q1+L(-k)*l1*prop)**2*AgqbZqb2(1,1)
     .            +cdabs(Q(-k)*q1+L(-k)*r1*prop)**2*AgqbZqb2(1,2)
     .            +cdabs(Q(-k)*q1+R(-k)*l1*prop)**2*AgqbZqb2(2,1)
     .            +cdabs(Q(-k)*q1+R(-k)*r1*prop)**2*AgqbZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      endif

 20   continue
      enddo
      enddo
      return
      end


      subroutine zgamps2(j1,j2,j3,j4,j5,za,zb,amps2)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex amps(2,2,2)
      double precision amps2(2,2)
      integer hq,hl,hg,j1,j2,j3,j4,j5
c-- amplitude helicities are amps(quark,lepton,gluon)

      amps(1,1,1)=za(j2,j3)/za(j1,j5)/za(j2,j5)
     .             *(za(j2,j1)*zb(j4,j1)+za(j2,j5)*zb(j4,j5))

      amps(1,1,2)=zb(j4,j1)/zb(j1,j5)/zb(j2,j5)
     .             *(za(j2,j3)*zb(j2,j1)+za(j3,j5)*zb(j1,j5))

      amps(1,2,1)=za(j2,j4)/za(j1,j5)/za(j2,j5)
     .             *(za(j2,j1)*zb(j3,j1)+za(j2,j5)*zb(j3,j5))

      amps(1,2,2)=zb(j3,j1)/zb(j1,j5)/zb(j2,j5)
     .             *(za(j2,j4)*zb(j2,j1)+za(j4,j5)*zb(j1,j5))

      do hl=1,2
      do hg=1,2
        amps(2,hl,hg)=-dconjg(amps(1,3-hl,3-hg))
      enddo
      enddo

      do hq=1,2
      do hl=1,2
        amps2(hq,hl)=cdabs(amps(hq,hl,1))**2+cdabs(amps(hq,hl,2))**2
      enddo
      enddo


      return
      end

