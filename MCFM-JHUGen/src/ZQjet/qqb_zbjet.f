      subroutine qqb_zbjet(p,msq)
      implicit none
      include 'types.f'

c--- J.Campbell,  7/4/05
c--- Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z + b(p5) + g(p6)
c                          |
c                          --> l(p3)+a(p4)
c
c--- all momenta are incoming
c--- Extended to include charm quark production via the variable "flav"
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'heavyflav.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_cs.f'
      include 'nflav.f'
      integer:: i,j,k,pq,pl,j1,j2,j3,icol
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,faclo,
     &   qgZqg2(2,2),gqZqg2(2,2),
     &   qgZqg2_cs(0:2,2,2),gqZqg2_cs(0:2,2,2),
     &   qbgZqbg2_cs(0:2,2,2),gqbZqbg2_cs(0:2,2,2)
      complex(dp):: a111,a112,a121,a211,a122,a212,a221,a222
c      complex(dp):: b111,b112,b121,b211,b122,b212,b221,b222
      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2),prop
      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)
      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)
      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)
      integer,parameter::swap(2)=(/2,1/),swap1(0:2)=(/0,2,1/)

c--- initialize the matrix element squared
      msq(:,:)=0._dp
      msq_cs(:,:,:)=0._dp

      call spinoru(6,p,za,zb)

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

c--- calculate the amplitudes which involve 2 heavy quarks
c--- and 2 gluons
      call z2jetsq(1,5,3,4,2,6,za,zb,qgZqg2)
      call storecsz(qgZqg2_cs)
      call z2jetsq(2,5,3,4,1,6,za,zb,gqZqg2)
      call storecsz(gqZqg2_cs)

      do j=1,2
      do k=1,2
      do i=0,2
      qbgZqbg2_cs(i,j,k)=qgZqg2_cs(swap1(i),swap(j),k)
      gqbZqbg2_cs(i,j,k)=gqZqg2_cs(swap1(i),swap(j),k)
      enddo
      enddo
      enddo

      fac=v*xn/four*(esq*gsq)**2
      do pq=1,2
      do pl=1,2
        do i=0,2
          gqZqg2_cs(i,pq,pl)  = aveqg*fac*gqZqg2_cs(i,pq,pl)
          qgZqg2_cs(i,pq,pl)  = aveqg*fac*qgZqg2_cs(i,pq,pl)
          gqbZqbg2_cs(i,pq,pl)= aveqg*fac*gqbZqbg2_cs(i,pq,pl)
          qbgZqbg2_cs(i,pq,pl)= aveqg*fac*qbgZqbg2_cs(i,pq,pl)
        enddo

        gqZqg2(pq,pl)  = gqZqg2_cs(1,pq,pl) +gqZqg2_cs(2,pq,pl)
     &                  +gqZqg2_cs(0,pq,pl)
        qgZqg2(pq,pl)  = qgZqg2_cs(1,pq,pl)  +qgZqg2_cs(2,pq,pl)
     &                  +qgZqg2_cs(0,pq,pl)
      enddo
      enddo


c--- calculate the amplitudes which involve 2 heavy quarks
c--- and 2 light quarks

c--- qRb->qRb
      call ampqqb_qqb(1,5,2,6,qRb_a,qRb_b)
c--- qR->qR
c instead of calling ampqqb_qqb(1,5,6,2,qR_a,qR_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qR_a(j1,j2,j3)=+qRb_a(j1,swap(j2),j3)
      qR_b(j1,j2,j3)=-qRb_b(j1,swap(j2),j3)
      enddo
      enddo
      enddo

c--- qbRb->qbRb
c instead of calling ampqqb_qqb(5,1,2,6,qbRb_a,qbRb_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbRb_a(j1,j2,j3)=-qRb_a(swap(j1),j2,j3)
      qbRb_b(j1,j2,j3)=+qRb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

c--- qqb->qqb
c-- changed from 1,2,5,6
      call ampqqb_qqb(1,6,2,5,qqb_a,qqb_b)
c--- qbq->qqb
c instead of calling ampqqb_qqb(2,1,5,6,qbq_a,qbq_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbq_a(j1,j2,j3)=-qqb_a(swap(j1),j2,j3)
      qbq_b(j1,j2,j3)=+qqb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

c--- qq->qq
      call ampqqb_qqb(1,6,5,2,qq_a,qq_b)
c--- qbqb->qbqb
c instead of calling ampqqb_qqb(6,1,2,5,qbqb_a,qbqb_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbqb_a(j1,j2,j3)=-qq_a(swap(j1),swap(j2),j3)
      qbqb_b(j1,j2,j3)=-qq_b(swap(j1),swap(j2),j3)
      enddo
      enddo
      enddo

c--- qbR->qbR
      call ampqqb_qqb(6,1,5,2,qbR_a,qbR_b)
c--- qbq->qbq
      call ampqqb_qqb(5,1,6,2,qbq_a,qbq_b)

      faclo=4._dp*V*gsq**2*esq**2*aveqq

c--- add up the various contributions
      do j=-nflav,nflav
      do k=-nflav,nflav

      if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 99
      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 99
c--- so that either abs(j) or abs(k) = flav (but not both).

c--- Q-G contribution
      if     ((j == +flav) .and. (k == 0)) then
        do icol=0,2
           msq_cs(icol,j,k)=
     &     +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg2_cs(icol,1,1)
     &     +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg2_cs(icol,2,2)
     &     +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg2_cs(icol,1,2)
     &     +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg2_cs(icol,2,1)
        enddo
c--- Qb-G contribution
      elseif ((j == -flav) .and. (k == 0)) then
        do icol=0,2
           msq_cs(icol,j,k)=
     &     +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg2_cs(icol,1,1)
     &     +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg2_cs(icol,2,2)
     &     +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg2_cs(icol,1,2)
     &     +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg2_cs(icol,2,1)
        enddo
c--- G-Q contribution
      elseif ((j == 0) .and. (k == +flav)) then
        do icol=0,2
           msq_cs(icol,j,k)=
     &     +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg2_cs(icol,1,1)
     &     +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg2_cs(icol,2,2)
     &     +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg2_cs(icol,1,2)
     &     +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg2_cs(icol,2,1)
        enddo
c--- G-Qb contribution
      elseif ((j == 0) .and. (k == -flav)) then
        do icol=0,2
           msq_cs(icol,j,k)=
     &     +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg2_cs(icol,1,1)
     &     +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg2_cs(icol,2,2)
     &     +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg2_cs(icol,1,2)
     &     +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg2_cs(icol,2,1)
        enddo
c--- Q-Q contribution
      elseif ((j > 0) .and. (k > 0)) then
        if     (j == +flav) then
          a111=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,1,1)
     &        +(Q(k)*q1+L(k)*l1*prop)*qR_b(1,1,1)
          a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     &        +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
          a112=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,1,2)
     &        +(Q(k)*q1+L(k)*r1*prop)*qR_b(1,1,2)
          a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     &        +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
          a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     &        +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
          a221=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,2,1)
     &        +(Q(k)*q1+R(k)*l1*prop)*qR_b(2,2,1)
          a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     &        +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
          a222=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,2,2)
     &        +(Q(k)*q1+R(k)*r1*prop)*qR_b(2,2,2)
        elseif (k == +flav) then
          a111=(Q(j)*q1+L(j)*l1*prop)*qq_a(1,1,1)
     &        +(Q(k)*q1+L(k)*l1*prop)*qq_b(1,1,1)
          a121=(Q(j)*q1+L(j)*l1*prop)*qq_a(1,2,1)
     &        +(Q(k)*q1+R(k)*l1*prop)*qq_b(1,2,1)
          a112=(Q(j)*q1+L(j)*r1*prop)*qq_a(1,1,2)
     &        +(Q(k)*q1+L(k)*r1*prop)*qq_b(1,1,2)
          a122=(Q(j)*q1+L(j)*r1*prop)*qq_a(1,2,2)
     &        +(Q(k)*q1+R(k)*r1*prop)*qq_b(1,2,2)
          a211=(Q(j)*q1+R(j)*l1*prop)*qq_a(2,1,1)
     &        +(Q(k)*q1+L(k)*l1*prop)*qq_b(2,1,1)
          a221=(Q(j)*q1+R(j)*l1*prop)*qq_a(2,2,1)
     &        +(Q(k)*q1+R(k)*l1*prop)*qq_b(2,2,1)
          a212=(Q(j)*q1+R(j)*r1*prop)*qq_a(2,1,2)
     &        +(Q(k)*q1+L(k)*r1*prop)*qq_b(2,1,2)
          a222=(Q(j)*q1+R(j)*r1*prop)*qq_a(2,2,2)
     &        +(Q(k)*q1+R(k)*r1*prop)*qq_b(2,2,2)
        endif
        msq_cs(0,j,k)=zip
        msq_cs(1,j,k)=
     &  +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &         +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
        msq_cs(2,j,k)=zip
c--- Qb-Qb contribution
      elseif ((j < 0) .and. (k < 0)) then
        if     (j == -flav) then
          a111=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(1,1,1)
          a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)

          a112=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(1,1,2)
          a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)

          a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
          a221=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(2,2,1)

          a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)
          a222=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(2,2,2)
        elseif (k == -flav) then
          a111=(Q(-j)*q1+L(-j)*l1*prop)*qbqb_a(1,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qbqb_b(1,1,1)
          a121=(Q(-j)*q1+L(-j)*l1*prop)*qbqb_a(1,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qbqb_b(1,2,1)

          a112=(Q(-j)*q1+L(-j)*r1*prop)*qbqb_a(1,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qbqb_b(1,1,2)
          a122=(Q(-j)*q1+L(-j)*r1*prop)*qbqb_a(1,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qbqb_b(1,2,2)

          a211=(Q(-j)*q1+R(-j)*l1*prop)*qbqb_a(2,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qbqb_b(2,1,1)
          a221=(Q(-j)*q1+R(-j)*l1*prop)*qbqb_a(2,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qbqb_b(2,2,1)

          a212=(Q(-j)*q1+R(-j)*r1*prop)*qbqb_a(2,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qbqb_b(2,1,2)
          a222=(Q(-j)*q1+R(-j)*r1*prop)*qbqb_a(2,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qbqb_b(2,2,2)
        endif
        msq_cs(0,j,k)=zip
        msq_cs(1,j,k)=
     &  +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &         +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
        msq_cs(2,j,k)=zip
c--- Q-Qb contribution
      elseif ((j > 0) .and. (k < 0)) then
        if     (j == +flav) then
          a111=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(1,1,1)
          a112=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(1,1,2)
          a221=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(2,2,1)
          a222=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(2,2,2)

          a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
          a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
          a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
          a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)
        elseif (k == -flav) then
          a111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(1,1,1)
          a112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(1,1,2)
          a221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(2,2,1)
          a222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(2,2,2)

          a121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     &        +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(1,2,1)
          a122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     &        +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(1,2,2)
          a211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     &        +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(2,1,1)
          a212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     &        +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(2,1,2)
        endif
        msq_cs(0,j,k)=zip
        msq_cs(1,j,k)=zip
        msq_cs(2,j,k)=
     &  +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &         +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
      elseif ((j < 0) .and. (k > 0)) then
c--- Qb-Q contribution
        if     (j == -flav) then
          a111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     &        +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(1,1,1)
          a121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     &        +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(1,2,1)
          a112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     &        +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(1,1,2)
          a122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     &        +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(1,2,2)
          a211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     &        +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(2,1,1)
          a221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     &        +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(2,2,1)
          a212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     &        +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(2,1,2)
          a222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     &        +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(2,2,2)
        elseif (k == +flav) then
          a111=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,1,1)
     &        +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(1,1,1)
          a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     &        +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
          a112=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,1,2)
     &        +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(1,1,2)
          a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     &        +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
          a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     &        +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
          a221=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,2,1)
     &        +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(2,2,1)
          a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     &        +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)
          a222=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,2,2)
     &        +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(2,2,2)
        endif
        msq_cs(0,j,k)=zip
        msq_cs(1,j,k)=zip
        msq_cs(2,j,k)=
     &  +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &         +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
      endif

c--- only the colour-ordered matrix elements are calculated above, so
c--- we have to calculate the total for "msq"
      msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)

   99 continue
      enddo
      enddo

      return
      end





