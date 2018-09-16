      subroutine gQ_zQ_g(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*    Authors: R.K. Ellis and John Campbell                             *
*    July, 2003.                                                       *
*    Matrix element for Z + heavy quark (of flavour "flav") production *
*    with gluon radiation (order alpha_s^3)                            *
*    averaged over initial colours and spins                           *
*     g(-p1)+Q(-p2)-->Z^+(l(p3)+a(p4))+Q(p5)+g(p6)                     *
*                                                                      *
*    isub=1 : particle 6 is a gluon or light quark                     *
*    isub=2 : particle 6 is also a heavy quark                         *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'heavyflav.f'
      integer:: j,k,pq,pl,j1,j2,j3,isub
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,faclo,
     & qgZqg2(2,2),qbgZqbg2(2,2),gqbZqbg2(2,2),gqZqg2(2,2),ggZqbq2(2,2)
      real(dp):: mqq(fn:nf,fn:nf)
C      real(dp):: tanni
      complex(dp):: a111,a112,a121,a211,a122,a212,a221,a222,prop
      complex(dp):: b111,b112,b121,b211,b122,b212,b221,b222

      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: Rqb_a(2,2,2),Rqb_b(2,2,2)

      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)
      complex(dp):: Rq_a(2,2,2),Rq_b(2,2,2)

      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)
      complex(dp):: Rbq_a(2,2,2),Rbq_b(2,2,2)

      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: Rbqb_a(2,2,2),Rbqb_b(2,2,2)

      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2)

c      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)

      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)

      common/isub/isub

      integer,parameter::swap(2)=(/2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      mqq(j,k)=0._dp
      enddo
      enddo

      call spinoru(6,p,za,zb)

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

c--- calculate 2-quark, 2-gluon amplitudes
      if (isub == 1) then
        call z2jetsq(1,5,3,4,2,6,za,zb,qgZqg2)
        call z2jetsq(2,5,3,4,1,6,za,zb,gqZqg2)

        do j=1,2
        do k=1,2
        qbgZqbg2(j,k)=qgZqg2(swap(j),k)
        gqbZqbg2(j,k)=gqZqg2(swap(j),k)
        enddo
        enddo
      endif

      if (isub == 2) then
c--- NB this is the matrix element for gg->Z qb(5) q(6)
        call z2jetsq(5,6,3,4,1,2,za,zb,ggZqbq2)
      endif

c--- Add the overall factor and average over colour, spins
      fac=v*xn/four*(esq*gsq)**2
      do pq=1,2
      do pl=1,2
        if (isub == 1) then
          gqZqg2(pq,pl)  = aveqg*fac*gqZqg2(pq,pl)
          qgZqg2(pq,pl)  = aveqg*fac*qgZqg2(pq,pl)
          gqbZqbg2(pq,pl)= aveqg*fac*gqbZqbg2(pq,pl)
          qbgZqbg2(pq,pl)= aveqg*fac*qbgZqbg2(pq,pl)
        endif
        if (isub == 2) then
          ggZqbq2(pq,pl) = avegg*fac*ggZqbq2(pq,pl)
        endif
      enddo
      enddo

c--- now calculate the four quark amplitudes

      if ((isub == 1) .or. (isub == 2)) then
c--- qR->qR (i.e. Q+q -> Q+q)
        call ampqqb_qqb(1,5,6,2,qR_a,qR_b)
c--- Rq->qR (i.e. q+Q -> Q+q)
        call ampqqb_qqb(1,6,5,2,Rq_a,Rq_b)

c--- qbRb->qbRb (i.e. Qb+qb -> Qb+qb)
        call ampqqb_qqb(5,1,2,6,qbRb_a,qbRb_b)
c--- Rbqb->qbRb (i.e. qb+Qb -> Qb+qb)
        call ampqqb_qqb(6,1,2,5,Rbqb_a,Rbqb_b)

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
      endif

      if (isub == 1) then
c--- qRb->qRb (i.e. Q+qb -> Q+qb)
        call ampqqb_qqb(1,5,2,6,qRb_a,qRb_b)
c--- Rqb->qbR (i.e. q+Qb -> Qb+q)
        call ampqqb_qqb(1,6,2,5,Rqb_a,Rqb_b)

c--- qR->qR
c instead of calling ampqqb_qqb(1,5,6,2,qR_a,qR_b)
c      do j1=1,2
c      do j2=1,2
c      do j3=1,2
c      qR_a(j1,j2,j3)=+qRb_a(j1,swap(j2),j3)
c      qR_b(j1,j2,j3)=-qRb_b(j1,swap(j2),j3)
c      enddo
c      enddo
c      enddo

c--- qbR->qbR (i.e. Qb+q -> Qb+q)
        call ampqqb_qqb(5,1,6,2,qbR_a,qbR_b)
c--- Rbq->qRb (i.e. qb+Q -> Q+qb)
        call ampqqb_qqb(6,1,5,2,Rbq_a,Rbq_b)

c--- qbRb->qbRb
c instead of calling ampqqb_qqb(5,1,2,6,qbRb_a,qbRb_b)
c      do j1=1,2
c      do j2=1,2
c      do j3=1,2
c      qbRb_a(j1,j2,j3)=-qRb_a(swap(j1),j2,j3)
c      qbRb_b(j1,j2,j3)=+qRb_b(swap(j1),j2,j3)
c      enddo
c      enddo
c      enddo

c--- qqb->qqb
        call ampqqb_qqb(1,2,5,6,qqb_a,qqb_b)

c      call ampqqb_qqb(2,1,6,5,qbq_a,qbq_b)
c      do j1=1,2
c      do j2=1,2
c      do j3=1,2
c        write(6,*) j1,j2,j3,qqb_b(j1,j2,j3),qbq_b(j1,j2,j3)
c        write(6,*) j1,j2,j3,qbq_a(j1,j2,j3)-
c     &             (-qqb_a(swap(j1),swap(j2),j3))
c        write(6,*) j1,j2,j3,qbq_b(j1,j2,j3)-
c     &             (-qqb_b(swap(j1),swap(j2),j3))
c      enddo
c      enddo
c      enddo
c      pause


c--- qbq->qbq
c instead of calling ampqqb_qqb(2,1,6,5,qbq_a,qbq_b)
c        do j1=1,2
c        do j2=1,2
c        do j3=1,2
c        qbq_a(j1,j2,j3)=-qqb_a(swap(j1),swap(j2),j3)
c        qbq_b(j1,j2,j3)=-qqb_b(swap(j1),swap(j2),j3)
c        enddo
c        enddo
c        enddo
      endif

      faclo=4._dp*V*gsq**2*esq**2*aveqq

c--- sum up the helicity amplitudes for the 2-quark, 2-gluon pieces
      do j=-flav,flav,flav
      do k=-flav,flav,flav
      if ((j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 19

      if     ((j == 0) .and. (k == 0)) then
          if (isub == 2) then
          msq(j,k)=
     &      +abs(Q(flav)*q1+L(flav)*l1*prop)**2*ggZqbq2(1,1)
     &      +abs(Q(flav)*q1+R(flav)*r1*prop)**2*ggZqbq2(2,2)
     &      +abs(Q(flav)*q1+L(flav)*r1*prop)**2*ggZqbq2(1,2)
     &      +abs(Q(flav)*q1+R(flav)*l1*prop)**2*ggZqbq2(2,1)
          endif
      elseif ((j > 0) .and. (k == 0)) then
          if (isub == 1) then
          msq(j,k)=
     &    +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg2(1,1)
     &    +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg2(2,2)
     &    +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg2(1,2)
     &    +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg2(2,1)
          endif
      elseif ((j < 0) .and. (k == 0)) then
          if (isub == 1) then
          msq(j,k)=
     &    +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg2(1,1)
     &    +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg2(2,2)
     &    +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg2(1,2)
     &    +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg2(2,1)
          endif
      elseif ((j == 0) .and. (k > 0)) then
          if (isub == 1) then
          msq(j,k)=
     &    +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg2(1,1)
     &    +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg2(2,2)
     &    +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg2(1,2)
     &    +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg2(2,1)
          endif
      elseif ((j == 0) .and. (k < 0)) then
          if (isub == 1) then
          msq(j,k)=
     &    +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg2(1,1)
     &    +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg2(2,2)
     &    +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg2(1,2)
     &    +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg2(2,1)
          endif
      endif

   19 continue
      enddo
      enddo

c--- sum up the helicity amplitudes for the 4-quark pieces
      do j=-nf,nf
      do k=-nf,nf
        if ((j > 0) .and. (k > 0)) then
c--- q-q case
          if ((j .ne. flav) .and. (k .ne. flav)) goto 20
          if (j .ne. k) then
            if (isub == 2) goto 20
            if (j == flav) then
            a111=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,1,1)
     &          +(Q(k)*q1+L(k)*l1*prop)*qR_b(1,1,1)
            a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     &          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
            a112=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,1,2)
     &          +(Q(k)*q1+L(k)*r1*prop)*qR_b(1,1,2)
            a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     &          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
            a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     &          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
            a221=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,2,1)
     &          +(Q(k)*q1+R(k)*l1*prop)*qR_b(2,2,1)
            a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     &          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
            a222=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,2,2)
     &          +(Q(k)*q1+R(k)*r1*prop)*qR_b(2,2,2)
            else
            a111=(Q(j)*q1+L(j)*l1*prop)*Rq_a(1,1,1)
     &          +(Q(k)*q1+L(k)*l1*prop)*Rq_b(1,1,1)
            a121=(Q(j)*q1+L(j)*l1*prop)*Rq_a(1,2,1)
     &          +(Q(k)*q1+R(k)*l1*prop)*Rq_b(1,2,1)
            a112=(Q(j)*q1+L(j)*r1*prop)*Rq_a(1,1,2)
     &          +(Q(k)*q1+L(k)*r1*prop)*Rq_b(1,1,2)
            a122=(Q(j)*q1+L(j)*r1*prop)*Rq_a(1,2,2)
     &          +(Q(k)*q1+R(k)*r1*prop)*Rq_b(1,2,2)
            a211=(Q(j)*q1+R(j)*l1*prop)*Rq_a(2,1,1)
     &          +(Q(k)*q1+L(k)*l1*prop)*Rq_b(2,1,1)
            a221=(Q(j)*q1+R(j)*l1*prop)*Rq_a(2,2,1)
     &          +(Q(k)*q1+R(k)*l1*prop)*Rq_b(2,2,1)
            a212=(Q(j)*q1+R(j)*r1*prop)*Rq_a(2,1,2)
     &          +(Q(k)*q1+L(k)*r1*prop)*Rq_b(2,1,2)
            a222=(Q(j)*q1+R(j)*r1*prop)*Rq_a(2,2,2)
     &          +(Q(k)*q1+R(k)*r1*prop)*Rq_b(2,2,2)
            endif

            mqq(j,k)=
     &      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

          elseif (j == k) then
c--- don't include the Q+Q->Z+Q+Q case here
c            goto 20
            if (isub == 1) goto 20
            a111=(Q(j)*q1+L(j)*l1*prop)*(qR_a(1,1,1)+qR_b(1,1,1))
            b111=(Q(j)*q1+L(j)*l1*prop)*(qq_a(1,1,1)+qq_b(1,1,1))
            a112=(Q(j)*q1+L(j)*r1*prop)*(qR_a(1,1,2)+qR_b(1,1,2))
            b112=(Q(j)*q1+L(j)*r1*prop)*(qq_a(1,1,2)+qq_b(1,1,2))
            a221=(Q(j)*q1+R(j)*l1*prop)*(qR_a(2,2,1)+qR_b(2,2,1))
            b221=(Q(j)*q1+R(j)*l1*prop)*(qq_a(2,2,1)+qq_b(2,2,1))
            a222=(Q(j)*q1+R(j)*r1*prop)*(qR_a(2,2,2)+qR_b(2,2,2))
            b222=(Q(j)*q1+R(j)*r1*prop)*(qq_a(2,2,2)+qq_b(2,2,2))

            a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     &          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
            b121=(Q(j)*q1+L(j)*l1*prop)*qq_a(1,2,1)
     &          +(Q(k)*q1+R(k)*l1*prop)*qq_b(1,2,1)
            a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     &          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
            b122=(Q(j)*q1+L(j)*r1*prop)*qq_a(1,2,2)
     &          +(Q(k)*q1+R(k)*r1*prop)*qq_b(1,2,2)
            a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     &          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
            b211=(Q(j)*q1+R(j)*l1*prop)*qq_a(2,1,1)
     &          +(Q(k)*q1+L(k)*l1*prop)*qq_b(2,1,1)
            a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     &          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
            b212=(Q(j)*q1+R(j)*r1*prop)*qq_a(2,1,2)
     &          +(Q(k)*q1+L(k)*r1*prop)*qq_b(2,1,2)

            mqq(j,k)=half*faclo*(
     &      +Dble(a111*conjg(b111))+Dble(a112*conjg(b112))
     &      +Dble(a221*conjg(b221))+Dble(a222*conjg(b222)))*two/xn
     &              +half*faclo*
     &      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
     &              +half*faclo*(
     &      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     &      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
          endif
        elseif ((j < 0) .and. (k < 0)) then
c--- qb-qb case
          if ((j .ne. -flav) .and. (k .ne. -flav)) goto 20
          if (j .ne. k) then
            if (isub == 2) goto 20
            if (j == -flav) then
            a111=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(1,1,1)
            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)

            a112=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(1,1,2)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)

            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
            a221=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(2,2,1)

            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)
            a222=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(2,2,2)
            else
            a111=(Q(-j)*q1+L(-j)*l1*prop)*Rbqb_a(1,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*Rbqb_b(1,1,1)
            a121=(Q(-j)*q1+L(-j)*l1*prop)*Rbqb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*Rbqb_b(1,2,1)

            a112=(Q(-j)*q1+L(-j)*r1*prop)*Rbqb_a(1,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*Rbqb_b(1,1,2)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*Rbqb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*Rbqb_b(1,2,2)

            a211=(Q(-j)*q1+R(-j)*l1*prop)*Rbqb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*Rbqb_b(2,1,1)
            a221=(Q(-j)*q1+R(-j)*l1*prop)*Rbqb_a(2,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*Rbqb_b(2,2,1)

            a212=(Q(-j)*q1+R(-j)*r1*prop)*Rbqb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*Rbqb_b(2,1,2)
            a222=(Q(-j)*q1+R(-j)*r1*prop)*Rbqb_a(2,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*Rbqb_b(2,2,2)
            endif

            mqq(j,k)=
     &      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
          elseif (j == k) then
c--- don't include the Qb+Qb->Z+Qb+Qb case here
c            goto 20
            if (isub == 1) goto 20

            a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbRb_a(1,1,1)+qbRb_b(1,1,1))
            b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbqb_a(1,1,1)+qbqb_b(1,1,1))
            a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbRb_a(1,1,2)+qbRb_b(1,1,2))
            b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbqb_a(1,1,2)+qbqb_b(1,1,2))
            a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbRb_a(2,2,1)+qbRb_b(2,2,1))
            b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbqb_a(2,2,1)+qbqb_b(2,2,1))
            a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbRb_a(2,2,2)+qbRb_b(2,2,2))
            b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbqb_a(2,2,2)+qbqb_b(2,2,2))


            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)

            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbqb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qbqb_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbqb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qbqb_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbqb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qbqb_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbqb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qbqb_b(2,1,2)


            mqq(j,k)=half*faclo*(
     &      +Dble(a111*conjg(b111))+Dble(a112*conjg(b112))
     &      +Dble(a221*conjg(b221))+Dble(a222*conjg(b222)))*two/xn
     &              +half*faclo*
     &      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
     &              +half*faclo*(
     &      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     &      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
          endif

        elseif ((j > 0) .and. (k < 0)) then
c--- q-qb case
          if (j .ne. -k) then
            if ((j .ne. flav) .and. (k .ne. -flav)) goto 20
            if (isub == 2) goto 20

            if (j == flav) then
            a111=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(1,1,1)
            a112=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(1,1,2)
            a221=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(2,2,1)
            a222=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(2,2,2)

            a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
            a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
            a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
            a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)
            else
            a111=(Q(+j)*q1+L(+j)*l1*prop)*Rqb_a(1,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*Rqb_b(1,1,1)
            a112=(Q(+j)*q1+L(+j)*r1*prop)*Rqb_a(1,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*Rqb_b(1,1,2)
            a221=(Q(+j)*q1+R(+j)*l1*prop)*Rqb_a(2,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*Rqb_b(2,2,1)
            a222=(Q(+j)*q1+R(+j)*r1*prop)*Rqb_a(2,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*Rqb_b(2,2,2)

            a121=(Q(+j)*q1+L(+j)*l1*prop)*Rqb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*Rqb_b(1,2,1)
            a122=(Q(+j)*q1+L(+j)*r1*prop)*Rqb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*Rqb_b(1,2,2)
            a211=(Q(+j)*q1+R(+j)*l1*prop)*Rqb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*Rqb_b(2,1,1)
            a212=(Q(+j)*q1+R(+j)*r1*prop)*Rqb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*Rqb_b(2,1,2)
            endif

            mqq(j,k)=
     &      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

          elseif (j == -k) then
c--- don't include the Q+Qb->Z+Q+Qb case here
            goto 20
c--case where final state from annihilation diagrams is the same quark
c            if (s(5,6) < 4._dp*mQsq) goto 20
c--- cut here to prevent the collinear divergence
            if ((abs(j) == flav) .and. (isub == 2)) then
              a111=(Q(j)*q1+L(j)*l1*prop)*(qRb_a(1,1,1)+qRb_b(1,1,1))
c              b111=(Q(j)*q1+L(j)*l1*prop)*(qqb_a(1,1,1)+qqb_b(1,1,1))

              a112=(Q(j)*q1+L(j)*r1*prop)*(qRb_a(1,1,2)+qRb_b(1,1,2))
c              b112=(Q(j)*q1+L(j)*r1*prop)*(qqb_a(1,1,2)+qqb_b(1,1,2))

              a221=(Q(j)*q1+R(j)*l1*prop)*(qRb_a(2,2,1)+qRb_b(2,2,1))
c              b221=(Q(j)*q1+R(j)*l1*prop)*(qqb_a(2,2,1)+qqb_b(2,2,1))

              a222=(Q(j)*q1+R(j)*r1*prop)*(qRb_a(2,2,2)+qRb_b(2,2,2))
c              b222=(Q(j)*q1+R(j)*r1*prop)*(qqb_a(2,2,2)+qqb_b(2,2,2))

              a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     &            +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
              a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     &            +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
              a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     &            +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
              a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     &            +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)

c              b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
c     &            +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(1,2,1)
c              b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
c     &            +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(1,2,2)
c              b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
c     &            +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(2,1,1)
c              b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
c     &            +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(2,1,2)

c--- the only contribution here is from scattering diagrams
              mqq(j,k)=faclo*
     &        (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

c              mqq(j,k)=faclo*(
c     &        +Dble(a111*conjg(b111))+Dble(a112*conjg(b112))
c     &        +Dble(a221*conjg(b221))+Dble(a222*conjg(b222)))*two/xn
c     &                +faclo*(
c     &        +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
c     &        +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
c     &                +faclo*
c     &        (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
c     &        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
            endif

c            if (isub == 2) then
c              b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
c     &            +(Q(flav)*q1+L(flav)*l1*prop)*qqb_b(1,1,1)
c              b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
c     &            +(Q(flav)*q1+L(flav)*r1*prop)*qqb_b(1,1,2)
c              b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
c     &            +(Q(flav)*q1+R(flav)*l1*prop)*qqb_b(2,2,1)
c              b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
c     &            +(Q(flav)*q1+R(flav)*r1*prop)*qqb_b(2,2,2)
c              b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
c     &            +(Q(flav)*q1+R(flav)*l1*prop)*qqb_b(1,2,1)
c              b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
c     &            +(Q(flav)*q1+R(flav)*r1*prop)*qqb_b(1,2,2)
c              b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
c     &            +(Q(flav)*q1+L(flav)*l1*prop)*qqb_b(2,1,1)
c              b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
c     &            +(Q(flav)*q1+L(flav)*r1*prop)*qqb_b(2,1,2)
c              tanni=faclo*(
c     &              abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
c     &             +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
c              mqq(j,k)=mqq(j,k)+tanni
c            endif

          endif

        elseif ((j < 0) .and. (k > 0)) then
C---qb-q case
          if (j .ne. -k) then
            if ((j .ne. -flav) .and. (k .ne. flav)) goto 20
            if (isub == 2) goto 20

            if (k == flav) then
            a111=(Q(-j)*q1+L(-j)*l1*prop)*Rbq_a(1,1,1)
     &          +(Q(+k)*q1+L(+k)*l1*prop)*Rbq_b(1,1,1)
            a121=(Q(-j)*q1+L(-j)*l1*prop)*Rbq_a(1,2,1)
     &          +(Q(+k)*q1+R(+k)*l1*prop)*Rbq_b(1,2,1)
            a112=(Q(-j)*q1+L(-j)*r1*prop)*Rbq_a(1,1,2)
     &          +(Q(+k)*q1+L(+k)*r1*prop)*Rbq_b(1,1,2)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*Rbq_a(1,2,2)
     &          +(Q(+k)*q1+R(+k)*r1*prop)*Rbq_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*Rbq_a(2,1,1)
     &          +(Q(+k)*q1+L(+k)*l1*prop)*Rbq_b(2,1,1)
            a221=(Q(-j)*q1+R(-j)*l1*prop)*Rbq_a(2,2,1)
     &          +(Q(+k)*q1+R(+k)*l1*prop)*Rbq_b(2,2,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*Rbq_a(2,1,2)
     &          +(Q(+k)*q1+L(+k)*r1*prop)*Rbq_b(2,1,2)
            a222=(Q(-j)*q1+R(-j)*r1*prop)*Rbq_a(2,2,2)
     &          +(Q(+k)*q1+R(+k)*r1*prop)*Rbq_b(2,2,2)
            else
            a111=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,1,1)
     &          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(1,1,1)
            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     &          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
            a112=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,1,2)
     &          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(1,1,2)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     &          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     &          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
            a221=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,2,1)
     &          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(2,2,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     &          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)
            a222=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,2,2)
     &          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(2,2,2)
            endif

            mqq(j,k)=
     &      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

          elseif (j == -k) then
c--- don't include the Qb+Q->Z+Qb+Q case here
            goto 20
c            if (s(5,6) < 4._dp*mQsq) goto 20
c--- cut here to prevent the collinear divergence
            if ((abs(j) == flav) .and. (isub == 2)) then
              a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbR_a(1,1,1)+qbR_b(1,1,1))
c              b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbq_a(1,1,1)+qbq_b(1,1,1))
              a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbR_a(1,1,2)+qbR_b(1,1,2))
c              b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbq_a(1,1,2)+qbq_b(1,1,2))
              a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbR_a(2,2,1)+qbR_b(2,2,1))
c              b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbq_a(2,2,1)+qbq_b(2,2,1))
              a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbR_a(2,2,2)+qbR_b(2,2,2))
c              b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbq_a(2,2,2)+qbq_b(2,2,2))

              a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     &            +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
              a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     &            +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
              a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     &            +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
              a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     &            +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)

c              b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
c     &            +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(1,2,1)
c              b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
c     &            +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(1,2,2)
c              b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
c     &            +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(2,1,1)
c              b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
c     &            +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(2,1,2)

c--- the only contribution here is from scattering diagrams
              mqq(j,k)=faclo*
     &        (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     &        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

c              mqq(j,k)=faclo*(
c     &        +Dble(a111*conjg(b111))+Dble(a112*conjg(b112))
c     &        +Dble(a221*conjg(b221))+Dble(a222*conjg(b222)))*two/xn
c     &                +faclo*
c     &        (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
c     &        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
c     &                +faclo*(
c     &        +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
c     &        +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
            endif

c            if (isub == 2) then
c              b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
c     &            +(Q(flav)*q1+L(flav)*l1*prop)*qbq_b(1,1,1)
c              b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
c     &            +(Q(flav)*q1+L(flav)*r1*prop)*qbq_b(1,1,2)
c              b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
c     &            +(Q(flav)*q1+R(flav)*l1*prop)*qbq_b(2,2,1)
c              b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
c     &            +(Q(flav)*q1+R(flav)*r1*prop)*qbq_b(2,2,2)
c              b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
c     &            +(Q(flav)*q1+R(flav)*l1*prop)*qbq_b(1,2,1)
c              b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
c     &            +(Q(flav)*q1+R(flav)*r1*prop)*qbq_b(1,2,2)
c              b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
c     &            +(Q(flav)*q1+L(flav)*l1*prop)*qbq_b(2,1,1)
c              b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
c     &            +(Q(flav)*q1+L(flav)*r1*prop)*qbq_b(2,1,2)
c              tanni=faclo*(
c     &              abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
c     &             +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
c              mqq(j,k)=mqq(j,k)+tanni
c            endif

          endif

        endif

   20 continue

      msq(j,k)=msq(j,k)+mqq(j,k)
      enddo
      enddo

      return
      end





