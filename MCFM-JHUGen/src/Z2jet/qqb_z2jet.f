      subroutine qqb_z2jet(p,msq)
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)
c
c--all momenta incoming
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_cs.f'
      include 'flags.f'
      include 'nflav.f'
      integer:: i,j,k,pq,pl,nquark,nup,ndo,j1,j2,j3,icol
      integer,parameter::swap(2)=(/2,1/),swap1(0:2)=(/0,2,1/)
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,faclo,
     &   qqbZgg2(2,2),qgZqg2(2,2),
c    .   qbqZgg2(2,2),qbgZqbg2(2,2),gqbZqbg2(2,2),
     &   gqZqg2(2,2),ggZqbq2(2,2),
     &   qqbZgg2_cs(0:2,2,2),qbqZgg2_cs(0:2,2,2),
     &   qgZqg2_cs(0:2,2,2),gqZqg2_cs(0:2,2,2),
     &   qbgZqbg2_cs(0:2,2,2),gqbZqbg2_cs(0:2,2,2),
     &   ggZqbq2_cs(0:2,2,2),ggtemp(0:2)
      real(dp):: tup,tdo
      complex(dp):: a111,a112,a121,a211,a122,a212,a221,a222
      complex(dp):: b111,b112,b121,b211,b122,b212,b221,b222

      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2),prop

      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)

      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)

      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)

      real(dp):: coupqe(nf,2,2)

      real(dp):: mqq(0:2,fn:nf,fn:nf)
      common/mqq/mqq
!$omp threadprivate(/mqq/)
      include 'cplx.h'

      msq(:,:)=zip

      call spinoru(6,p,za,zb)

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

c--- calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then
        call z2jetsq(1,2,3,4,5,6,za,zb,qqbZgg2)
        call storecsz(qqbZgg2_cs)
        call z2jetsq(1,5,3,4,2,6,za,zb,qgZqg2)
        call storecsz(qgZqg2_cs)
        call z2jetsq(2,5,3,4,1,6,za,zb,gqZqg2)
        call storecsz(gqZqg2_cs)

        do j=1,2
        do k=1,2
c        qbqZgg2(j,k)=qqbZgg2(swap(j),k)
c        qbgZqbg2(j,k)=qgZqg2(swap(j),k)
c        gqbZqbg2(j,k)=gqZqg2(swap(j),k)
        do i=0,2
        qbqZgg2_cs(i,j,k)=qqbZgg2_cs(swap1(i),swap(j),k)
        qbgZqbg2_cs(i,j,k)=qgZqg2_cs(swap1(i),swap(j),k)
        gqbZqbg2_cs(i,j,k)=gqZqg2_cs(swap1(i),swap(j),k)
        enddo
        enddo
        enddo

c        call z2jetsq(2,1,3,4,5,6,za,zb,qbqZgg2)
c        call storecsz(qbqZgg2_cs)
c        call z2jetsq(5,1,3,4,2,6,za,zb,qbgZqbg2)
c        call storecsz(qbgZqbg2_cs)
c        call z2jetsq(5,2,3,4,1,6,za,zb,gqbZqbg2)
c        call storecsz(gqbZqbg2_cs)

C --NB this is the matrix element for gg->Z qb(5) q(6)
        call z2jetsq(5,6,3,4,1,2,za,zb,ggZqbq2)
        call storecsz(ggZqbq2_cs)

        fac=v*xn/four*(esq*gsq)**2
        do pq=1,2
        do pl=1,2
        do i=0,2
          qqbZgg2_cs(i,pq,pl) = half*aveqq*fac*qqbZgg2_cs(i,pq,pl)
          qbqZgg2_cs(i,pq,pl) = half*aveqq*fac*qbqZgg2_cs(i,pq,pl)
          gqZqg2_cs(i,pq,pl)  = aveqg*fac*gqZqg2_cs(i,pq,pl)
          qgZqg2_cs(i,pq,pl)  = aveqg*fac*qgZqg2_cs(i,pq,pl)
          gqbZqbg2_cs(i,pq,pl)= aveqg*fac*gqbZqbg2_cs(i,pq,pl)
          qbgZqbg2_cs(i,pq,pl)= aveqg*fac*qbgZqbg2_cs(i,pq,pl)
          ggZqbq2_cs(i,pq,pl) = avegg*fac*ggZqbq2_cs(i,pq,pl)
       enddo

        qqbZgg2(pq,pl) = qqbZgg2_cs(1,pq,pl)+qqbZgg2_cs(2,pq,pl)
     &                  +qqbZgg2_cs(0,pq,pl)
        gqZqg2(pq,pl)  = gqZqg2_cs(1,pq,pl) +gqZqg2_cs(2,pq,pl)
     &                  +gqZqg2_cs(0,pq,pl)
        qgZqg2(pq,pl)  = qgZqg2_cs(1,pq,pl)  +qgZqg2_cs(2,pq,pl)
     &                  +qgZqg2_cs(0,pq,pl)
c        qbqZgg2(pq,pl) = qbqZgg2_cs(1,pq,pl)+qbqZgg2_cs(2,pq,pl)
c     &                  +qbqZgg2_cs(0,pq,pl)
c        gqbZqbg2(pq,pl)= gqbZqbg2_cs(1,pq,pl)+gqbZqbg2_cs(2,pq,pl)
c     &                  +gqbZqbg2_cs(0,pq,pl)
c        qbgZqbg2(pq,pl)= qbgZqbg2_cs(1,pq,pl)+qbgZqbg2_cs(2,pq,pl)
c     &                  +qbgZqbg2_cs(0,pq,pl)
        ggZqbq2(pq,pl) = ggZqbq2_cs(1,pq,pl) +ggZqbq2_cs(2,pq,pl)
     &                  +ggZqbq2_cs(0,pq,pl)
        enddo
        enddo
      endif


      if (Qflag) then
      call spinoru(6,p,za,zb)

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
c--- qbR->qbR
      call ampqqb_qqb(6,1,5,2,qbR_a,qbR_b)

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
      call ampqqb_qqb(1,2,5,6,qqb_a,qqb_b)
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

      faclo=four*V*gsq**2*esq**2*aveqq
      endif



      if (Gflag) then

c--- initialize couplings of quarks to leptons
      do nquark=1,nf
        coupqe(nquark,1,1)=real((Q(nquark)*q1+L(nquark)*l1*prop)
     &                    *conjg(Q(nquark)*q1+L(nquark)*l1*prop),dp)
        coupqe(nquark,2,2)=real((Q(nquark)*q1+R(nquark)*r1*prop)
     &                    *conjg(Q(nquark)*q1+R(nquark)*r1*prop),dp)
        coupqe(nquark,1,2)=real((Q(nquark)*q1+L(nquark)*r1*prop)
     &                    *conjg(Q(nquark)*q1+L(nquark)*r1*prop),dp)
        coupqe(nquark,2,1)=real((Q(nquark)*q1+R(nquark)*l1*prop)
     &                    *conjg(Q(nquark)*q1+R(nquark)*l1*prop),dp)
      enddo

      do j=-nf,nf
      do k=-nf,nf

      do icol=0,2
      msq_cs(icol,j,k)=zip
      enddo

      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then

          do icol=0,2
          ggtemp(icol)=zip
          do nquark=1,nflav
           ggtemp(icol)=ggtemp(icol)
     &      +coupqe(nquark,1,1)*ggZqbq2_cs(icol,1,1)
     &      +coupqe(nquark,2,2)*ggZqbq2_cs(icol,2,2)
     &      +coupqe(nquark,1,2)*ggZqbq2_cs(icol,1,2)
     &      +coupqe(nquark,2,1)*ggZqbq2_cs(icol,2,1)
          enddo
          msq_cs(icol,j,k)=ggtemp(icol)
          enddo
      elseif ((j > 0) .and. (k < 0)) then
          do icol=0,2
             msq_cs(icol,j,k)=
     &       +coupqe(j,1,1)*qqbZgg2_cs(icol,1,1)
     &       +coupqe(j,2,2)*qqbZgg2_cs(icol,2,2)
     &       +coupqe(j,1,2)*qqbZgg2_cs(icol,1,2)
     &       +coupqe(j,2,1)*qqbZgg2_cs(icol,2,1)
          enddo
c---Statistical factor already included above
      elseif ((j < 0) .and. (k > 0)) then
          do icol=0,2
             msq_cs(icol,j,k)=
     &       +coupqe(k,1,1)*qbqZgg2_cs(icol,1,1)
     &       +coupqe(k,2,2)*qbqZgg2_cs(icol,2,2)
     &       +coupqe(k,1,2)*qbqZgg2_cs(icol,1,2)
     &       +coupqe(k,2,1)*qbqZgg2_cs(icol,2,1)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          do icol=0,2
             msq_cs(icol,j,k)=
     &       +coupqe(j,1,1)*qgZqg2_cs(icol,1,1)
     &       +coupqe(j,2,2)*qgZqg2_cs(icol,2,2)
     &       +coupqe(j,1,2)*qgZqg2_cs(icol,1,2)
     &       +coupqe(j,2,1)*qgZqg2_cs(icol,2,1)
          enddo
      elseif ((j < 0) .and. (k == 0)) then
          do icol=0,2
             msq_cs(icol,j,k)=
     &       +coupqe(-j,1,1)*qbgZqbg2_cs(icol,1,1)
     &       +coupqe(-j,2,2)*qbgZqbg2_cs(icol,2,2)
     &       +coupqe(-j,1,2)*qbgZqbg2_cs(icol,1,2)
     &       +coupqe(-j,2,1)*qbgZqbg2_cs(icol,2,1)
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          do icol=0,2
             msq_cs(icol,j,k)=
     &       +coupqe(k,1,1)*gqZqg2_cs(icol,1,1)
     &       +coupqe(k,2,2)*gqZqg2_cs(icol,2,2)
     &       +coupqe(k,1,2)*gqZqg2_cs(icol,1,2)
     &       +coupqe(k,2,1)*gqZqg2_cs(icol,2,1)
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          do icol=0,2
             msq_cs(icol,j,k)=
     &       +coupqe(-k,1,1)*gqbZqbg2_cs(icol,1,1)
     &       +coupqe(-k,2,2)*gqbZqbg2_cs(icol,2,2)
     &       +coupqe(-k,1,2)*gqbZqbg2_cs(icol,1,2)
     &       +coupqe(-k,2,1)*gqbZqbg2_cs(icol,2,1)
          enddo
      endif
      msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)

   19 continue
      enddo
      enddo
      endif

      if (Qflag) then

      do j=-nf,nf
      do k=-nf,nf

      do icol=0,2
      mqq(icol,j,k)=zip
      enddo

          if ((j > 0) .and. (k > 0)) then
c----QQ case
            if (j .ne. k) then
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
            mqq(0,j,k)=zip
            mqq(1,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=zip
            elseif (j == k) then
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

            mqq(0,j,k)=half*faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(1,j,k)=half*faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=half*faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

            endif
          elseif ((j < 0) .and. (k < 0)) then
c----QbQb case
            if (j .ne. k) then
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
            mqq(0,j,k)=zip
            mqq(1,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=zip
            elseif (j == k) then

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


            mqq(0,j,k)=half*faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(1,j,k)=half*faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=half*faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

            endif
C---q-qb case
         elseif ((j > 0) .and. (k < 0)) then
             if (j .ne. -k) then
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
            mqq(0,j,k)=zip
            mqq(1,j,k)=zip
            mqq(2,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))

            elseif (j == -k) then
c--case where final state from annihilation diagrams is the same quark
            a111=(Q(j)*q1+L(j)*l1*prop)*(qRb_a(1,1,1)+qRb_b(1,1,1))
            b111=(Q(j)*q1+L(j)*l1*prop)*(qqb_a(1,1,1)+qqb_b(1,1,1))

            a112=(Q(j)*q1+L(j)*r1*prop)*(qRb_a(1,1,2)+qRb_b(1,1,2))
            b112=(Q(j)*q1+L(j)*r1*prop)*(qqb_a(1,1,2)+qqb_b(1,1,2))

            a221=(Q(j)*q1+R(j)*l1*prop)*(qRb_a(2,2,1)+qRb_b(2,2,1))
            b221=(Q(j)*q1+R(j)*l1*prop)*(qqb_a(2,2,1)+qqb_b(2,2,1))

            a222=(Q(j)*q1+R(j)*r1*prop)*(qRb_a(2,2,2)+qRb_b(2,2,2))
            b222=(Q(j)*q1+R(j)*r1*prop)*(qqb_a(2,2,2)+qqb_b(2,2,2))

            a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
            a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
            a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
            a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)

            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     &          +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(1,2,1)
            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     &          +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(1,2,2)
            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     &          +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(2,1,1)
            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     &          +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(2,1,2)

            mqq(0,j,k)=faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(1,j,k)=faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

       if ((j==1).or.(j==3).or.(j==5)) then
           nup=2
           ndo=nf-3
       else
           nup=1
           ndo=nf-2
       endif
       if (nflav <= 4) ndo=ndo-1
       if (nflav <= 3) nup=nup-1
            b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     &          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(1,1,1)
            b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     &          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(1,1,2)
            b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     &          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(2,2,1)
            b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     &          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(2,2,2)
            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     &          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(1,2,1)
            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     &          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(1,2,2)
            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     &          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(2,1,1)
            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     &          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(2,1,2)

      tdo=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

            b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     &          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(1,1,1)
            b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     &          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(1,1,2)
            b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     &          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(2,2,1)
            b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     &          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(2,2,2)
            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     &          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(1,2,1)
            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     &          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(1,2,2)
            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     &          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(2,1,1)
            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     &          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(2,1,2)

      tup=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

      mqq(1,j,k)=mqq(1,j,k)+real(nup,dp)*tup+real(ndo,dp)*tdo
      endif
      elseif ((j < 0) .and. (k > 0)) then
C---Qb-q case
            if (j .ne. -k) then
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

            mqq(0,j,k)=zip
            mqq(1,j,k)=zip
            mqq(2,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            elseif (j == -k) then

            a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbR_a(1,1,1)+qbR_b(1,1,1))
            b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbq_a(1,1,1)+qbq_b(1,1,1))
            a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbR_a(1,1,2)+qbR_b(1,1,2))
            b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbq_a(1,1,2)+qbq_b(1,1,2))
            a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbR_a(2,2,1)+qbR_b(2,2,1))
            b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbq_a(2,2,1)+qbq_b(2,2,1))
            a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbR_a(2,2,2)+qbR_b(2,2,2))
            b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbq_a(2,2,2)+qbq_b(2,2,2))

            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     &          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     &          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     &          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     &          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)

            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     &          +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     &          +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     &          +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     &          +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(2,1,2)

            mqq(0,j,k)=faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(2,j,k)=faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(1,j,k)=faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

c--Here we must also add the contribution of other final state quarks
c  unequal to initial annihilating quarks
       if ((k==1).or.(k==3).or.(k==5)) then
           nup=2
           ndo=nf-3
       else
           nup=1
           ndo=nf-2
       endif
       if (nflav <= 4) ndo=ndo-1
       if (nflav <= 3) nup=nup-1
            b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     &          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(1,1,1)
            b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     &          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(1,1,2)
            b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     &          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(2,2,1)
            b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     &          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(2,2,2)
            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     &          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     &          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     &          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     &          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(2,1,2)
      tdo=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

            b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     &          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(1,1,1)
            b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     &          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(1,1,2)
            b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     &          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(2,2,1)
            b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     &          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(2,2,2)
            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     &          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(1,2,1)
            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     &          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(1,2,2)
            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     &          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(2,1,1)
            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     &          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(2,1,2)
      tup=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

      mqq(1,j,k)=mqq(1,j,k)+real(nup,dp)*tup+real(ndo,dp)*tdo

          endif
          endif
      msq(j,k)=msq(j,k)+mqq(0,j,k)+mqq(1,j,k)+mqq(2,j,k)
      enddo
      enddo
      endif
      return
      end





