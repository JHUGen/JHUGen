      subroutine qqb_z2jet_gvec(p,n,in,msqv)
      implicit none
      include 'types.f'

C-----Nov 10 99 --- checked that gives the right matrix element
C     when summed over polarizations.

c----Matrix element for Z+2jet production
C----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
c     q(-p1)+qbar(-p2)--> g(p5)+ g(p6)+Z(f(p3)+af(p4))

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'msqv_cs.f'
      include 'flags.f'

C ip is the label of the emitter
C in is the label of the contracted line
      integer:: j,k,pq,pl,nquark,in,icol
      real(dp):: fac,n(4),tmp
      complex(dp):: prop,zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),ggtemp(0:2),
     & qqbZgg_cs(0:2,2,2),qbqZgg_cs(0:2,2,2),qgZqg_cs(0:2,2,2),
     & gqZqg_cs(0:2,2,2),qbgZqbg_cs(0:2,2,2),gqbZqbg_cs(0:2,2,2),
     & ggZqbq_cs(0:2,2,2),
     & qqbZgg(2,2),qbqZgg(2,2),qgZqg(2,2),gqZqg(2,2),
     & qbgZqbg(2,2),gqbZqbg(2,2),ggZqbq(2,2)
      include 'cplx.h'


      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      do icol=0,2
        msqv_cs(icol,j,k)=0._dp
      enddo
      enddo
      enddo

      call spinoru(6,p,za,zb)
C   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector
c---Conventions of Bern, Dixon, Kosower, Weinzierl,
c---ie, za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)

      do icol=0,2
      do pq=1,2
      do pl=1,2
      qqbZgg_cs(icol,pq,pl) =0._dp
      gqZqg_cs(icol,pq,pl)  =0._dp
      qgZqg_cs(icol,pq,pl)  =0._dp
      qbqZgg_cs(icol,pq,pl) =0._dp
      gqbZqbg_cs(icol,pq,pl)=0._dp
      qbgZqbg_cs(icol,pq,pl)=0._dp
      ggZqbq_cs(icol,pq,pl) =0._dp
      enddo
      enddo
      enddo

      if (in == 1) then
Cargument 1-4 represent (1) incoming quark line
C                       (2) incoming anti-quark line
C                       (3) outgoing gluon line
C                       (4) outgoing gluon line contracted with n
           call z2jetsqn(5,6,2,1,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv(ggZqbq_cs)
           call z2jetsqn(2,5,6,1,p,n,za,zb,zab,zba,gqZqg)
           call storezcsv(gqZqg_cs)
           call z2jetsqn(5,2,6,1,p,n,za,zb,zab,zba,gqbZqbg)
           call storezcsv(gqbZqbg_cs)
      elseif (in == 2)  then
           call z2jetsqn(5,6,1,2,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv(ggZqbq_cs)
c--- since we have interchanged 1 and 2 to get the gg matrix element,
c---  the colour structures should be interchanged too
           do pq=1,2
           do pl=1,2
           tmp=ggZqbq_cs(1,pq,pl)
           ggZqbq_cs(1,pq,pl)=ggZqbq_cs(2,pq,pl)
           ggZqbq_cs(2,pq,pl)=tmp
           enddo
           enddo
           call z2jetsqn(5,1,6,2,p,n,za,zb,zab,zba,qbgZqbg)
           call storezcsv(qbgZqbg_cs)
           call z2jetsqn(1,5,6,2,p,n,za,zb,zab,zba,qgZqg)
           call storezcsv(qgZqg_cs)
      elseif (in == 5) then
          call z2jetsqn(1,2,6,5,p,n,za,zb,zab,zba,qqbZgg)
          call storezcsv(qqbZgg_cs)
          call z2jetsqn(2,1,6,5,p,n,za,zb,zab,zba,qbqZgg)
          call storezcsv(qbqZgg_cs)
          call z2jetsqn(1,6,2,5,p,n,za,zb,zab,zba,qgZqg)
          call storezcsv(qgZqg_cs)
          call z2jetsqn(2,6,1,5,p,n,za,zb,zab,zba,gqZqg)
          call storezcsv(gqZqg_cs)
          call z2jetsqn(6,1,2,5,p,n,za,zb,zab,zba,qbgZqbg)
          call storezcsv(qbgZqbg_cs)
          call z2jetsqn(6,2,1,5,p,n,za,zb,zab,zba,gqbZqbg)
          call storezcsv(gqbZqbg_cs)
      elseif (in == 6) then
          call z2jetsqn(1,2,5,6,p,n,za,zb,zab,zba,qqbZgg)
          call storezcsv(qqbZgg_cs)
          call z2jetsqn(2,1,5,6,p,n,za,zb,zab,zba,qbqZgg)
          call storezcsv(qbqZgg_cs)
          call z2jetsqn(1,5,2,6,p,n,za,zb,zab,zba,qgZqg)
          call storezcsv(qgZqg_cs)
          call z2jetsqn(2,5,1,6,p,n,za,zb,zab,zba,gqZqg)
          call storezcsv(gqZqg_cs)
          call z2jetsqn(5,1,2,6,p,n,za,zb,zab,zba,qbgZqbg)
          call storezcsv(qbgZqbg_cs)
          call z2jetsqn(5,2,1,6,p,n,za,zb,zab,zba,gqbZqbg)
          call storezcsv(gqbZqbg_cs)
      endif

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      fac=v*xn/four*(esq*gsq)**2*two

      do icol=0,2
      do pq=1,2
      do pl=1,2
      qqbZgg_cs(icol,pq,pl) =half*aveqq*fac*qqbZgg_cs(icol,pq,pl)
      qbqZgg_cs(icol,pq,pl) =half*aveqq*fac*qbqZgg_cs(icol,pq,pl)
      gqZqg_cs(icol,pq,pl)  =aveqg*fac*gqZqg_cs(icol,pq,pl)
      qgZqg_cs(icol,pq,pl)  =aveqg*fac*qgZqg_cs(icol,pq,pl)
      gqbZqbg_cs(icol,pq,pl)=aveqg*fac*gqbZqbg_cs(icol,pq,pl)
      qbgZqbg_cs(icol,pq,pl)=aveqg*fac*qbgZqbg_cs(icol,pq,pl)
      ggZqbq_cs(icol,pq,pl) =avegg*fac*ggZqbq_cs(icol,pq,pl)
      enddo
      enddo
      enddo

      if (Gflag) then
      do j=-nf,nf
      do k=-nf,nf
      if( (j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 19

      if     ((j == 0) .and. (k == 0)) then

           do icol=0,2
           ggtemp(icol)=0._dp
             do nquark=1,nf
             ggtemp(icol)=ggtemp(icol)
     &       +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggZqbq_cs(icol,1,1)
     &       +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggZqbq_cs(icol,2,2)
     &       +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggZqbq_cs(icol,1,2)
     &       +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggZqbq_cs(icol,2,1)
             enddo
           msqv_cs(icol,j,k)=ggtemp(icol)
           enddo
      elseif ((j > 0) .and. (k < 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgg_cs(icol,1,1)
     &       +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgg_cs(icol,2,2)
     &       +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgg_cs(icol,1,2)
     &       +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgg_cs(icol,2,1)
          enddo
      elseif ((j < 0) .and. (k > 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgg_cs(icol,1,1)
     &       +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgg_cs(icol,2,2)
     &       +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgg_cs(icol,1,2)
     &       +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgg_cs(icol,2,1)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg_cs(icol,1,1)
     &       +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg_cs(icol,2,2)
     &       +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg_cs(icol,1,2)
     &       +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg_cs(icol,2,1)
          enddo
      elseif ((j < 0) .and. (k == 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg_cs(icol,1,1)
     &       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg_cs(icol,2,2)
     &       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg_cs(icol,1,2)
     &       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg_cs(icol,2,1)
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg_cs(icol,1,1)
     &       +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg_cs(icol,2,2)
     &       +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg_cs(icol,1,2)
     &       +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg_cs(icol,2,1)
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg_cs(icol,1,1)
     &       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg_cs(icol,2,2)
     &       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg_cs(icol,1,2)
     &       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg_cs(icol,2,1)
          enddo
      endif
      msqv(j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)

   19 continue
      enddo
      enddo
      endif

      return
      end


      subroutine storezcsv(p1p2)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the Z2jet_gvec matrix elements into elements of p1p2

      include 'mmsqv_cs.f'
      integer:: icol,pq,pl
      real(dp):: p1p2(0:2,2,2)

      do pq=1,2
      do pl=1,2
      do icol=0,2
        p1p2(icol,pq,pl)=mmsqv_cs(icol,pq,pl)
      enddo
      enddo
      enddo

      return
      end

