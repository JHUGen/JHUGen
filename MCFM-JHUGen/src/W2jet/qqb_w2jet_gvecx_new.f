      subroutine qqb_w2jet_gvecx_new(p,n,in,msq,msqv_cs,ppmsqvx)
      implicit none
      include 'types.f'

c----Matrix element for W+2jet production
C----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
C For nwz=+1
c     u(-p1)+dbar(-p2)--> g(p5)+ g(p6)+W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)--> g(p5)+ g(p6)+W^-(e^-(p3)+nbar(p4))
c---It has been checked that this gives the right matrix element
c---squared when n is replaced by two physical polarizations,
c---but this routine contains one factor of 1/2 for identical gluons
c   in the final state.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'ppwp2j.f'
C ip is the label of the emitter
      integer:: j,k,in,i,n1,n2
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,prop,n(4),Vfac
      real(dp):: p1p2(0:2,-1:1,-1:1),q1q2(0:2,-1:1,-1:1)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msqv_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: ppmsqvx(ppmax),tmpmsqvx(0:ppmax)
      common/p1p2/p1p2
      common/q1q2/q1q2
!$omp threadprivate(/p1p2/,/q1q2/)
c--- note that we will use the first index of p1p2 to label
c--- the colour structure of the squared matrix element:
c---     0 --> -ninth*(qed piece)
c---     1 --> +(qcd ordering 1)
c---     2 --> +(qcd ordering 1)
c---     3 --> Total (0+1+2)
c--- In this (x) version, the array q1q2 is supposed to provide
c--- the same information as p1p2, but with the identities of the
c--- two outgoing particles switched (eg. gg->qqb rather than gg->qbq)

      msq(:,:)=zip
      msqv_cs(:,:,:)=zip

      p1p2(:,:,:)=zip
      q1q2(:,:,:)=zip

      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      fac=V*xn/eight*(gsq*gwsq)**2*prop
      call spinoru(6,p,za,zb)
C   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)

      tmpmsqvx(:)=zip

      if (in == 1) then
C--initial-initial
        call w2jetnx(2,5,3,4,6,1,p,n,za,zb,zab,zba)
        call storecsv_px(0,+1)
        call w2jetnx(5,2,3,4,6,1,p,n,za,zb,zab,zba)
        call storecsv_px(0,-1)
c--- g g -> qb q
        call w2jetnx(5,6,3,4,2,1,p,n,za,zb,zab,zba)
        call storecsv_px(0,0)
c--- g g -> q qb
        call w2jetnx(6,5,3,4,2,1,p,n,za,zb,zab,zba)
        call storecsv_qx(0,0)
        do i=0,2
          p1p2(i,0,+1)=aveqg*fac*p1p2(i,0,+1)
          p1p2(i,0,-1)=aveqg*fac*p1p2(i,0,-1)
          p1p2(i,0, 0)=avegg*fac*p1p2(i,0,0)
          q1q2(i,0, 0)=avegg*fac*q1q2(i,0,0)
        enddo
c      p1p2(0,+1)=aveqg*fac*w2jetn(2,6,3,4,5,1,p,n,za,zb,zab,zba)
c      p1p2(0,-1)=aveqg*fac*w2jetn(6,2,3,4,5,1,p,n,za,zb,zab,zba)
      elseif (in == 2) then
        call w2jetnx(1,5,3,4,6,2,p,n,za,zb,zab,zba)
        call storecsv_px(+1,0)
        call w2jetnx(5,1,3,4,6,2,p,n,za,zb,zab,zba)
        call storecsv_px(-1,0)
c--- g g -> qb q
        call w2jetnx(5,6,3,4,1,2,p,n,za,zb,zab,zba)
        call storecsv_px(0,0)
c--- g g -> q qb
        call w2jetnx(6,5,3,4,1,2,p,n,za,zb,zab,zba)
        call storecsv_qx(0,0)
        do i=0,2
          p1p2(i,+1,0)=aveqg*fac*p1p2(i,+1,0)
          p1p2(i,-1,0)=aveqg*fac*p1p2(i,-1,0)
          p1p2(i,0, 0)=avegg*fac*p1p2(i,0,0)
          q1q2(i,0, 0)=avegg*fac*q1q2(i,0,0)
        enddo
c      p1p2(+1,0)=aveqg*fac*w2jetn(1,6,3,4,5,2,p,n,za,zb,zab,zba)
c      p1p2(-1,0)=aveqg*fac*w2jetn(6,1,3,4,5,2,p,n,za,zb,zab,zba)
      elseif (in == 5) then
        call w2jetnx(1,2,3,4,6,5,p,n,za,zb,zab,zba)
        call storecsv_px(1,-1)
        call w2jetnx(2,1,3,4,6,5,p,n,za,zb,zab,zba)
        call storecsv_px(-1,1)
        call w2jetnx(1,6,3,4,2,5,p,n,za,zb,zab,zba)
        call storecsv_px(+1,0)
        call w2jetnx(6,1,3,4,2,5,p,n,za,zb,zab,zba)
        call storecsv_px(-1,0)
        call w2jetnx(2,6,3,4,1,5,p,n,za,zb,zab,zba)
        call storecsv_px(0,+1)
        call w2jetnx(6,2,3,4,1,5,p,n,za,zb,zab,zba)
        call storecsv_px(0,-1)
        do i=0,2
          p1p2(i,1,-1)=half*aveqq*fac*p1p2(i,1,-1)
          p1p2(i,-1,1)=half*aveqq*fac*p1p2(i,-1,1)
          p1p2(i,+1,0)=aveqg*fac*p1p2(i,+1,0)
          p1p2(i,-1,0)=aveqg*fac*p1p2(i,-1,0)
          p1p2(i,0,+1)=aveqg*fac*p1p2(i,0,+1)
          p1p2(i,0,-1)=aveqg*fac*p1p2(i,0,-1)
        enddo
c      p1p2(1,-1)=aveqq*fac*w2jetn(1,2,3,4,6,5,p,n,za,zb,zab,zba)
c      p1p2(-1,1)=aveqq*fac*w2jetn(2,1,3,4,6,5,p,n,za,zb,zab,zba)
c      p1p2(+1,0)=aveqg*fac*w2jetn(1,6,3,4,2,5,p,n,za,zb,zab,zba)
c      p1p2(-1,0)=aveqg*fac*w2jetn(6,1,3,4,2,5,p,n,za,zb,zab,zba)
c      p1p2(0,+1)=aveqg*fac*w2jetn(2,6,3,4,1,5,p,n,za,zb,zab,zba)
c      p1p2(0,-1)=aveqg*fac*w2jetn(6,2,3,4,1,5,p,n,za,zb,zab,zba)
      elseif (in == 6) then
        call w2jetnx(1,2,3,4,5,6,p,n,za,zb,zab,zba)
        call storecsv_px(1,-1)
        call w2jetnx(2,1,3,4,5,6,p,n,za,zb,zab,zba)
        call storecsv_px(-1,1)
        call w2jetnx(1,5,3,4,2,6,p,n,za,zb,zab,zba)
        call storecsv_px(+1,0)
        call w2jetnx(5,1,3,4,2,6,p,n,za,zb,zab,zba)
        call storecsv_px(-1,0)
        call w2jetnx(2,5,3,4,1,6,p,n,za,zb,zab,zba)
        call storecsv_px(0,+1)
        call w2jetnx(5,2,3,4,1,6,p,n,za,zb,zab,zba)
        call storecsv_px(0,-1)
        do i=0,2
          p1p2(i,1,-1)=half*aveqq*fac*p1p2(i,1,-1)
          p1p2(i,-1,1)=half*aveqq*fac*p1p2(i,-1,1)
          p1p2(i,+1,0)=aveqg*fac*p1p2(i,+1,0)
          p1p2(i,-1,0)=aveqg*fac*p1p2(i,-1,0)
          p1p2(i,0,+1)=aveqg*fac*p1p2(i,0,+1)
          p1p2(i,0,-1)=aveqg*fac*p1p2(i,0,-1)
        enddo
c      p1p2(1,-1)=aveqq*fac*w2jetn(1,2,3,4,5,6,p,n,za,zb,zab,zba)
c      p1p2(-1,1)=aveqq*fac*w2jetn(2,1,3,4,5,6,p,n,za,zb,zab,zba)
c      p1p2(+1,0)=aveqg*fac*w2jetn(1,5,3,4,2,6,p,n,za,zb,zab,zba)
c      p1p2(-1,0)=aveqg*fac*w2jetn(5,1,3,4,2,6,p,n,za,zb,zab,zba)
c      p1p2(0,+1)=aveqg*fac*w2jetn(2,5,3,4,1,6,p,n,za,zb,zab,zba)
c      p1p2(0,-1)=aveqg*fac*w2jetn(5,2,3,4,1,6,p,n,za,zb,zab,zba)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          do i=0,2
            msqv_cs(i,j,k)=Vsq(j,k)*p1p2(i,1,-1)
          enddo
      elseif ((j < 0) .and. (k > 0)) then
          do i=0,2
            msqv_cs(i,j,k)=Vsq(j,k)*p1p2(i,-1,1)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*p1p2(i,+1,0)
          enddo
          if (j <= 4) then
            do n1=1,4
            tmpmsqvx(pp(j,k,n1,0))=Vsq(j,-n1)*(
     &                           p1p2(0,+1,0)+p1p2(1,+1,0)+p1p2(2,+1,0))
            tmpmsqvx(pp(j,k,0,n1))=tmpmsqvx(pp(j,k,n1,0))
            enddo
          endif
      elseif ((j < 0) .and. (k == 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*p1p2(i,-1,0)
          enddo
          if (j >= -4) then
            do n1=1,4
            tmpmsqvx(pp(j,k,-n1,0))=Vsq(j,n1)*(
     &                            p1p2(0,-1,0)+p1p2(1,-1,0)+p1p2(2,-1,0))
            tmpmsqvx(pp(j,k,0,-n1))=tmpmsqvx(pp(j,k,-n1,0))
            enddo
          endif
      elseif ((j == 0) .and. (k > 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*p1p2(i,0,+1)
          enddo
          if (k <= 4) then
            do n1=1,4
            tmpmsqvx(pp(j,k,n1,0))=Vsq(-n1,k)*(
     &                            p1p2(0,0,+1)+p1p2(1,0,+1)+p1p2(2,0,+1))
            tmpmsqvx(pp(j,k,0,n1))=tmpmsqvx(pp(j,k,n1,0))
            enddo
          endif
      elseif ((j == 0) .and. (k < 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*p1p2(i,0,-1)
          enddo
          if (k >= -4) then
            do n1=1,4
              tmpmsqvx(pp(j,k,-n1,0))=Vsq(n1,k)*(
     &                              p1p2(0,0,-1)+p1p2(1,0,-1)+p1p2(2,0,-1))
              tmpmsqvx(pp(j,k,0,-n1))=tmpmsqvx(pp(j,k,-n1,0))
            enddo
          endif
      elseif ((j == 0) .and. (k == 0)) then
          Vfac=0._dp
          do n1=1,4
            do n2=-4,-1
              Vfac=Vfac+Vsq(n1,n2)
              tmpmsqvx(pp(j,k,-n1,-n2))=Vsq(n1,n2)*(
     &                            p1p2(0,0,0)+p1p2(1,0,0)+p1p2(2,0,0))
              tmpmsqvx(pp(j,k,-n2,-n1))=Vsq(n1,n2)*(
     &                            q1q2(0,0,0)+q1q2(1,0,0)+q1q2(2,0,0))
            enddo
          enddo
          do i=0,2
            msqv_cs(i,j,k)=Vfac*p1p2(i,0,0)
          enddo
      endif
      msq(j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
      enddo
      enddo

c--- fill ppmsqx array to be returned
      ppmsqvx(1:ppmax)=tmpmsqvx(1:ppmax)

      return
      end





