      subroutine qqb_w_cjet_g(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + cbar(p5) + f(p6)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5) + f(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      include 'nwz.f'
      include 'nflav.f'
      include 'ckm.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     &                 facgg,prop
      real(dp):: qgWqg2,qbgWqbg2,gqbWqbg2,gqWqg2,ggWqqb2,ggWqbq2
      real(dp):: qgWqg2_cs(0:2),qbgWqbg2_cs(0:2),
     &                 gqbWqbg2_cs(0:2),gqWqg2_cs(0:2),
     &                 ggWqbq2_cs(0:2),
     &                 ggWqqb2_cs(0:2)
      real(dp):: sbqWcbq,sqWcq,qsWcq,qsbWcbq,qqbWcsb,qqbWcbs
      real(dp):: s34
      common/facgg/facgg
c--- we label the amplitudes by helicity (qqb1 ... qqb4)
c--- and by type of contribution qqb(1) ... qqb(n)
      integer:: i,j,k
!$omp threadprivate(/facgg/)

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--- Note that the charm quark mass is calculated dynamically inside
c--- the routines wqq_sc and w2jetsq_mass

c--- calculate 4-quark contributions
c--- note that these are symmetric under interchange of q and qb
      call wqq_sc(2,5,3,4,6,1,p,qsWcq)
      call wqq_sc(1,5,3,4,6,2,p,sqWcq)

      call wqq_sc(1,5,4,3,6,2,p,sbqWcbq)
      call wqq_sc(2,5,4,3,6,1,p,qsbWcbq)

      call wqq_sc(6,5,3,4,2,1,p,qqbWcsb)
      call wqq_sc(6,5,4,3,1,2,p,qqbWcbs)

      s34=2._dp*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      prop=s34**2/((s34-wmass**2)**2+wmass**2*wwidth**2)
      facgg=V*xn/four*(gwsq/2._dp*gsq)**2*prop

c--- calculate 2-quark, 2-gluon amplitudes

      if     (nwz == +1) then
        call w2jetsq_mass(1,5,4,3,2,6,p,qbgWqbg2)
        call storecs(qbgWqbg2_cs)

        call w2jetsq_mass(2,5,4,3,1,6,p,gqbWqbg2)
        call storecs(gqbWqbg2_cs)

        call w2jetsq_mass(6,5,4,3,1,2,p,ggWqbq2)
        call storecs(ggWqbq2_cs)
        do i=0,2
          gqbWqbg2_cs(i)= aveqg*facgg*gqbWqbg2_cs(i)
          qbgWqbg2_cs(i)= aveqg*facgg*qbgWqbg2_cs(i)
          ggWqbq2_cs(i) = avegg*facgg*ggWqbq2_cs(i)
        enddo
        gqbWqbg2= gqbWqbg2_cs(1)+gqbWqbg2_cs(2)+gqbWqbg2_cs(0)
        qbgWqbg2= qbgWqbg2_cs(1)+qbgWqbg2_cs(2)+qbgWqbg2_cs(0)
        ggWqbq2 = ggWqbq2_cs(1) +ggWqbq2_cs(2) +ggWqbq2_cs(0)
      elseif (nwz == -1) then
        call w2jetsq_mass(2,5,3,4,1,6,p,gqWqg2)
        call storecs(gqWqg2_cs)

        call w2jetsq_mass(1,5,3,4,2,6,p,qgWqg2)
        call storecs(qgWqg2_cs)

        call w2jetsq_mass(6,5,3,4,1,2,p,ggWqqb2)
        call storecs(ggWqqb2_cs)
        do i=0,2
          gqWqg2_cs(i)  = aveqg*facgg*gqWqg2_cs(i)
          qgWqg2_cs(i)  = aveqg*facgg*qgWqg2_cs(i)
          ggWqqb2_cs(i) = avegg*facgg*ggWqqb2_cs(i)
        enddo
        gqWqg2  = gqWqg2_cs(1)  +gqWqg2_cs(2)  +gqWqg2_cs(0)
        qgWqg2  = qgWqg2_cs(1)  +qgWqg2_cs(2)  +qgWqg2_cs(0)
        ggWqqb2 = ggWqqb2_cs(1) +ggWqqb2_cs(2) +ggWqqb2_cs(0)
      endif


      do j=-nflav,nflav
      do k=-nflav,nflav

c--- 2-quark, 2-gluon contribution to matrix elements
      if     (((j == +1).or.(j == +3)) .and. (k == 0)) then
          msq(j,k)=Vsq(j,-4)*qgWqg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(j,-4)*qgWqg2_cs(i)
          enddo
      elseif (((j == -1).or.(j == -3)) .and. (k == 0)) then
          msq(j,k)=Vsq(j,4)*qbgWqbg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(j,4)*qbgWqbg2_cs(i)
          enddo
      elseif ((j == 0) .and. ((k == +1).or.(k == +3))) then
          msq(j,k)=Vsq(k,-4)*gqWqg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(k,-4)*gqWqg2_cs(i)
          enddo
      elseif ((j == 0) .and. ((k == -1).or.(k == -3))) then
          msq(j,k)=Vsq(k,4)*gqbWqbg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(k,4)*gqbWqbg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k == 0)) then
          if     (nwz == +1) then
            msq(j,k)=Vsum(+4)*ggWqbq2
            do i=0,2
              msq_cs(i,j,k)=Vsum(+4)*ggWqbq2_cs(i)
            enddo
          elseif (nwz == -1) then
            msq(j,k)=Vsum(-4)*ggWqqb2
            do i=0,2
              msq_cs(i,j,k)=Vsum(-4)*ggWqqb2_cs(i)
            enddo
          endif
      endif

c--- 4-quark contribution to matrix elements
      if     ((j > 0) .and. (k > 0)) then
        msq(j,k)=Vsq(j,-4)*sqWcq+Vsq(k,-4)*qsWcq
      elseif ((j < 0) .and. (k < 0)) then
        msq(j,k)=Vsq(j,+4)*sbqWcbq+Vsq(k,+4)*qsbWcbq
      elseif ((j > 0) .and. (k < 0)) then
        msq(j,k)=Vsq(j,-4)*sqWcq+Vsq(k,+4)*qsbWcbq
        if (j == -k) then
          msq(j,k)=msq(j,k)+Vsum(-4)*qqbWcsb+Vsum(+4)*qqbWcbs
        endif
      elseif ((j < 0) .and. (k > 0)) then
        msq(j,k)=Vsq(j,+4)*sbqWcbq+Vsq(k,-4)*qsWcq
        if (j == -k) then
          msq(j,k)=msq(j,k)+Vsum(-4)*qqbWcsb+Vsum(+4)*qqbWcbs
        endif
      endif

      enddo
      enddo

      return
      end

