      subroutine qqb_w2jet(p,msq)
      implicit none
      include 'types.f'

c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + qbar(-p2) --> W + f(p5) + f(p6)
c                            |
c                            --> nu(p3) + e^+(p4)
c     where the fermions are either q(p5) and qbar(p6) [Qflag = .true.]
c                                or g(p5) and g(p6)    [Gflag = .true.]
c--- all momenta are incoming
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      include 'kpart.f'
      include 'first.f'
      include 'mpicommon.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     &                 facgg,facqq,prop,Vfac
      real(dp):: qqbWgg2,qbqWgg2,qgWqg2,qbgWqbg2,
     &                 gqbWqbg2,gqWqg2,ggWqbq2
      real(dp):: qqbWgg2_cs(0:2),qbqWgg2_cs(0:2),qgWqg2_cs(0:2),
     &                 qbgWqbg2_cs(0:2),gqbWqbg2_cs(0:2),
     &                 gqWqg2_cs(0:2),ggWqbq2_cs(0:2)
      real(dp)::
     & qqb_ijkk(0:2),qqb_ijii(0:2),qqb_ijjj(0:2),qqb_ijkj(0:2),
     & qqb_ijik(0:2),qqb_ijkl(0:2),qqb_iiij(0:2),qqb_iiji(0:2),
     & qbq_ijkk(0:2),qbq_ijii(0:2),qbq_ijjj(0:2),qbq_ijkj(0:2),
     & qbq_ijik(0:2),qbq_ijkl(0:2),qbq_iiij(0:2),qbq_iiji(0:2),
     & qq_iiji(0:2),qq_ijkj(0:2),qq_ijik(0:2),
     & qq_ijjj(0:2),qq_ijii(0:2),
     & qbqb_iiji(0:2),qbqb_ijkj(0:2),qbqb_ijik(0:2),
     & qbqb_ijjj(0:2),qbqb_ijii(0:2)
      real(dp):: mqq(0:2,fn:nf,fn:nf)
      complex(dp):: qqb1(3),qqb2(3),qqb3(3),qqb4(3),
     &               qq1(4),qq2(4),qq3(4),qq4(4),
     &               qbq1(3),qbq2(3),qbq3(3),qbq4(3),
     &               qbqb1(4),qbqb2(4),qbqb3(4),qbqb4(4)
      integer:: rcolourchoice
      common/mqq/mqq
c--- we label the amplitudes by helicity (qqb1 ... qqb4)
c--- and by type of contribution qqb(1) ... qqb(n)
      integer:: i,j,k,n1,n2
!$omp threadprivate(/mqq/)

      if (first) then
        first=.false.
!        if (rank == 0) then
!        if (Gflag) then
!          write(*,*) 'Using QQGG matrix elements'
!          write(*,*) '[LC is   N ]'
!          write(*,*) '[SLC is 1/N]'
!        endif
!        if (Qflag) then
!          write(*,*) 'Using QQBQQB matrix elements'
!          write(*,*) '[LC is   1 ]'
!          write(*,*) '[SLC is 1/N]'
!        endif
!        if (kpart==klord) then
!          if     (colourchoice == 1) then
!            write(*,*) 'Leading colour only'
!          elseif (colourchoice == 2) then
!            write(*,*) 'Sub-leading colour only'
!          elseif (colourchoice == 0) then
!            write(*,*) 'Total of both colour structures'
!          else
!            write(*,*) 'Bad colourchoice'
!            stop
!          endif
!        else
!          write(*,*) 'Calculating all colour structures in LO'
!        endif
!        endif
      endif

c--- if we're calculating the REAL or VIRT matrix elements, we
c--- need all the colour structures, but want to preserve
c--- the actual value of colourchoice
      if ((kpart==kreal) .or. (kpart==kvirt)) then
        rcolourchoice=colourchoice
        colourchoice=0
      endif

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zip
      enddo
      enddo

c--- set up spinors
      call spinoru(6,p,za,zb)
      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      facqq=four*V*gsq**2*(gwsq/two)**2*aveqq*prop
      facgg=V*xn/four*(gwsq/two)**2*gsq**2*prop


c--- calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then
        call w2jetsq(1,2,3,4,5,6,za,zb,qqbWgg2)
        call storecs(qqbWgg2_cs)
        call w2jetsq(2,1,3,4,5,6,za,zb,qbqWgg2)
        call storecs(qbqWgg2_cs)
        call w2jetsq(1,5,3,4,2,6,za,zb,qgWqg2)
        call storecs(qgWqg2_cs)
        call w2jetsq(2,5,3,4,1,6,za,zb,gqWqg2)
        call storecs(gqWqg2_cs)
        call w2jetsq(5,1,3,4,2,6,za,zb,qbgWqbg2)
        call storecs(qbgWqbg2_cs)
        call w2jetsq(5,2,3,4,1,6,za,zb,gqbWqbg2)
        call storecs(gqbWqbg2_cs)
        call w2jetsq(5,6,3,4,1,2,za,zb,ggWqbq2)
        call storecs(ggWqbq2_cs)
        do i=0,2
          qqbWgg2_cs(i) = half*aveqq*facgg*qqbWgg2_cs(i)
          qbqWgg2_cs(i) = half*aveqq*facgg*qbqWgg2_cs(i)
          gqWqg2_cs(i)  = aveqg*facgg*gqWqg2_cs(i)
          qgWqg2_cs(i)  = aveqg*facgg*qgWqg2_cs(i)
          gqbWqbg2_cs(i)= aveqg*facgg*gqbWqbg2_cs(i)
          qbgWqbg2_cs(i)= aveqg*facgg*qbgWqbg2_cs(i)
          ggWqbq2_cs(i) = avegg*facgg*ggWqbq2_cs(i)
        enddo
        qqbWgg2 = qqbWgg2_cs(1) +qqbWgg2_cs(2) +qqbWgg2_cs(0)
        qbqWgg2 = qbqWgg2_cs(1) +qbqWgg2_cs(2) +qbqWgg2_cs(0)
        gqWqg2  = gqWqg2_cs(1)  +gqWqg2_cs(2)  +gqWqg2_cs(0)
        qgWqg2  = qgWqg2_cs(1)  +qgWqg2_cs(2)  +qgWqg2_cs(0)
        gqbWqbg2= gqbWqbg2_cs(1)+gqbWqbg2_cs(2)+gqbWqbg2_cs(0)
        qbgWqbg2= qbgWqbg2_cs(1)+qbgWqbg2_cs(2)+qbgWqbg2_cs(0)
        ggWqbq2 = ggWqbq2_cs(1) +ggWqbq2_cs(2) +ggWqbq2_cs(0)
      endif

c--- calculate four-quark amplitudes
      if (Qflag) then
c--- basic amplitudes - q qb --> W + g* (--> q qb) (amps 1 and 3)
c---                and q qb --> g* --> q (--> W q) qb (amps 2 and 4)
        call amp_q_QbQ_qb(1,2,5,6,qqb1(1),qqb2(1),qqb3(1),qqb4(1))
c--- crossed - q qb --> q qb ( --> W qb)
        call amp_q_QbQ_qb(1,5,2,6,qqb1(2),qqb2(2),qqb3(2),qqb4(2))
c--- crossed - q qb --> q ( --> W q) qb
        call amp_q_QbQ_qb(6,2,5,1,qqb1(3),qqb2(3),qqb3(3),qqb4(3))

c--- now the qb q amplitudes
        call amp_q_QbQ_qb(2,1,5,6,qbq1(1),qbq2(1),qbq3(1),qbq4(1))
c--- crossed - qb q --> qb q ( --> W q)
        call amp_q_QbQ_qb(2,5,1,6,qbq1(2),qbq2(2),qbq3(2),qbq4(2))
c--- crossed - qb q --> qb ( --> W qb) q
        call amp_q_QbQ_qb(6,1,5,2,qbq1(3),qbq2(3),qbq3(3),qbq4(3))

c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,5,6,2,qq1(1),qq2(1),qq3(1),qq4(1))
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,5,6,1,qq1(2),qq2(2),qq3(2),qq4(2))
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,6,5,1,qq1(3),qq2(3),qq3(3),qq4(3))
c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,6,5,2,qq1(4),qq2(4),qq3(4),qq4(4))

c--- crossed qb qb --> qb ( --> W qb) qb
        call amp_q_QbQ_qb(5,1,2,6,qbqb1(1),qbqb2(1),qbqb3(1),qbqb4(1))
c--- crossed qb qb --> qb qb ( --> W qb)
        call amp_q_QbQ_qb(5,2,1,6,qbqb1(2),qbqb2(2),qbqb3(2),qbqb4(2))
c--- crossed qb qb --> qb ( --> W qb) qb
        call amp_q_QbQ_qb(6,2,1,5,qbqb1(3),qbqb2(3),qbqb3(3),qbqb4(3))
c--- crossed qb qb --> qb qb ( --> W qb)
        call amp_q_QbQ_qb(6,1,2,5,qbqb1(4),qbqb2(4),qbqb3(4),qbqb4(4))

c--- now square these amplitudes separating into color structures
c   1) Amplitude
c   2) Amplitude with (5<-->6)
c   0) Interference between above
c
c--- q(i) qb(j) --> W + g* (--> q(k) qb(k)) with k != i,j
        qqb_ijkk(1)=abs(qqb1(1))**2+abs(qqb3(1))**2
        qqb_ijkk(2)=zip
        qqb_ijkk(0)=zip

        qbq_ijkk(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijkk(2)=zip
        qbq_ijkk(0)=zip
c--- q(i) qb(j) --> W + g* (--> q(i) qb(i)) i.e. k = i
        qqb_ijii(1)=abs(qqb1(1))**2+abs(qqb3(1))**2
        qqb_ijii(2)=abs(qqb2(2))**2+abs(qqb4(2))**2
        qqb_ijii(0)=+two/xn*real(qqb1(1)*conjg(qqb2(2)))

        qbq_ijii(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijii(2)=abs(qbq2(2))**2+abs(qbq4(2))**2
        qbq_ijii(0)=+two/xn*real(qbq1(1)*conjg(qbq2(2)))
c--- q(i) qb(j) --> W + g* (--> q(j) qb(j)) i.e. k = j
        qqb_ijjj(1)=abs(qqb1(1))**2+abs(qqb3(1))**2
        qqb_ijjj(2)=abs(qqb2(3))**2+abs(qqb4(3))**2
        qqb_ijjj(0)=two/xn*real(qqb1(1)*conjg(qqb2(3)))

        qbq_ijjj(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijjj(2)=abs(qbq2(3))**2+abs(qbq4(3))**2
        qbq_ijjj(0)=two/xn*real(qbq1(1)*conjg(qbq2(3)))
c--- q (i) qb(j) --> q(i) qb(j) ( --> W qb(k)) with k != i,j
        qqb_ijik(1)=zip
        qqb_ijik(2)=abs(qqb2(2))**2+abs(qqb4(2))**2
        qqb_ijik(0)=zip

        qbq_ijik(2)=abs(qbq2(2))**2+abs(qbq4(2))**2
        qbq_ijik(1)=zip
        qbq_ijik(0)=zip
c--- q (i) qb(j) --> q(i) ( --> W q(k)) qb(j) with k != i,j
        qqb_ijkj(1)=zip
        qqb_ijkj(2)=abs(qqb2(3))**2+abs(qqb4(3))**2
        qqb_ijkj(0)=zip

        qbq_ijkj(2)=abs(qbq2(3))**2+abs(qbq4(3))**2
        qbq_ijkj(1)=zip
        qbq_ijkj(0)=zip
c--- q(i) qb(j) --> g* --> q(l) (--> W q(k)) qb(l)
        qqb_ijkl(1)=abs(qqb2(1))**2+abs(qqb4(1))**2
        qqb_ijkl(2)=zip
        qqb_ijkl(0)=zip

        qbq_ijkl(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_ijkl(2)=zip
        qbq_ijkl(0)=zip
c--- q(i) qb(i) --> g* --> q(j) (--> W q(i)) qb(j)
        qqb_iiij(1)=abs(qqb2(1))**2+abs(qqb4(1))**2
        qqb_iiij(2)=abs(qqb1(3))**2+abs(qqb3(3))**2
        qqb_iiij(0)=two/xn*real(qqb2(1)*conjg(qqb1(3)))

        qbq_iiij(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_iiij(2)=abs(qbq1(3))**2+abs(qbq3(3))**2
        qbq_iiij(0)=two/xn*real(qbq2(1)*conjg(qbq1(3)))
c--- q(i) qb(i) --> g* --> q(j) qb(j) (--> W qb(i))
        qqb_iiji(1)=abs(qqb2(1))**2+abs(qqb4(1))**2
        qqb_iiji(2)=abs(qqb1(2))**2+abs(qqb3(2))**2
        qqb_iiji(0)=two/xn*real(qqb2(1)*conjg(qqb1(2)))

        qbq_iiji(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_iiji(2)=abs(qbq1(2))**2+abs(qbq3(2))**2
        qbq_iiji(0)=two/xn*real(qbq2(1)*conjg(qbq1(2)))

c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i)
        qq_iiji(1)=abs(qq1(1))**2+abs(qq3(1))**2
        qq_iiji(2)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_iiji(0)=two/xn*real(qq1(1)*conjg(qq1(2)))
c--- q(i) q(j) --> q(i) ( --> W q(k) ) q(j)
        qq_ijkj(1)=abs(qq1(1))**2+abs(qq3(1))**2
        qq_ijkj(2)=zip
        qq_ijkj(0)=zip
c--- q(i) q(j) --> q(i) q(j) ( --> W q(k) )
        qq_ijik(1)=abs(qq1(3))**2+abs(qq3(3))**2
        qq_ijik(2)=zip
        qq_ijik(0)=zip
c--- q(i) q(j) --> q(i) ( --> W q(j) ) q(j)
        qq_ijjj(1)=abs(qq1(1))**2+abs(qq3(1))**2
        qq_ijjj(2)=abs(qq1(4))**2+abs(qq3(4))**2
        qq_ijjj(0)=two/xn*real(qq1(1)*conjg(qq1(4)))
c--- q(i) q(j) --> q(i) q(j) ( --> W q(i) )
        qq_ijii(1)=abs(qq1(3))**2+abs(qq3(3))**2
        qq_ijii(2)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_ijii(0)=two/xn*real(qq1(3)*conjg(qq1(2)))

c--- qb(i) qb(i) --> qb(i) ( --> W qb(j) ) qb(i)
        qbqb_iiji(1)=abs(qbqb1(1))**2+abs(qbqb3(1))**2
        qbqb_iiji(2)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_iiji(0)=two/xn*real(qbqb1(1)*conjg(qbqb1(2)))
c--- qb(i) qb(j) --> qb(i) ( --> W qb(k) ) qb(j)
        qbqb_ijkj(1)=abs(qbqb1(1))**2+abs(qbqb3(1))**2
        qbqb_ijkj(2)=zip
        qbqb_ijkj(0)=zip
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(k) )
        qbqb_ijik(1)=abs(qbqb1(3))**2+abs(qbqb3(3))**2
        qbqb_ijik(2)=zip
        qbqb_ijik(0)=zip
c--- qb(i) qb(j) --> qb(i) ( --> W qb(j) ) qb(j)
        qbqb_ijjj(1)=abs(qbqb1(1))**2+abs(qbqb3(1))**2
        qbqb_ijjj(2)=abs(qbqb1(4))**2+abs(qbqb3(4))**2
        qbqb_ijjj(0)=two/xn*real(qbqb1(1)*conjg(qbqb1(4)))
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(i) )
        qbqb_ijii(2)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_ijii(1)=abs(qbqb1(3))**2+abs(qbqb3(3))**2
        qbqb_ijii(0)=two/xn*real(qbqb1(2)*conjg(qbqb1(3)))
      endif


c--- 4-quark contribution to matrix elements
      if (Qflag) then

      do j=-nf,nf
      do k=-nf,nf
      mqq(0,j,k)=zip
      mqq(1,j,k)=zip
      mqq(2,j,k)=zip
      if ((j > 0) .and. (k < 0)) then
        if (j .ne. -k) then
c--- Q QBAR - different flavours
           mqq(0,j,k)=facqq*Vsq(j,k)*(qqb_ijii(0)+qqb_ijjj(0))
           mqq(1,j,k)=facqq*Vsq(j,k)*(
     &             (nf-2)*qqb_ijkk(1)+qqb_ijii(1)+qqb_ijjj(1))
           mqq(2,j,k)=facqq*Vsq(j,k)*(qqb_ijii(2)+qqb_ijjj(2))
     &            +facqq*(Vsum(j)-Vsq(j,k))*qqb_ijkj(2)
     &            +facqq*(Vsum(k)-Vsq(j,k))*qqb_ijik(2)
        else
c--- Q QBAR - same flavours
          Vfac=zip
          do n1=1,nf
          do n2=-nf,-1
          if ((n1 .ne. j) .and. (n2 .ne. k)) then
            Vfac=Vfac+Vsq(n1,n2)
          endif
          enddo
          enddo
          mqq(0,j,k)=
     &             +facqq*Vsum(k)*qqb_iiij(0)
     &             +facqq*Vsum(j)*qqb_iiji(0)
          mqq(1,j,k)=
     &             +facqq*Vsum(k)*qqb_iiij(1)
     &             +facqq*Vsum(j)*qqb_iiji(1)
     &             +facqq*Vfac*qqb_ijkl(1)
          mqq(2,j,k)=
     &             +facqq*Vsum(k)*qqb_iiij(2)
     &             +facqq*Vsum(j)*qqb_iiji(2)
        endif
      elseif ((j < 0) .and. (k > 0)) then
        if (j .ne. -k) then
c--- QBAR Q - different flavours
          mqq(0,j,k)=facqq*Vsq(j,k)*(
     &             +qbq_ijii(0)+qbq_ijjj(0))
          mqq(1,j,k)=facqq*Vsq(j,k)*(
     &             (nf-2)*qbq_ijkk(1)
     &             +qbq_ijii(1)+qbq_ijjj(1))
          mqq(2,j,k)=facqq*Vsq(j,k)*(qbq_ijii(2)+qbq_ijjj(2))
     &             +facqq*(Vsum(k)-Vsq(j,k))*qbq_ijkj(2)
     &             +facqq*(Vsum(j)-Vsq(j,k))*qbq_ijik(2)
        else
c--- QBAR Q - same flavours
          Vfac=zip
          do n1=-nf,-1
          do n2=1,nf
          if ((n1 .ne. j) .and. (n2 .ne. k)) then
            Vfac=Vfac+Vsq(n1,n2)
          endif
          enddo
          enddo
          mqq(0,j,k)=
     &             +facqq*Vsum(j)*qbq_iiij(0)
     &             +facqq*Vsum(k)*qbq_iiji(0)
          mqq(1,j,k)=
     &             +facqq*Vsum(j)*qbq_iiij(1)
     &             +facqq*Vsum(k)*qbq_iiji(1)
     &             +facqq*Vfac*qbq_ijkl(1)
          mqq(2,j,k)=
     &             +facqq*Vsum(j)*qbq_iiij(2)
     &             +facqq*Vsum(k)*qbq_iiji(2)
        endif
      elseif ((j > 0) .and. (k > 0)) then
        if (j .ne. k) then
c--- Q Q - different flavours
          mqq(0,j,k)=
     &          +facqq*half*Vsq(j,-k)*qq_ijjj(0)
     &          +facqq*half*Vsq(k,-j)*qq_ijii(0)
          mqq(1,j,k)=
     &          +facqq*(Vsum(j)-Vsq(j,-k))*qq_ijkj(1)
     &          +facqq*(Vsum(k)-Vsq(k,-j))*qq_ijik(1)
     &          +facqq*half*Vsq(j,-k)*qq_ijjj(1)
     &          +facqq*half*Vsq(k,-j)*qq_ijii(1)
          mqq(2,j,k)=
     &          +facqq*half*Vsq(j,-k)*qq_ijjj(2)
     &          +facqq*half*Vsq(k,-j)*qq_ijii(2)
        else
c--- Q Q - same flavours
          mqq(0,j,k)=facqq*Vsum(j)*qq_iiji(0)
          mqq(1,j,k)=facqq*Vsum(j)*qq_iiji(1)
          mqq(2,j,k)=facqq*Vsum(j)*qq_iiji(2)
        endif
      elseif ((j < 0) .and. (k < 0)) then
        if (j .ne. k) then
c--- QBAR QBAR - different flavours
          mqq(0,j,k)=
     .+facqq*half*Vsq(j,-k)*qbqb_ijjj(0)
     .+facqq*half*Vsq(k,-j)*qbqb_ijii(0)
          mqq(1,j,k)=facqq*(Vsum(j)-Vsq(j,-k))*qbqb_ijkj(1)
     &              +facqq*(Vsum(k)-Vsq(k,-j))*qbqb_ijik(1)
     .+facqq*half*Vsq(j,-k)*qbqb_ijjj(1)
     .+facqq*half*Vsq(k,-j)*qbqb_ijii(1)
          mqq(2,j,k)=
     .+facqq*half*Vsq(j,-k)*qbqb_ijjj(2)
     .+facqq*half*Vsq(k,-j)*qbqb_ijii(2)
        else
c--- QBAR QBAR - same flavours
          mqq(0,j,k)=facqq*Vsum(j)*qbqb_iiji(0)
          mqq(1,j,k)=facqq*Vsum(j)*qbqb_iiji(1)
          mqq(2,j,k)=facqq*Vsum(j)*qbqb_iiji(2)
        endif
      endif
      if     (colourchoice == 1) then
        mqq(0,j,k)=zip
      elseif (colourchoice == 2) then
        mqq(1,j,k)=zip
        mqq(2,j,k)=zip
      endif
      msq(j,k)=mqq(0,j,k)+mqq(1,j,k)+mqq(2,j,k)
      enddo
      enddo

      endif

c--- 2-quark, 2-gluon contribution to matrix elements
      if (Gflag) then

      do j=-nf,nf
      do k=-nf,nf

      do i=0,2
      msq_cs(i,j,k)=zip
      enddo

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*qqbWgg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(j,k)*qqbWgg2_cs(i)
          enddo
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*qbqWgg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(j,k)*qbqWgg2_cs(i)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWqg2
          do i=0,2
            msq_cs(i,j,k)=
     &(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWqg2_cs(i)
          enddo
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqbg2
          do i=0,2
            msq_cs(i,j,k)=
     &(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqbg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWqg2
          do i=0,2
            msq_cs(i,j,k)=
     &(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWqg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqbg2
          do i=0,2
            msq_cs(i,j,k)=
     &(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqbg2_cs(i)
          enddo
      elseif ((j == 0) .and. (k == 0)) then
          Vfac=zip
          do n1=1,nf
            do n2=-nf,-1
              Vfac=Vfac+Vsq(n1,n2)
            enddo
          enddo
          msq(j,k)=msq(j,k)+Vfac*ggWqbq2
          do i=0,2
            msq_cs(i,j,k)=Vfac*ggWqbq2_cs(i)
          enddo
      endif

      enddo
      enddo

      endif

c--- restore proper colourchoice if necessary
      if ((kpart==kreal) .or. (kpart==kvirt)) then
        colourchoice=rcolourchoice
      endif

      return
      end

      subroutine amp_q_QbQ_qb(i1,i2,i5,i6,amp1,amp2,amp3,amp4)
      implicit none
      include 'types.f'

      integer:: i1,i2,i5,i6
      complex(dp):: aqqb_zbb_new,amp1,amp2,amp3,amp4
c--- Amplitudes for q(i1) + qb(i2) --> qb(i6) + q(i5) + W (-> 3+4)
c
c--- the form of this function is taken from the subroutine msqzbb
c--- in qqb_zbb.f. See there for comments.
c--- We return four amplitudes, in pairs corresponding to diagrams
c--- where the W couples to both quark lines

c--- quark i5 is left-handed
      amp1=+aqqb_zbb_new(i1,i6,i5,i2,3,4)
      amp2=-conjg(aqqb_zbb_new(i5,i2,i1,i6,4,3))

c--- quark i5 is right-handed
      amp3=-aqqb_zbb_new(i1,i5,i6,i2,3,4)
      amp4=-conjg(aqqb_zbb_new(i5,i1,i2,i6,4,3))

      return
      end

      subroutine storecs(mcs)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the W2jet matrix elements into separate arrays for each
c-- incoming parton case

      include 'mmsq_cs.f'
      integer:: i
      real(dp):: mcs(0:2)

      do i=0,2
        mcs(i)=mmsq_cs(i,+1,+1)
      enddo

      return
      end


