      subroutine qqb_WH1jet_g(p,msq)
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+q(-p2)->W(p34)+H(p56)+q(p7)+q(p8)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'flags.f'
      include 'ckm.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'hbbparams.f'
      include 'lc.f'
      integer j,k,j1,j2,j3,j4,n1,n2,ig1,ig2
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),s4,s1278,prop,
     & WHqqbgg,facgg,Vfac,facqq,hdecay,
     & qbqWHgg2,qqbWHgg2,qgWHqg2,gqWHqg2,qbgWHqbg2,gqbWHqbg2,ggWHqbq2,
     & msqhtautau,msqhbb,msqhgamgam,sH
      real(dp)::
     & qqb_ijkk(0:2),qqb_ijii(0:2),qqb_ijjj(0:2),qqb_ijkj(0:2),
     & qqb_ijik(0:2),qqb_ijkl(0:2),qqb_iiij(0:2),qqb_iiji(0:2),
     & qbq_ijkk(0:2),qbq_ijii(0:2),qbq_ijjj(0:2),qbq_ijkj(0:2),
     & qbq_ijik(0:2),qbq_ijkl(0:2),qbq_iiij(0:2),qbq_iiji(0:2),
     & qq_iiji(0:2),qq_ijkj(0:2),qq_ijik(0:2),
     & qq_ijjj(0:2),qq_ijii(0:2),
     & qbqb_iiji(0:2),qbqb_ijkj(0:2),qbqb_ijik(0:2),
     & qbqb_ijjj(0:2),qbqb_ijii(0:2)
      complex(dp):: qqb1(3),qqb2(3),qqb3(3),qqb4(3),
     &               qq1(4),qq2(4),qq3(4),qq4(4),
     &               qbq1(3),qbq2(3),qbq3(3),qbq4(3),
     &               qbqb1(4),qbqb2(4),qbqb3(4),qbqb4(4)

      real(dp),parameter::stat=0.5_dp
      real(dp):: mqq(0:2,fn:nf,fn:nf)
      common/mqq/mqq
!$omp threadprivate(/mqq/)

!----begin statement functions
      s4(j1,j2,j3,j4)=
     & s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
!----end statement functions

      if (hdecaymode == 'wpwm') then
        ig1=9
        ig2=10
      else
        ig1=7
        ig2=8
      endif

      call spinoru(ig2,p,za,zb)

      s1278=s4(1,2,ig1,ig2)
      prop=     ((s1278-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(sH)
c--- adjust for fixed H->bb BR if necessary
          if (FixBrHbb) then
            hdecay=hdecay*GamHbb/GamHbb0
          endif
      elseif (hdecaymode == 'gaga') then
          sH=s(5,6)
          hdecay=msqhgamgam(sH)
      elseif (hdecaymode == 'wpwm') then
          sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
          call hwwdecay(p,5,6,7,8,hdecay)
      else
          write(6,*) 'Unimplemented decay mode in qqb_WH1jet_g'
          stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)

      facqq=aveqq*V*gsq**2*gwsq**3*wmass**2*(s(3,4)**2/prop)*hdecay
      facgg=V*xn*gsq**2*gwsq**3*wmass**2/prop*hdecay

c--- calculate four-quark amplitudes
      if (Qflag) then
c--- basic amplitudes - q qb --> W + g* (--> q qb) (amps 1 and 3)
c---                and q qb --> g* --> q (--> W q) qb (amps 2 and 4)
        call amp_q_QbQ_qb(1,2,ig1,ig2,qqb1(1),qqb2(1),qqb3(1),qqb4(1))
c--- crossed - q qb --> q qb ( --> W qb)
        call amp_q_QbQ_qb(1,ig1,2,ig2,qqb1(2),qqb2(2),qqb3(2),qqb4(2))
c--- crossed - q qb --> q ( --> W q) qb
        call amp_q_QbQ_qb(ig2,2,ig1,1,qqb1(3),qqb2(3),qqb3(3),qqb4(3))

c--- now the qb q amplitudes
        call amp_q_QbQ_qb(2,1,ig1,ig2,qbq1(1),qbq2(1),qbq3(1),qbq4(1))
c--- crossed - qb q --> qb q ( --> W q)
        call amp_q_QbQ_qb(2,ig1,1,ig2,qbq1(2),qbq2(2),qbq3(2),qbq4(2))
c--- crossed - qb q --> qb ( --> W qb) q
        call amp_q_QbQ_qb(ig2,1,ig1,2,qbq1(3),qbq2(3),qbq3(3),qbq4(3))

c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,ig1,ig2,2,qq1(1),qq2(1),qq3(1),qq4(1))
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,ig1,ig2,1,qq1(2),qq2(2),qq3(2),qq4(2))
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,ig2,ig1,1,qq1(3),qq2(3),qq3(3),qq4(3))
c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,ig2,ig1,2,qq1(4),qq2(4),qq3(4),qq4(4))

c--- crossed qb qb --> qb ( --> W qb) qb
        call amp_q_QbQ_qb(ig1,1,2,ig2,qbqb1(1),qbqb2(1),qbqb3(1),qbqb4(1))
c--- crossed qb qb --> qb qb ( --> W qb)
        call amp_q_QbQ_qb(ig1,2,1,ig2,qbqb1(2),qbqb2(2),qbqb3(2),qbqb4(2))
c--- crossed qb qb --> qb ( --> W qb) qb
        call amp_q_QbQ_qb(ig2,2,1,ig1,qbqb1(3),qbqb2(3),qbqb3(3),qbqb4(3))
c--- crossed qb qb --> qb qb ( --> W qb)
        call amp_q_QbQ_qb(ig2,1,2,ig1,qbqb1(4),qbqb2(4),qbqb3(4),qbqb4(4))

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
        qqb_ijii(0)=+two/xn*real(qqb1(1)*conjg(qqb2(2)),dp)

        qbq_ijii(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijii(2)=abs(qbq2(2))**2+abs(qbq4(2))**2
        qbq_ijii(0)=+two/xn*real(qbq1(1)*conjg(qbq2(2)),dp)
c--- q(i) qb(j) --> W + g* (--> q(j) qb(j)) i.e. k = j
        qqb_ijjj(1)=abs(qqb1(1))**2+abs(qqb3(1))**2
        qqb_ijjj(2)=abs(qqb2(3))**2+abs(qqb4(3))**2
        qqb_ijjj(0)=two/xn*real(qqb1(1)*conjg(qqb2(3)),dp)

        qbq_ijjj(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijjj(2)=abs(qbq2(3))**2+abs(qbq4(3))**2
        qbq_ijjj(0)=two/xn*real(qbq1(1)*conjg(qbq2(3)),dp)
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
        qqb_iiij(0)=two/xn*real(qqb2(1)*conjg(qqb1(3)),dp)

        qbq_iiij(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_iiij(2)=abs(qbq1(3))**2+abs(qbq3(3))**2
        qbq_iiij(0)=two/xn*real(qbq2(1)*conjg(qbq1(3)),dp)
c--- q(i) qb(i) --> g* --> q(j) qb(j) (--> W qb(i))
        qqb_iiji(1)=abs(qqb2(1))**2+abs(qqb4(1))**2
        qqb_iiji(2)=abs(qqb1(2))**2+abs(qqb3(2))**2
        qqb_iiji(0)=two/xn*real(qqb2(1)*conjg(qqb1(2)),dp)

        qbq_iiji(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_iiji(2)=abs(qbq1(2))**2+abs(qbq3(2))**2
        qbq_iiji(0)=two/xn*real(qbq2(1)*conjg(qbq1(2)),dp)

c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i)
        qq_iiji(1)=abs(qq1(1))**2+abs(qq3(1))**2
        qq_iiji(2)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_iiji(0)=two/xn*real(qq1(1)*conjg(qq1(2)),dp)
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
        qq_ijjj(0)=two/xn*real(qq1(1)*conjg(qq1(4)),dp)
c--- q(i) q(j) --> q(i) q(j) ( --> W q(i) )
        qq_ijii(1)=abs(qq1(3))**2+abs(qq3(3))**2
        qq_ijii(2)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_ijii(0)=two/xn*real(qq1(3)*conjg(qq1(2)),dp)

c--- qb(i) qb(i) --> qb(i) ( --> W qb(j) ) qb(i)
        qbqb_iiji(1)=abs(qbqb1(1))**2+abs(qbqb3(1))**2
        qbqb_iiji(2)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_iiji(0)=two/xn*real(qbqb1(1)*conjg(qbqb1(2)),dp)
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
        qbqb_ijjj(0)=two/xn*real(qbqb1(1)*conjg(qbqb1(4)),dp)
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(i) )
        qbqb_ijii(2)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_ijii(1)=abs(qbqb1(3))**2+abs(qbqb3(3))**2
        qbqb_ijii(0)=two/xn*real(qbqb1(2)*conjg(qbqb1(3)),dp)
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

      qqbWHgg2=stat*aveqq*facgg*WHqqbgg(1,2,ig1,ig2)
      qbqWHgg2=stat*aveqq*facgg*WHqqbgg(2,1,ig1,ig2)
      qgWHqg2=  aveqg*facgg*WHqqbgg(1,ig1,2,ig2)
      gqWHqg2=  aveqg*facgg*WHqqbgg(2,ig1,1,ig2)
      qbgWHqbg2=aveqg*facgg*WHqqbgg(ig1,1,2,ig2)
      gqbWHqbg2=aveqg*facgg*WHqqbgg(ig1,2,1,ig2)
      ggWHqbq2=avegg*facgg*WHqqbgg(ig1,ig2,1,2)

      do j=-nf,nf
      do k=-nf,nf
c      do i=0,2
c      msq_cs(i,j,k)=0._dp
c      enddo

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*qqbWHgg2
!          do i=0,2
!            msq_cs(i,j,k)=Vsq(j,k)*qqbWHgg2_cs(i)
!          enddo
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*qbqWHgg2
!          do i=0,2
!            msq_cs(i,j,k)=Vsq(j,k)*qbqWHgg2_cs(i)
!          enddo
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWHqg2
!          do i=0,2
!            msq_cs(i,j,k)=
!     &(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWHqg2_cs(i)
!          enddo
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWHqbg2
!          do i=0,2
!            msq_cs(i,j,k)=
!     &(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWHqbg2_cs(i)
!          enddo
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWHqg2
!          do i=0,2
!            msq_cs(i,j,k)=
!     &(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWHqg2_cs(i)
!          enddo
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWHqbg2
!          do i=0,2
!            msq_cs(i,j,k)=
!     &(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWHqbg2_cs(i)
!          enddo
      elseif ((j == 0) .and. (k == 0)) then
          Vfac=0._dp
          do n1=1,nf
            do n2=-nf,-1
              Vfac=Vfac+Vsq(n1,n2)
            enddo
          enddo
          msq(j,k)=msq(j,k)+Vfac*ggWHqbq2
!          do i=0,2
!            msq_cs(i,j,k)=Vfac*ggWHqbq2_cs(i)
!          enddo
      endif

      enddo
      enddo

      endif

      return
      end
