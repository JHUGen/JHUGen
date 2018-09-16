      subroutine qqb_w2jetx(p,msq,mqq,msqx,msqx_cs)
      implicit none
      include 'types.f'

c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + qbar(-p2) --> W + f(p5) + f(p6)
c                            |
c                            --> nu(p3) + e^+(p4)
c     where the fermions are either q(p5) and qbar(p6) [Qflag = .true.]
c                                or g(p5) and g(p6)    [Gflag = .true.]
c--- all momenta are incoming
c--- This routine also calculates msq for each specific final state
c--- for the 4-quark piece, ie. qi + qj --> qk + ql
c--- and returns the result in msqx(col,i,j,k,l)
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
      include 'lc.f'
      include 'kpart.f'
      include 'first.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     &                 facgg,facqq,prop,Vfac
      real(dp):: qqbWgg2,qbqWgg2,qgWqg2,qbgWqbg2,
     &                 gqbWqbg2,gqWqg2,ggWqbq2,ggWqqb2,
     &                 gqWgq2,gqbWgqb2,qgWgq2,qbgWgqb2
      real(dp):: qqbWgg2_cs(0:2),qbqWgg2_cs(0:2),qgWqg2_cs(0:2),
     &                 qbgWqbg2_cs(0:2),gqbWqbg2_cs(0:2),
     &                 gqWqg2_cs(0:2),ggWqbq2_cs(0:2),ggWqqb2_cs(0:2),
     &                 gqWgq2_cs(0:2),gqbWgqb2_cs(0:2),
     &                 qgWgq2_cs(0:2),qbgWgqb2_cs(0:2)
      real(dp)::
     & qqb_ijkk(0:2),qqb_ijii(0:2),qqb_ijjj(0:2),qqb_ijkj(0:2),
     & qqb_ijik(0:2),qqb_ijkl(0:2),qqb_iiij(0:2),qqb_iiji(0:2),
     & qbq_ijkk(0:2),qbq_ijii(0:2),qbq_ijjj(0:2),qbq_ijkj(0:2),
     & qbq_ijik(0:2),qbq_ijkl(0:2),qbq_iiij(0:2),qbq_iiji(0:2),
     & qq_iiji(0:2),qq_ijkj(0:2),qq_ijik(0:2),
     & qq_ijjj(0:2),qq_ijii(0:2),
     & qbqb_iiji(0:2),qbqb_ijkj(0:2),qbqb_ijik(0:2),
     & qbqb_ijjj(0:2),qbqb_ijii(0:2),
     & qq_iiij(0:2),qq_ijjk(0:2),qq_ijki(0:2),
     & qbqb_iiij(0:2),qbqb_ijjk(0:2),qbqb_ijki(0:2),
     & qqbXijkl(0:2),qbqXijkl(0:2),qqbXiiij(0:2),qbqXiiij(0:2),
     & qqbXiiji(0:2),qbqXiiji(0:2),qqbXijii(0:2),qbqXijii(0:2),
     & qqbXijjj(0:2),qbqXijjj(0:2),qqbxijik(0:2),qbqxijik(0:2),
     & qqbxijkj(0:2),qbqxijkj(0:2)
      complex(dp):: qqb1(6),qqb2(6),qqb3(6),qqb4(6),
     &               qq1(4),qq2(4),qq3(4),qq4(4),
     &               qbq1(6),qbq2(6),qbq3(6),qbq4(6),
     &               qbqb1(4),qbqb2(4),qbqb3(4),qbqb4(4)
      integer:: rcolourchoice
      real(dp):: mqq(0:2,fn:nf,fn:nf)
      real(dp):: msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)
c--- we label the amplitudes by helicity (qqb1 ... qqb4)
c--- and by type of contribution qqb(1) ... qqb(n)
      integer:: i,j,k,l,m,n1,n2
      integer, parameter :: sw(0:2)=(/0,2,1/)
      logical:: rGflag

      if (first) then
      first=.false.
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG matrix elements'
          write(*,*) '[LC is   N ]'
          write(*,*) '[SLC is 1/N]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB matrix elements'
          write(*,*) '[LC is   1 ]'
          write(*,*) '[SLC is 1/N]'
        endif
c        if (kpart==klord) then
          if     (colourchoice == 1) then
            write(*,*) 'Leading colour only'
          elseif (colourchoice == 2) then
            write(*,*) 'Sub-leading colour only'
          elseif (colourchoice == 0) then
            write(*,*) 'Total of both colour structures'
          else
            write(*,*) 'Bad colourchoice'
            stop
          endif
c        else
c          write(*,*) 'Calculating all colour structures in LO'
c        endif
      endif

c--- if we're calculating the REAL or VIRT matrix elements, we
c--- need all the colour structures, but want to preserve
c--- the actual value of colourchoice
      if ((kpart==kreal) .or. (kpart==kvirt)) then
        rcolourchoice=colourchoice
c--- 21/5/09: we only need all the information for Gflag, not Qflag
        if (Gflag) then
          colourchoice=0
        endif
      endif
c--- if we're calculating the REAL matrix elements with Qflag=TRUE,
c    the subtraction terms involve the (Gflag=TRUE) matrix elements
      if ((kpart==kreal) .and. (Qflag .eqv. .true.)) then
        rGflag=Gflag
        Gflag=.true.
      endif

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo
      mqq=0._dp
      msqx=0._dp
      msqx_cs=0._dp

c--- set up spinors
      call spinoru(6,p,za,zb)
      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      facqq=4._dp*V*gsq**2*(gwsq/2._dp)**2*aveqq*prop
      facgg=V*xn/four*(gwsq/2._dp)**2*gsq**2*prop


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
        call w2jetsq(6,5,3,4,1,2,za,zb,ggWqqb2)
        call storecs(ggWqqb2_cs)
        call w2jetsq(2,6,3,4,1,5,za,zb,gqWgq2)
        call storecs(gqWgq2_cs)
        call w2jetsq(6,2,3,4,1,5,za,zb,gqbWgqb2)
        call storecs(gqbWgqb2_cs)
        call w2jetsq(1,6,3,4,2,5,za,zb,qgWgq2)
        call storecs(qgWgq2_cs)
        call w2jetsq(6,1,3,4,2,5,za,zb,qbgWgqb2)
        call storecs(qbgWgqb2_cs)
        do i=0,2
          qqbWgg2_cs(i) = half*aveqq*facgg*qqbWgg2_cs(i)
          qbqWgg2_cs(i) = half*aveqq*facgg*qbqWgg2_cs(i)
          gqWqg2_cs(i)  = aveqg*facgg*gqWqg2_cs(i)
          qgWqg2_cs(i)  = aveqg*facgg*qgWqg2_cs(i)
          gqbWqbg2_cs(i)= aveqg*facgg*gqbWqbg2_cs(i)
          qbgWqbg2_cs(i)= aveqg*facgg*qbgWqbg2_cs(i)
          ggWqbq2_cs(i) = avegg*facgg*ggWqbq2_cs(i)
          ggWqqb2_cs(i) = avegg*facgg*ggWqqb2_cs(i)
          gqWgq2_cs(i)  = aveqg*facgg*gqWgq2_cs(i)
          gqbWgqb2_cs(i)= aveqg*facgg*gqbWgqb2_cs(i)
          qgWgq2_cs(i)  = aveqg*facgg*qgWgq2_cs(i)
          qbgWgqb2_cs(i)= aveqg*facgg*qbgWgqb2_cs(i)
        enddo
        qqbWgg2 = qqbWgg2_cs(1) +qqbWgg2_cs(2) +qqbWgg2_cs(0)
        qbqWgg2 = qbqWgg2_cs(1) +qbqWgg2_cs(2) +qbqWgg2_cs(0)
        gqWqg2  = gqWqg2_cs(1)  +gqWqg2_cs(2)  +gqWqg2_cs(0)
        qgWqg2  = qgWqg2_cs(1)  +qgWqg2_cs(2)  +qgWqg2_cs(0)
        gqbWqbg2= gqbWqbg2_cs(1)+gqbWqbg2_cs(2)+gqbWqbg2_cs(0)
        qbgWqbg2= qbgWqbg2_cs(1)+qbgWqbg2_cs(2)+qbgWqbg2_cs(0)
        ggWqbq2 = ggWqbq2_cs(1) +ggWqbq2_cs(2) +ggWqbq2_cs(0)
        ggWqqb2 = ggWqqb2_cs(1) +ggWqqb2_cs(2) +ggWqqb2_cs(0)
        gqWgq2  = gqWgq2_cs(1)  +gqWgq2_cs(2)  +gqWgq2_cs(0)
        gqbWgqb2= gqbWgqb2_cs(1)+gqbWgqb2_cs(2)+gqbWgqb2_cs(0)
        qgWgq2  = qgWgq2_cs(1)  +qgWgq2_cs(2)  +qgWgq2_cs(0)
        qbgWgqb2= qbgWgqb2_cs(1)+qbgWgqb2_cs(2)+qbgWgqb2_cs(0)
      endif

c--- calculate four-quark amplitudes
c--- note: amplitudes that have been commented out have been
c---       replaced by implementing symmetries
      if (Qflag) then
c--- basic amplitudes - q qb --> W + g* (--> q qb) (amps 1 and 3)
c---                and q qb --> g* --> q (--> W q) qb (amps 2 and 4)
        call amp_q_QbQ_qb(1,2,5,6,qqb1(1),qqb2(1),qqb3(1),qqb4(1))
c--- crossed - q qb --> q qb ( --> W qb)
c        call amp_q_QbQ_qb(1,5,2,6,qqb1(2),qqb2(2),qqb3(2),qqb4(2))
c--- crossed - q qb --> q ( --> W q) qb
c        call amp_q_QbQ_qb(6,2,5,1,qqb1(3),qqb2(3),qqb3(3),qqb4(3))
c--- NEW: ADDED 7/17/01
        call amp_q_QbQ_qb(1,2,6,5,qqb1(4),qqb2(4),qqb3(4),qqb4(4))
c        call amp_q_QbQ_qb(1,6,2,5,qqb1(5),qqb2(5),qqb3(5),qqb4(5))
c        call amp_q_QbQ_qb(5,2,6,1,qqb1(6),qqb2(6),qqb3(6),qqb4(6))

c--- now the qb q amplitudes
        call amp_q_QbQ_qb(2,1,5,6,qbq1(1),qbq2(1),qbq3(1),qbq4(1))
c--- crossed - qb q --> qb q ( --> W q)
        call amp_q_QbQ_qb(2,5,1,6,qbq1(2),qbq2(2),qbq3(2),qbq4(2))
c--- crossed - qb q --> qb ( --> W qb) q
c        call amp_q_QbQ_qb(6,1,5,2,qbq1(3),qbq2(3),qbq3(3),qbq4(3))
c--- NEW: ADDED 7/17/01
c        call amp_q_QbQ_qb(2,1,6,5,qbq1(4),qbq2(4),qbq3(4),qbq4(4))
c        call amp_q_QbQ_qb(2,6,1,5,qbq1(5),qbq2(5),qbq3(5),qbq4(5))
c        call amp_q_QbQ_qb(5,1,6,2,qbq1(6),qbq2(6),qbq3(6),qbq4(6))

c--- crossed q q --> q ( --> W q) q
        call amp_q_QbQ_qb(1,5,6,2,qq1(1),qq2(1),qq3(1),qq4(1))
c--- crossed q q --> q q ( --> W q)
        call amp_q_QbQ_qb(2,5,6,1,qq1(2),qq2(2),qq3(2),qq4(2))
c--- crossed q q --> q q ( --> W q)
c        call amp_q_QbQ_qb(2,6,5,1,qq1(3),qq2(3),qq3(3),qq4(3))
c--- crossed q q --> q ( --> W q) q
c        call amp_q_QbQ_qb(1,6,5,2,qq1(4),qq2(4),qq3(4),qq4(4))

c--- crossed qb qb --> qb ( --> W qb) qb
        call amp_q_QbQ_qb(5,1,2,6,qbqb1(1),qbqb2(1),qbqb3(1),qbqb4(1))
c--- crossed qb qb --> qb qb ( --> W qb)
        call amp_q_QbQ_qb(5,2,1,6,qbqb1(2),qbqb2(2),qbqb3(2),qbqb4(2))
c--- crossed qb qb --> qb ( --> W qb) qb
c        call amp_q_QbQ_qb(6,2,1,5,qbqb1(3),qbqb2(3),qbqb3(3),qbqb4(3))
c--- crossed qb qb --> qb qb ( --> W qb)
c        call amp_q_QbQ_qb(6,1,2,5,qbqb1(4),qbqb2(4),qbqb3(4),qbqb4(4))

c--- implement symmetries

        qqb1(2)=-qq3(1)
        qqb2(2)=+qbqb4(1)
        qqb3(2)=-qq1(1)
        qqb4(2)=+qbqb2(1)

        qqb1(3)=+qbqb4(1)
        qqb2(3)=-qq3(1)
        qqb3(3)=-qbqb2(1)
        qqb4(3)=+qq1(1)

        qqb1(5)=+qq4(2)
        qqb2(5)=-qbqb3(2)
        qqb3(5)=-qq2(2)
        qqb4(5)=+qbqb1(2)

        qqb1(6)=-qbqb3(2)
        qqb2(6)=+qq4(2)
        qqb3(6)=-qbqb1(2)
        qqb4(6)=+qq2(2)

        qbq1(3)=+qbq2(2)
        qbq2(3)=+qbq1(2)
        qbq3(3)=-qbqb2(2)
        qbq4(3)=+qq1(2)

        qbq1(4)=-qbq3(1)
        qbq2(4)=+qqb4(4)
        qbq3(4)=-qbq1(1)
        qbq4(4)=+qqb2(4)

        qbq1(5)=+qq4(1)
        qbq2(5)=-qbqb3(1)
        qbq3(5)=-qq2(1)
        qbq4(5)=+qbqb1(1)

        qbq1(6)=-qbqb3(1)
        qbq2(6)=+qq4(1)
        qbq3(6)=-qbqb1(1)
        qbq4(6)=+qq2(1)

        qq1(3)=+qq2(1)
        qq2(3)=+qq1(1)
        qq3(3)=-qq4(1)
        qq4(3)=-qq3(1)

        qq1(4)=+qq2(2)
        qq2(4)=+qq1(2)
        qq3(4)=-qq4(2)
        qq4(4)=-qq3(2)

        qbqb1(3)=+qbqb2(1)
        qbqb2(3)=+qbqb1(1)
        qbqb3(3)=-qbqb4(1)
        qbqb4(3)=-qbqb3(1)

        qbqb1(4)=+qbqb2(2)
        qbqb2(4)=+qbqb1(2)
        qbqb3(4)=-qbqb4(2)
        qbqb4(4)=-qbqb3(2)

c--- code for detecting symmetries
c        do i=1,4
c        write(*,*) 'i=',i
c        write(*,81) ' qqb',qqb1(i),qqb2(i),qqb3(i),qqb4(i)
c        write(*,81) ' qbq',qbq1(i),qbq2(i),qbq3(i),qbq4(i)
c        write(*,81) '  qq',qq1(i),qq2(i),qq3(i),qq4(i)
c        write(*,81) 'qbqb',qbqb1(i),qbqb2(i),qbqb3(i),qbqb4(i)
c        enddo
c        write(*,*) 'i= 5'
c        write(*,81) ' qqb',qqb1(5),qqb2(5),qqb3(5),qqb4(5)
c        write(*,81) ' qbq',qbq1(5),qbq2(5),qbq3(5),qbq4(5)
c        write(*,*) 'i= 6'
c        write(*,81) ' qqb',qqb1(6),qqb2(6),qqb3(6),qqb4(6)
c        write(*,81) ' qbq',qbq1(6),qbq2(6),qbq3(6),qbq4(6)
c        pause
c   81   format(a5,4('(',e12.5,',',e13.5,') '))

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
        qqb_ijii(0)=+2._dp/xn*real(qqb1(1)*conjg(qqb2(2)))

        qbq_ijii(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijii(2)=abs(qbq2(2))**2+abs(qbq4(2))**2
        qbq_ijii(0)=+2._dp/xn*real(qbq1(1)*conjg(qbq2(2)))
c--- q(i) qb(j) --> W + g* (--> q(j) qb(j)) i.e. k = j
        qqb_ijjj(1)=abs(qqb1(1))**2+abs(qqb3(1))**2
        qqb_ijjj(2)=abs(qqb2(3))**2+abs(qqb4(3))**2
        qqb_ijjj(0)=2._dp/xn*real(qqb1(1)*conjg(qqb2(3)))

        qbq_ijjj(1)=abs(qbq1(1))**2+abs(qbq3(1))**2
        qbq_ijjj(2)=abs(qbq2(3))**2+abs(qbq4(3))**2
        qbq_ijjj(0)=2._dp/xn*real(qbq1(1)*conjg(qbq2(3)))
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
        qqb_iiij(0)=2._dp/xn*real(qqb2(1)*conjg(qqb1(3)))

        qbq_iiij(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_iiij(2)=abs(qbq1(3))**2+abs(qbq3(3))**2
        qbq_iiij(0)=2._dp/xn*real(qbq2(1)*conjg(qbq1(3)))
c--- q(i) qb(i) --> g* --> q(j) qb(j) (--> W qb(i))
        qqb_iiji(1)=abs(qqb2(1))**2+abs(qqb4(1))**2
        qqb_iiji(2)=abs(qqb1(2))**2+abs(qqb3(2))**2
        qqb_iiji(0)=2._dp/xn*real(qqb2(1)*conjg(qqb1(2)))

        qbq_iiji(1)=abs(qbq2(1))**2+abs(qbq4(1))**2
        qbq_iiji(2)=abs(qbq1(2))**2+abs(qbq3(2))**2
        qbq_iiji(0)=2._dp/xn*real(qbq2(1)*conjg(qbq1(2)))

c--- NEW: ADDED 7/17/01
c--- q(i) qb(j) --> g* --> q(l) (--> W q(k)) qb(l)
        qqbXijkl(1)=abs(qqb2(4))**2+abs(qqb4(4))**2
        qqbXijkl(2)=zip
        qqbXijkl(0)=zip

        qbqXijkl(1)=abs(qbq2(4))**2+abs(qbq4(4))**2
        qbqXijkl(2)=zip
        qbqXijkl(0)=zip
c--- q(i) qb(i) --> g* --> q(j) (--> W q(i)) qb(j)
        qqbXiiij(1)=abs(qqb2(4))**2+abs(qqb4(4))**2
        qqbXiiij(2)=abs(qqb1(6))**2+abs(qqb3(6))**2
        qqbXiiij(0)=2._dp/xn*real(qqb2(4)*conjg(qqb1(6)))

        qbqXiiij(1)=abs(qbq2(4))**2+abs(qbq4(4))**2
        qbqXiiij(2)=abs(qbq1(6))**2+abs(qbq3(6))**2
        qbqXiiij(0)=2._dp/xn*real(qbq2(4)*conjg(qbq1(6)))
c--- q(i) qb(i) --> g* --> q(j) qb(j) (--> W qb(i))
        qqbXiiji(1)=abs(qqb2(4))**2+abs(qqb4(4))**2
        qqbXiiji(2)=abs(qqb1(5))**2+abs(qqb3(5))**2
        qqbXiiji(0)=2._dp/xn*real(qqb2(4)*conjg(qqb1(5)))

        qbqXiiji(1)=abs(qbq2(4))**2+abs(qbq4(4))**2
        qbqXiiji(2)=abs(qbq1(5))**2+abs(qbq3(5))**2
        qbqXiiji(0)=2._dp/xn*real(qbq2(4)*conjg(qbq1(5)))
c--- q(i) qb(j) --> W + g* (--> q(i) qb(i)) i.e. k = i
        qqbxijii(1)=abs(qqb1(4))**2+abs(qqb3(4))**2
        qqbxijii(2)=abs(qqb2(5))**2+abs(qqb4(5))**2
        qqbxijii(0)=+2._dp/xn*real(qqb1(4)*conjg(qqb2(5)))

        qbqxijii(1)=abs(qbq1(4))**2+abs(qbq3(4))**2
        qbqxijii(2)=abs(qbq2(5))**2+abs(qbq4(5))**2
        qbqxijii(0)=+2._dp/xn*real(qbq1(4)*conjg(qbq2(5)))
c--- q(i) qb(j) --> W + g* (--> q(j) qb(j)) i.e. k = j
        qqbxijjj(1)=abs(qqb1(4))**2+abs(qqb3(4))**2
        qqbxijjj(2)=abs(qqb2(6))**2+abs(qqb4(6))**2
        qqbxijjj(0)=2._dp/xn*real(qqb1(4)*conjg(qqb2(6)))

        qbqxijjj(1)=abs(qbq1(4))**2+abs(qbq3(4))**2
        qbqxijjj(2)=abs(qbq2(6))**2+abs(qbq4(6))**2
        qbqxijjj(0)=2._dp/xn*real(qbq1(4)*conjg(qbq2(6)))
c--- q (i) qb(j) --> q(i) qb(j) ( --> W qb(k)) with k != i,j
        qqbxijik(1)=zip
        qqbxijik(2)=abs(qqb2(5))**2+abs(qqb4(5))**2
        qqbxijik(0)=zip

        qbqxijik(2)=abs(qbq2(5))**2+abs(qbq4(5))**2
        qbqxijik(1)=zip
        qbqxijik(0)=zip
c--- q (i) qb(j) --> q(i) ( --> W q(k)) qb(j) with k != i,j
        qqbxijkj(1)=zip
        qqbxijkj(2)=abs(qqb2(6))**2+abs(qqb4(6))**2
        qqbxijkj(0)=zip

        qbqxijkj(2)=abs(qbq2(6))**2+abs(qbq4(6))**2
        qbqxijkj(1)=zip
        qbqxijkj(0)=zip

c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i)
        qq_iiji(1)=abs(qq1(1))**2+abs(qq3(1))**2
        qq_iiji(2)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_iiji(0)=2._dp/xn*real(qq1(1)*conjg(qq1(2)))
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
        qq_ijjj(0)=2._dp/xn*real(qq1(1)*conjg(qq1(4)))
c--- q(i) q(j) --> q(i) q(j) ( --> W q(i) )
        qq_ijii(1)=abs(qq1(3))**2+abs(qq3(3))**2
        qq_ijii(2)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_ijii(0)=2._dp/xn*real(qq1(3)*conjg(qq1(2)))
c--- NEW: ADDED 7/17/01
c--- q(i) q(i) --> q(i) ( --> W q(j) ) q(i)
        qq_iiij(1)=abs(qq1(4))**2+abs(qq3(4))**2
        qq_iiij(2)=abs(qq1(3))**2+abs(qq3(3))**2
        qq_iiij(0)=2._dp/xn*real(qq1(4)*conjg(qq1(3)))
c--- q(i) q(j) --> q(i) q(j) ( --> W q(k) )
        qq_ijki(1)=abs(qq1(2))**2+abs(qq3(2))**2
        qq_ijki(2)=zip
        qq_ijki(0)=zip
c--- q(i) q(j) --> q(i) ( --> W q(k) ) q(j)
        qq_ijjk(1)=abs(qq1(4))**2+abs(qq3(4))**2
        qq_ijjk(2)=zip
        qq_ijjk(0)=zip

c--- qb(i) qb(i) --> qb(i) ( --> W qb(j) ) qb(i)
        qbqb_iiji(1)=abs(qbqb1(1))**2+abs(qbqb3(1))**2
        qbqb_iiji(2)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_iiji(0)=2._dp/xn*real(qbqb1(1)*conjg(qbqb1(2)))
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
        qbqb_ijjj(0)=2._dp/xn*real(qbqb1(1)*conjg(qbqb1(4)))
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(i) )
        qbqb_ijii(2)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_ijii(1)=abs(qbqb1(3))**2+abs(qbqb3(3))**2
        qbqb_ijii(0)=2._dp/xn*real(qbqb1(2)*conjg(qbqb1(3)))
c--- NEW: ADDED 7/17/01
c--- qb(i) qb(i) --> qb(i) ( --> W qb(j) ) qb(i)
        qbqb_iiij(1)=abs(qbqb1(4))**2+abs(qbqb3(4))**2
        qbqb_iiij(2)=abs(qbqb1(3))**2+abs(qbqb3(3))**2
        qbqb_iiij(0)=2._dp/xn*real(qbqb1(4)*conjg(qbqb1(3)))
c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(k) )
        qbqb_ijki(1)=abs(qbqb1(2))**2+abs(qbqb3(2))**2
        qbqb_ijki(2)=zip
        qbqb_ijki(0)=zip
c--- qb(i) qb(j) --> qb(i) ( --> W qb(k) ) qb(j)
        qbqb_ijjk(1)=abs(qbqb1(4))**2+abs(qbqb3(4))**2
        qbqb_ijjk(2)=zip
        qbqb_ijjk(0)=zip
      endif

c--- remove all subleading colour pieces if they are not required
      if (colourchoice == 1) then
        qqb_ijkk(0)=zip
        qqb_ijii(0)=zip
        qqb_ijjj(0)=zip
        qqb_ijkj(0)=zip
        qqb_ijik(0)=zip
        qqb_ijkl(0)=zip
        qqb_iiij(0)=zip
        qqb_iiji(0)=zip
        qbq_ijkk(0)=zip
        qbq_ijii(0)=zip
        qbq_ijjj(0)=zip
        qbq_ijkj(0)=zip
        qbq_ijik(0)=zip
        qbq_ijkl(0)=zip
        qbq_iiij(0)=zip
        qbq_iiji(0)=zip
        qq_iiji(0)=zip
        qq_ijkj(0)=zip
        qq_ijik(0)=zip
        qq_ijjj(0)=zip
        qq_ijii(0)=zip
        qbqb_iiji(0)=zip
        qbqb_ijkj(0)=zip
        qbqb_ijik(0)=zip
        qbqb_ijjj(0)=zip
        qbqb_ijii(0)=zip
        qqbxiiij(0)=zip
        qqbxiiji(0)=zip
        qqbxijii(0)=zip
        qqbxijjj(0)=zip
        qbqxiiij(0)=zip
        qbqxiiji(0)=zip
        qbqxijii(0)=zip
        qbqxijjj(0)=zip
      endif


c--- 4-quark contribution to matrix elements
      if (Qflag) then

      do j=-nf,nf
      do k=-nf,nf

      do i=0,2
      mqq(i,j,k)=0._dp
      enddo

      if ((j > 0) .and. (k < 0)) then
        if (j .ne. -k) then
c--- Q QBAR - different flavours
           mqq(0,j,k)=facqq*Vsq(j,k)*(qqb_ijii(0)+qqb_ijjj(0))
           mqq(1,j,k)=facqq*Vsq(j,k)*(
     &             real(nf-2,dp)*qqb_ijkk(1)+qqb_ijii(1)+qqb_ijjj(1))
           mqq(2,j,k)=facqq*Vsq(j,k)*(qqb_ijii(2)+qqb_ijjj(2))
     &            +facqq*(Vsum(j)-Vsq(j,k))*qqb_ijkj(2)
     &            +facqq*(Vsum(k)-Vsq(j,k))*qqb_ijik(2)
           do i=0,2
           msqx(i,j,k,j,-j)=facqq*Vsq(j,k)*qqb_ijii(i)
           msqx(i,j,k,-k,k)=facqq*Vsq(j,k)*qqb_ijjj(i)
           msqx(i,j,k,-j,j)=facqq*Vsq(j,k)*qqbxijii(i)
           msqx(i,j,k,k,-k)=facqq*Vsq(j,k)*qqbxijjj(i)
           enddo
           do l=1,nf
           if ((l .ne. +j) .and. (l .ne. -k)) then
            msqx(1,j,k,l,-l)=facqq*Vsq(j,k)*qqb_ijkk(1)
            msqx(1,j,k,-l,l)=msqx(1,j,k,l,-l)
           endif
           if (l .ne. +j) then
             msqx(2,j,k,j,-l)=facqq*Vsq(k,l)*qqb_ijik(2)
             msqx(2,j,k,-l,j)=facqq*Vsq(k,l)*qqbxijik(2)
           endif
           if (l .ne. -k) then
             msqx(2,j,k,+l,k)=facqq*Vsq(j,-l)*qqb_ijkj(2)
             msqx(2,j,k,k,+l)=facqq*Vsq(j,-l)*qqbxijkj(2)
           endif
           enddo

         else
c--- Q QBAR - same flavours
          Vfac=0._dp
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
          do l=1,nf
          do i=0,2
          msqx(i,j,k,j,-l)=facqq*Vsq(k,l)*qqb_iiij(i)
          msqx(i,j,k,l,k)=facqq*Vsq(j,-l)*qqb_iiji(i)
          msqx(i,j,k,-l,j)=facqq*Vsq(k,l)*qqbXiiij(i)
          msqx(i,j,k,k,l)=facqq*Vsq(j,-l)*qqbXiiji(i)
          enddo
          do m=-nf,1
            if ((l .ne. j) .and. (m .ne. k)) then
            msqx(1,j,k,l,m)=facqq*Vsq(-l,-m)*qqb_ijkl(1)
            msqx(1,j,k,m,l)=facqq*Vsq(-l,-m)*qqbXijkl(1)
            endif
          enddo
          enddo
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
           do i=0,2
           msqx(i,j,k,-j,j)=facqq*Vsq(j,k)*qbq_ijjj(i)
           msqx(i,j,k,k,-k)=facqq*Vsq(j,k)*qbq_ijii(i)
           msqx(sw(i),j,k,j,-j)=facqq*Vsq(j,k)*qbqxijjj(i)
           msqx(sw(i),j,k,-k,k)=facqq*Vsq(j,k)*qbqxijii(i)
           enddo
           do l=-nf,-1
           if ((l .ne. +j) .and. (l .ne. -k)) then
            msqx(1,j,k,l,-l)=facqq*Vsq(j,k)*qbq_ijkk(1)
            msqx(sw(1),j,k,-l,l)=msqx(1,j,k,l,-l)
           endif
           if (l .ne. +j) then
             msqx(2,j,k,-l,j)=facqq*Vsq(k,l)*qbq_ijkj(2)
             msqx(sw(2),j,k,j,-l)=facqq*Vsq(k,l)*qbqxijkj(2)
           endif
           if (l .ne. -k) then
             msqx(2,j,k,k,+l)=facqq*Vsq(j,-l)*qbq_ijik(2)
             msqx(sw(2),j,k,+l,k)=facqq*Vsq(j,-l)*qbqxijik(2)
           endif
           enddo
        else
c--- QBAR Q - same flavours
          Vfac=0._dp
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
          do l=-nf,-1
          do i=0,2
          msqx(i,j,k,-l,j)=facqq*Vsq(k,l)*qbq_iiji(i)
          msqx(i,j,k,k,l)=facqq*Vsq(j,-l)*qbq_iiij(i)
          msqx(sw(i),j,k,j,-l)=facqq*Vsq(k,l)*qbqXiiji(i)
          msqx(sw(i),j,k,l,k)=facqq*Vsq(j,-l)*qbqXiiij(i)
          enddo
          do m=1,nf
            if ((l .ne. j) .and. (m .ne. k)) then
            msqx(1,j,k,m,l)=facqq*Vsq(-l,-m)*qbq_ijkl(1)
            msqx(sw(1),j,k,l,m)=facqq*Vsq(-l,-m)*qbqXijkl(1)
            endif
          enddo
          enddo
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
        do i=0,2
        msqx(i,j,k,j,j)=facqq*half*Vsq(k,-j)*qq_ijii(i)
        msqx(i,j,k,k,k)=facqq*half*Vsq(j,-k)*qq_ijjj(i)
        enddo
        do l=-nf,nf
         if (l .ne. k) then
           msqx(1,j,k,l,k)=facqq*Vsq(j,-l)*qq_ijkj(1)
           msqx(1,j,k,k,l)=facqq*Vsq(j,-l)*qq_ijjk(1)
         endif
         if (l .ne. j) then
           msqx(1,j,k,j,l)=facqq*Vsq(k,-l)*qq_ijik(1)
           msqx(1,j,k,l,j)=facqq*Vsq(k,-l)*qq_ijki(1)
         endif
        enddo
        else
c--- Q Q - same flavours
          mqq(0,j,k)=facqq*Vsum(j)*qq_iiji(0)
          mqq(1,j,k)=facqq*Vsum(j)*qq_iiji(1)
          mqq(2,j,k)=facqq*Vsum(j)*qq_iiji(2)
          do l=-nf,nf
            if (Vsq(j,-l) .ne. 0._dp) then
            do i=0,2
            msqx(i,j,k,l,j)=facqq*Vsq(j,-l)*qq_iiji(i)
            msqx(i,j,k,j,l)=facqq*Vsq(j,-l)*qq_iiij(i)
            enddo
            endif
          enddo
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
        do i=0,2
        msqx(i,j,k,j,j)=facqq*half*Vsq(k,-j)*qbqb_ijii(i)
        msqx(i,j,k,k,k)=facqq*half*Vsq(j,-k)*qbqb_ijjj(i)
        enddo
        do l=-nf,nf
         if (l .ne. k) then
           msqx(1,j,k,l,k)=facqq*Vsq(j,-l)*qbqb_ijkj(1)
           msqx(1,j,k,k,l)=facqq*Vsq(j,-l)*qbqb_ijjk(1)
         endif
         if (l .ne. j) then
           msqx(1,j,k,j,l)=facqq*Vsq(k,-l)*qbqb_ijik(1)
           msqx(1,j,k,l,j)=facqq*Vsq(k,-l)*qbqb_ijki(1)
         endif
        enddo
        else
c--- QBAR QBAR - same flavours
          mqq(0,j,k)=facqq*Vsum(j)*qbqb_iiji(0)
          mqq(1,j,k)=facqq*Vsum(j)*qbqb_iiji(1)
          mqq(2,j,k)=facqq*Vsum(j)*qbqb_iiji(2)
          do l=-nf,nf
            if (Vsq(j,-l) .ne. 0._dp) then
            do i=0,2
            msqx(i,j,k,l,j)=facqq*Vsq(j,-l)*qbqb_iiji(i)
            msqx(i,j,k,j,l)=facqq*Vsq(j,-l)*qbqb_iiij(i)
            enddo
            endif
          enddo
        endif
      endif
      if     (colourchoice == 1) then
        mqq(0,j,k)=0._dp
      elseif (colourchoice == 2) then
        mqq(1,j,k)=0._dp
        mqq(2,j,k)=0._dp
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
        msqx_cs(i,j,k)=0._dp
      enddo

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*qqbWgg2
          do i=0,2
            msqx_cs(i,j,k)=Vsq(j,k)*qqbWgg2_cs(i)
            msqx(i,j,k,0,0)=Vsq(j,k)*qqbWgg2_cs(i)
          enddo
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*qbqWgg2
          do i=0,2
            msqx_cs(i,j,k)=Vsq(j,k)*qbqWgg2_cs(i)
            msqx(i,j,k,0,0)=Vsq(j,k)*qbqWgg2_cs(i)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWqg2
          do i=0,2
            msqx_cs(i,j,k)=
     &(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWqg2_cs(i)
            do n1=1,nf
            msqx(i,j,k,n1,0)=Vsq(j,-n1)*qgWqg2_cs(i)
            msqx(i,j,k,0,n1)=Vsq(j,-n1)*qgWgq2_cs(i)
            enddo
           enddo
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqbg2
          do i=0,2
            msqx_cs(i,j,k)=
     &(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqbg2_cs(i)
            do n1=1,nf
            msqx(i,j,k,-n1,0)=Vsq(j,n1)*qbgWqbg2_cs(i)
            msqx(i,j,k,0,-n1)=Vsq(j,n1)*qbgWgqb2_cs(i)
            enddo
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWqg2
          do i=0,2
            msqx_cs(i,j,k)=
     &(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWqg2_cs(i)
            do n1=1,nf
            msqx(i,j,k,n1,0)=Vsq(-n1,k)*gqWqg2_cs(i)
            msqx(i,j,k,0,n1)=Vsq(-n1,k)*gqWgq2_cs(i)
            enddo
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=msq(j,k)+
     &(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqbg2
          do i=0,2
            msqx_cs(i,j,k)=
     &(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqbg2_cs(i)
            do n1=1,nf
            msqx(i,j,k,-n1,0)=Vsq(n1,k)*gqbWqbg2_cs(i)
            msqx(i,j,k,0,-n1)=Vsq(n1,k)*gqbWgqb2_cs(i)
            enddo
           enddo
      elseif ((j == 0) .and. (k == 0)) then
          Vfac=0._dp
          do n1=1,nf
            do n2=-nf,-1
              Vfac=Vfac+Vsq(n1,n2)
              do i=0,2
              msqx(i,j,k,-n1,-n2)=Vsq(n1,n2)*ggWqbq2_cs(i)
              msqx(i,j,k,-n2,-n1)=Vsq(n1,n2)*ggWqqb2_cs(i)
              enddo
            enddo
          enddo
          msq(j,k)=msq(j,k)+Vfac*ggWqbq2
          do i=0,2
            msqx_cs(i,j,k)=Vfac*ggWqbq2_cs(i)
          enddo
      endif

      enddo
      enddo

      endif

c--- restore proper colourchoice if necessary
      if ((kpart==kreal) .or. (kpart==kvirt)) then
        colourchoice=rcolourchoice
      endif
c--- restore proper parton sub-process selection, if necessary
      if ((kpart==kreal) .and. (Qflag .eqv. .true.)) then
        Gflag=rGflag
      endif

      return
      end

