      subroutine findmind(p,pjet,pjetmin,pjetmax,dijmin,nmin1,nmin2,
     &                    ipow)
      implicit none
      include 'types.f'
c--- this finds the minimum dij for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),pjet(mxpart,4),dijmin,dij,d
      integer:: pjetmin,pjetmax,nmin1,nmin2,i,j,ipow

      do i=pjetmin,pjetmax
        do j=i+1,pjetmax
          d=dij(p,pjet,i,j,ipow)
          if ((i == pjetmin) .and. (j == i+1)) then
            dijmin=d
            nmin1=i
            nmin2=j
          elseif (d < dijmin) then
            dijmin=d
            nmin1=i
            nmin2=j
          endif
        enddo
      enddo

      return
      end

      subroutine findminet(p,pjet,pjetmin,pjetmax,dkmin,nk,ipow)
      implicit none
      include 'types.f'
c--- this finds the minimum dkmin for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
C--- calculate the beam proto-jet separation see NPB406(1993)187, Eqn. 7
C---  S.~Catani, Y.~L.~Dokshitzer, M.~H.~Seymour and B.~R.~Webber
C--- in practice this is just the minimum ptsq of protojets
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),pjet(mxpart,4),dkmin,dk,pt
      integer:: pjetmin,pjetmax,nk,i,ipow

      dkmin=pt(pjetmin,pjet)
      if (ipow .ne. 1) dkmin=dkmin**(ipow)
      nk=pjetmin

c--- if only one entry, this must be the minimum
      if (pjetmin+1 > pjetmax) return

      do i=pjetmin+1,pjetmax
        dk=pt(i,pjet)
        if (ipow .ne. 1) dk=dk**(ipow)
        if (dk < dkmin) then
          dkmin=dk
          nk=i
        endif
      enddo

      return
      end

      function dij(p,pjet,i,j,ipow)
      implicit none
      include 'types.f'
      real(dp):: dij
C---calculate the proto-jet separation see NPB406(1993)187, Eqn. 7

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i,j,ipow
      real(dp):: p(mxpart,4),pjet(mxpart,4),pti,ptj,pt,r
c      real(dp):: etarap,yi,yj,phii,phij

      pti=pt(i,pjet)
      ptj=pt(j,pjet)

c--- old method - bad because (phii-phij) can be > pi
c      yi=etarap(i,pjet)
c      yj=etarap(j,pjet)

c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))

c      dij=sqrt((yi-yj)**2+(phii-phij)**2)

c--- new method - r() calculates true value of 0 < (phi-phij) < pi
      dij=r(pjet,i,j)

      if (ipow .ne. 1) then
        pti=pti**(ipow)
        ptj=ptj**(ipow)
      endif
      dij=dij*min(pti,ptj)

      return
      end

      subroutine combine(pjet,i,j)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      include 'kprocess.f'
      integer:: i,j
      real(dp):: pjet(mxpart,4)

c--Run II prescription
      pjet(i,1)=pjet(i,1)+pjet(j,1)
      pjet(i,2)=pjet(i,2)+pjet(j,2)
      pjet(i,3)=pjet(i,3)+pjet(j,3)
      pjet(i,4)=pjet(i,4)+pjet(j,4)

c--- special combination tag for W+heavy quarks
      if ((kcase==kWbbmas) .or. (kcase==kW_bjet)
     & .or.(kcase==kWHbbar) .or. (kcase==kZHbbar)
     &) then
      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'ba'))
     ..or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'ba'))) then
        jetlabel(i)='bb'
        return
      endif
      if (((jetlabel(i) == 'bb') .and. (jetlabel(j) == 'pp'))
     ..or.((jetlabel(j) == 'bb') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='bb'
        return
      endif
      endif

      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'pp'))
     ..or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) == 'ba') .and. (jetlabel(j) == 'pp'))
     ..or.((jetlabel(j) == 'ba') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'ba'))
     ..or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'ba'))) then
        jetlabel(i)='pp'
        return
      endif
      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'qj'))
     ..or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) == 'ba') .and. (jetlabel(j) == 'qj'))
     ..or.((jetlabel(j) == 'ba') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) == 'qj') .and. (jetlabel(j) == 'pp'))
     ..or.((jetlabel(j) == 'pp') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='qj'
        return
      endif

      return
      end

c      subroutine combine_snowmass(p,pjet,i,j)
c      implicit none
c      include 'types.f'
c
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'jetlabel.f'
c      integer:: i,j
c      real(dp):: p(mxpart,4),pjet(mxpart,4),ptjetij,yjet,phijet,
c     & ejet,pt,etarap,pti,ptj,yi,yj,phii,phij

C----Snowmass style prescripton
c      pti=pt(i,pjet)
c      ptj=pt(j,pjet)

c      yi=etarap(i,pjet)
c      yj=etarap(j,pjet)

c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))
c

c      ptjetij=pti+ptj
c      yjet=(pti*yi+ptj*yj)/ptjetij
c      phijet=(pti*phii+ptj*phij)/ptjetij
c      ejet=exp(yjet)

c      pjet(i,1)=ptjetij*sin(phijet)
c      pjet(i,2)=ptjetij*cos(phijet)
c      pjet(i,3)=ptjetij*(ejet-1._dp/ejet)/2._dp
c      pjet(i,4)=ptjetij*(ejet+1._dp/ejet)/2._dp


c      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'pp'))
c     ..or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'pp'))) then
c        jetlabel(i)='bq'
c        return
c      endif
c      if (((jetlabel(i) == 'ba') .and. (jetlabel(j) == 'pp'))
c     ..or.((jetlabel(j) == 'ba') .and. (jetlabel(i) == 'pp'))) then
c        jetlabel(i)='ba'
c        return
c      endif
c      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'ba'))
c     ..or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'ba'))) then
c        jetlabel(i)='pp'
c        return
c      endif
c
c      return
c      end

c      subroutine shuffle(pjet,nmin,nmax)
c      implicit none
c      include 'types.f'
c--- shuffles jets nmin..nmax-1 in pjet down by 1 index
c
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      integer:: i,j,nmin,nmax
c      real(dp):: pjet(mxpart,4)

c      if (nmin == nmax) return

c      do i=nmin,nmax-1
c        do j=1,4
c          pjet(i,j)=pjet(i+1,j)
c        enddo
c      enddo

c      return
c      end

      subroutine swap(pjet,i,j)
      implicit none
      include 'types.f'
c--- swaps jets i..j in pjet

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,j,k
      real(dp):: pjet(mxpart,4),tmp
      character*2 chartmp

      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo

      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp

      return
      end

c      function ptjet(j,p,pjet)
c      implicit none
c      include 'types.f'
c      real(dp):: ptjet
c
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      integer:: j
c      real(dp):: p(mxpart,4),pjet(mxpart,4)
c--- This is the formula for pt
c      ptjet=sqrt(pjet(j,1)**2+pjet(j,2)**2)
c--- This is the formula for Et
c      ptjet=sqrt(pjet(j,1)**2+pjet(j,2)**2)
c     & *pjet(j,4)/sqrt(pjet(j,1)**2+pjet(j,2)**2+pjet(j,3)**2)
c      return
c      end

      function dotjet(p,i,pjet,j)
      implicit none
      include 'types.f'
      real(dp):: dotjet
C---Dot the ith vector p with the jth vector pjet

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i,j
      real(dp):: p(mxpart,4),pjet(mxpart,4)

      dotjet=p(i,4)*pjet(j,4)-p(i,1)*pjet(j,1)
     &      -p(i,2)*pjet(j,2)-p(i,3)*pjet(j,3)

      return
      end

      function bclustmass(pjet)
      implicit none
      include 'types.f'
      real(dp):: bclustmass

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,nbq,nba
      real(dp):: pjet(mxpart,4)

c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      bclustmass=0._dp
      nbq=0
      nba=0

      do i=1,jets
        if (jetlabel(i) == 'bq') nbq=i+4
        if (jetlabel(i) == 'ba') nba=i+4
      enddo

      if ((nbq == 0) .or. (nba == 0)) return

      bclustmass=(pjet(nbq,4)+pjet(nba,4))**2
      do i=1,3
        bclustmass=bclustmass-(pjet(nbq,i)+pjet(nba,i))**2
      enddo

      return
      end

