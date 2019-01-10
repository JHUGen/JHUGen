      subroutine findmind(p,pjet,pjetmin,pjetmax,dijmin,nmin1,nmin2,
     .                    ipow)
c--- this finds the minimum dij for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),dijmin,dij,d
      integer pjetmin,pjetmax,nmin1,nmin2,i,j,ipow

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
c--- this finds the minimum dkmin for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
C--- calculate the beam proto-jet separation see NPB406(1993)187, Eqn. 7
C---  S.~Catani, Y.~L.~Dokshitzer, M.~H.~Seymour and B.~R.~Webber
C--- in  practice this is just the minimum ptsq of protojets
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),dkmin,dk,pt
      integer pjetmin,pjetmax,nk,i,ipow
      logical dkerror

      dkmin=1d9
      dkerror=.true.

      do i=pjetmin,pjetmax
        dk=pt(i,pjet)
        if (ipow .ne. 1) dk=dk**(ipow)
        if (dk .lt. dkmin) then
          dkmin=dk
          nk=i
          dkerror=.false.
        endif
      enddo

      if (dkerror) then
        write(*,*) 'Error in dk minimum-finding routine'
        stop
      endif

      return
      end

      double precision function dij(p,pjet,i,j,ipow)
C---calculate the proto-jet separation see NPB406(1993)187, Eqn. 7
      implicit none
      include 'constants.f'
      integer i,j,ipow
      double precision p(mxpart,4),pjet(mxpart,4),pti,ptj,pt,r
c      double precision etarap,yi,yj,phii,phij

      pti=pt(i,pjet)
      ptj=pt(j,pjet)

c--- old method - bad because (phii-phij) can be > pi
c      yi=etarap(i,pjet)
c      yj=etarap(j,pjet)

c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))

c      dij=dsqrt((yi-yj)**2+(phii-phij)**2)

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
      include 'constants.f'
      include 'jetlabel.f'
      include 'process.f'
      integer i,j
      double precision pjet(mxpart,4)

c--Run II prescription
      pjet(i,1)=pjet(i,1)+pjet(j,1)
      pjet(i,2)=pjet(i,2)+pjet(j,2)
      pjet(i,3)=pjet(i,3)+pjet(j,3)
      pjet(i,4)=pjet(i,4)+pjet(j,4)

c--- special combination tag for W+heavy quarks
      if ((case .eq. 'Wbbmas') .or. (case .eq. 'W_bjet')) then
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
        jetlabel(i)='bb'
        return
      endif
      if (((jetlabel(i) .eq. 'bb') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'bb') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='bb'
        return
      endif
      endif

      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
        jetlabel(i)='pp'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'qj'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'qj'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) .eq. 'qj') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'pp') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='qj'
        return
      endif

      return
      end

c      subroutine combine_snowmass(p,pjet,i,j)
c      implicit none
c      include 'constants.f'
c      include 'jetlabel.f'
c      integer i,j
c      double precision p(mxpart,4),pjet(mxpart,4),ptjetij,yjet,phijet,
c     . ejet,pt,etarap,pti,ptj,yi,yj,phii,phij

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

c      pjet(i,1)=ptjetij*dsin(phijet)
c      pjet(i,2)=ptjetij*dcos(phijet)
c      pjet(i,3)=ptjetij*(ejet-1d0/ejet)/2d0
c      pjet(i,4)=ptjetij*(ejet+1d0/ejet)/2d0


c      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'pp'))
c     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'pp'))) then
c        jetlabel(i)='bq'
c        return
c      endif
c      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'pp'))
c     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'pp'))) then
c        jetlabel(i)='ba'
c        return
c      endif
c      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
c     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
c        jetlabel(i)='pp'
c        return
c      endif
c
c      return
c      end

c      subroutine shuffle(pjet,nmin,nmax)
c--- shuffles jets nmin..nmax-1 in pjet down by 1 index
c      implicit none
c      include 'constants.f'
c      integer i,j,nmin,nmax
c      double precision pjet(mxpart,4)

c      if (nmin .eq. nmax) return

c      do i=nmin,nmax-1
c        do j=1,4
c          pjet(i,j)=pjet(i+1,j)
c        enddo
c      enddo

c      return
c      end

      subroutine swap(pjet,i,j)
c--- swaps jets i..j in pjet
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      integer i,j,k
      double precision pjet(mxpart,4),tmp
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

c      double precision function ptjet(j,p,pjet)
c      implicit none
c      include 'constants.f'
c      integer j
c      double precision p(mxpart,4),pjet(mxpart,4)
c--- This is the formula for pt
c      ptjet=dsqrt(pjet(j,1)**2+pjet(j,2)**2)
c--- This is the formula for Et
c      ptjet=dsqrt(pjet(j,1)**2+pjet(j,2)**2)
c     . *pjet(j,4)/dsqrt(pjet(j,1)**2+pjet(j,2)**2+pjet(j,3)**2)
c      return
c      end

      double precision function dotjet(p,i,pjet,j)
C---Dot the ith vector p with the jth vector pjet
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),pjet(mxpart,4)

      dotjet=p(i,4)*pjet(j,4)-p(i,1)*pjet(j,1)
     .      -p(i,2)*pjet(j,2)-p(i,3)*pjet(j,3)

      return
      end

      double precision function bclustmass(pjet)
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      integer i,nbq,nba
      double precision pjet(mxpart,4)

c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      bclustmass=0d0
      nbq=0
      nba=0

      do i=1,jets
        if (jetlabel(i) .eq. 'bq') nbq=i+4
        if (jetlabel(i) .eq. 'ba') nba=i+4
      enddo

      if ((nbq .eq. 0) .or. (nba .eq. 0)) return

      bclustmass=(pjet(nbq,4)+pjet(nba,4))**2
      do i=1,3
        bclustmass=bclustmass-(pjet(nbq,i)+pjet(nba,i))**2
      enddo

      return
      end

