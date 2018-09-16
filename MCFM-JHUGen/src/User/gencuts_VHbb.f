!---- CW. Oct 2015
      function gencuts_VHbb(p)
!=====cuts for VH=>l_1,l_2 + b + bbar

!======for the NNLO calculation the process is setup to treat the
!=====bb pair as EW particles (i.e. no clustering).
!=====however for cuts we want to treat the b as jets.
!=====This routine therefore calls a second clustering routine, and defines
!==== new jet phase space, and observables to be histogrammed, therefore there
!==== is an intimate relationship between the observable and the cutting routine

      implicit none
      include 'types.f'
      include 'mpicommon.f'
      logical gencuts_VHbb
      include 'constants.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'nwz.f'
      include 'jetcuts.f'
      include 'clustering.f'
      include 'jetlabel.f'
      include 'first.f'
      include 'energy.f'
      real(kind=dp):: p(mxpart,4),pjet(mxpart,4)
      logical failed
      real(dp) :: pt,yrap
      real(dp) :: mtrans,mt34max,ptVmin,ptlmin,metmin
      real(dp) :: ylmax,m34,m34min,m34max,ptV,ayrap
      real(dp) :: pttwo,yraptwo,R
      real(dp) :: mbb,mVbb,ptbb,ptbh,ptbs,ybb
      integer lepid,metid,j,njets
      integer nj,nbj,jets_orig,noj
      integer njetmax,nbreq,idbmax
      integer :: bid(mxpart)
      real(dp) :: ptbmax,rbbmin,mtW,mtWH
      common/njetsVH/nj,nbj,noj
      common/observables_VHbb/ptV,mbb,mVbb,ptbb,ptbh,ptbs,mtW,mtWH
!$omp threadprivate(/njetsVH/,/observables_VHbb/)
      real(dp):: yr4_ptH(0:3), yr4_ptl1(0:3),yr4_ptl2(0:3)
      real(dp):: yr4_yH(0:3), yr4_yl1(0:3),yr4_yl2(0:3)
      common/yr4_observables/yr4_ptH,yr4_ptl1,yr4_ptl2,yr4_yH
     &     ,yr4_yl1,yr4_yl2
!$omp threadprivate(/yr4_observables/)
      integer itag

      jets_orig=jets
      noj=jets_orig
      nj=jets
      nbj=0
      gencuts_VHbb=.false.
!      goto 999
      bid(:)=0

!======= lepton cuts
      ptVmin=0._dp
      ptlmin=15._dp
      ylmax=2.5_dp
      m34min=75._dp
      m34max=105._dp
      metmin=15._dp
      mt34max=sqrts
      mtW=-1._dp
      mtWH=-1._dp

      if(first) then
         first=.false.
         call read_jetcuts(ptjetmin,etajetmin,etajetmax)
         if (rank == 0) then
      write(6,*)  '****************** VH=>bb cuts ********************'
      write(6,99) '*        pt(lepton)      >   ',ptlmin,
     &                ' GeV            *'
      write(6,99) '*        |eta(lepton)|      <   ',ylmax,
     &     '            *'
      write(6,99) '*        pt(jet)      >   ',ptjetmin,
     &                ' GeV            *'
      write(6,99) '*        |eta(jet)|      <   ',etajetmax,
     &     '            *'
      write(6,99) '*        pt(Vector Boson)   >   ',ptVmin,
     &     ' GeV            *'
      if((kcase.eq.kWHbbar).or.(kcase.eq.kWH1jet)) then
         write(6,99) '*        mtrans)   <   ',mt34max,
     &        ' GeV            *'
         write(6,99) '*        MET      >   ',metmin,
     &                ' GeV            *'

      endif
      write(6,*)  '****************************** ********************'
         endif
      endif


!======jet cuts

!=   min # of b's required (<0 = inclusive)
      nbreq=2
!=   # maximum number of jets (  >= 4 = inclusive at NNLO)
      njetmax=4


!=====lepton cuts W
      if((kcase.eq.kWHbbar).or.(kcase.eq.kWH1jet)
     &     .or.(kcase.eq.kWHbbdk)) then
         if(nwz==1) then
!=========lepton is p4
            lepid=4
            metid=3
         elseif(nwz==-1) then
!=========lepton is p3
            lepid=4
            metid=3
         else
            write(6,*) 'unrecognized nwz =',nwz
            stop
         endif

!=========lepton is p4
        if((pt(lepid,p).lt.ptlmin).or.(abs(yrap(lepid,p)).gt.ylmax))then
            gencuts_VHbb=.true.
            goto 999
         endif
!=========met is p3
         if(pt(metid,p).lt.metmin) then
            gencuts_VHbb=.true.
            goto 999
         endif


!===== transverse mass cut
         mtrans=
     &   (p(3,1)*p(4,1)+p(3,2)*p(4,2))
     &   /sqrt((p(3,1)**2+p(3,2)**2)
     &         *(p(4,1)**2+p(4,2)**2))
        mtrans=2._dp*sqrt(p(3,1)**2+p(3,2)**2)
     &   *sqrt(p(4,1)**2+p(4,2)**2)*(1._dp-mtrans)
        mtrans=sqrt(max(mtrans,zip))
        if (mtrans > abs(mt34max)) then
            gencuts_VHbb=.true.
            goto 999
         endif
         mtW=mtrans

!======transverse mass of H,V system, not used for cuts but for histo

      elseif((kcase.eq.kZHbbar).or.(kcase.eq.kZH1jet)
     &       .or.(kcase.eq.kZHbbdk)) then
!======check both leptons and m34
         do lepid=3,4
         if((pt(lepid,p).lt.ptlmin).or.(abs(yrap(lepid,p)).gt.ylmax))then
            gencuts_VHbb=.true.
            goto 999
         endif
         enddo
         m34=(p(3,4)+p(4,4))**2
         do j=1,3
            m34=m34-(p(3,j)+p(4,j))**2
         enddo
         m34=sqrt(m34)
         if((m34.lt.m34min).or.(m34.gt.m34max)) then
            gencuts_VHbb=.true.
            goto 999
         endif
      endif

!=====PT_V cut
      ptV=zip
      do j=1,2
         ptV=ptV+(p(3,j)+p(4,j))**2
      enddo
      ptV=sqrt(ptV)
      if(ptV.lt.ptVmin) then
         gencuts_VHbb=.true.
         goto 999
      endif


      nbj=0
      ptbmax=0._dp
      idbmax=-1
!======jet cuts, first count bjets and collect hardest pt
!      if(bbclust) then
!==== bs are in jet label, so do normal cuts
         do j=1,jets
            if((jetlabel(j)=='bq'.or.jetlabel(j)=='ba')
     &     .or.(jetlabel(j)=='qb'.or.jetlabel(j)=='ab')) then
               nbj=nbj+1
               bid(nbj)=j+4
               if(pt(j+4,p).gt.ptbmax) then
                  ptbmax=pt(j+4,p)
                  idbmax=j+4
               endif
            endif
         enddo
         nj=jets

!      else
!==== jets are either 0,1,2 from QCD radiation, dont cluster bs and simply check
!==== pts of bs and Rbb against input jet parameters they are always in 5 and 6

!==== are the b's close enough to be clustered?
!         if(R(p,5,6).lt.rbbmin) then
!=========treat as one jet
!            if(pttwo(5,6,p).gt.ptjetmin
!     &           .and.(abs(yraptwo(5,6,p)).lt.etajetmax))then
!               nbj=nbj+1
!               ptbmax=pttwo(5,6,p)
!            endif
!         else
!            do j=5,6
!               if((pt(j,p).gt.ptjetmin).and.(ayrap(j,p).lt.etajetmax))
!     &              then
!                  nbj=nbj+1
!                  bid(nbj)=j
!                  if(pt(j,p).gt.ptbmax) then
!                     ptbmax=pt(j,p)
!                     idbmax=j
!                  endif
!      endif
!            enddo
!         endif
!         nj=jets+nbj
!      endif



!      goto 999
!========cut on b-jets require exactly nbrec number of b jets
!========and require that hardest p > ptmin_bjet (read from input.DAT)
!====== to make the code inclusive in b-jets set nbreq < 0
      if((nbreq.ge.0).and.(nbj.lt.nbreq)) then
         gencuts_VHbb=.true.
         goto 999
      endif
      if(ptbmax < ptbjetmin) then
         gencuts_VHbb=.true.
         goto 999
      endif

!=====optionally veto nj > nmaxj additional jets
      if(nj > njetmax) then
         gencuts_VHbb=.true.
         goto 999
      endif

!===== next define b-jet observables for histogramming

!======mbb
      mbb=zip
      mvbb=zip

      if(nbj==0) then
         ptbh=-1._dp
         ptbs=-1._dp
      endif

      if(nbj==1) then
         ptbh=ptbmax
         ptbs=-1._dp
      endif

      if(nbj==2) then
         mbb=(p(bid(1),4)+p(bid(2),4))**2
         mVbb=(p(bid(1),4)+p(bid(2),4)+p(3,4)+p(4,4))**2
         do j=1,3
            mbb=mbb-(p(bid(1),j)+p(bid(2),j))**2
            mVbb=mVbb
     &           -(p(bid(1),j)+p(bid(2),j)+p(3,4)+p(4,j))**2
         enddo
         mbb=sqrt(mbb)
         mVbb=sqrt(mVbb)

!     = ptbb
         ptbb= (p(bid(1),1)+p(bid(2),1))**2
     &        +(p(bid(1),2)+p(bid(2),2))**2
         ptbb=sqrt(ptbb)
         if(pt(bid(1),p).gt.pt(bid(2),p)) then
            ptbh=pt(bid(1),p)
            ptbs=pt(bid(2),p)
         else
            ptbh=pt(bid(2),p)
            ptbs=pt(bid(1),p)
         endif
         ybb=yraptwo(bid(1),bid(2),p)
      else
!======= set out of range for histogramming
         mbb=-1._dp
         mvbb=-1._dp
         ptbb=-1._dp
         ybb=-1000._dp
      endif


!-----initialize
      yr4_ptH(:) = -1._dp
      yr4_yH(:) = -1000._dp
      yr4_ptl1(:) = -1._dp
      yr4_ptl2(:) = -1._dp
      yr4_yl1(:) = -10000._dp

!-------now calualte observables for histograms for YR4
!------determine ptV
      if(ptV.lt.150._dp) then
         itag=1
      elseif((ptV.gt.150._dp).and.(ptV.lt.250._dp)) then
         itag=2
      elseif(ptV.gt.250._dp) then
         itag=3
      endif

!=======ptH
      yr4_ptH(0)=ptbb
      yr4_ptH(itag)=ptbb
!======= yH
      yr4_yH(0)=ybb
      yr4_yH(0)=ybb
!=======ptl1
      yr4_ptl1(0)=pt(3,p)
      yr4_ptl1(itag)=pt(3,p)
!=======ptl2
      yr4_ptl2(0)=pt(4,p)
      yr4_ptl2(itag)=pt(4,p)
!=======yl1
      yr4_yl1(0)=yrap(3,p)
      yr4_yl1(itag)=yrap(3,p)
!=======yl2
      yr4_yl2(0)=yrap(4,p)
      yr4_yl2(itag)=yrap(4,p)




!=====reset jets to not disrubt any future routines
      jets=jets_orig
      return
 99   format(1x,a29,f6.2,a17)
 999  continue
      jets=jets_orig
      return
      end

