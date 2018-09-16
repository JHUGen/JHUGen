      function gencuts_VHWW(p)
      implicit none
      include 'types.f'
      logical:: gencuts_VHWW
!======C.W Nov 2015, basic cut to do triboson style VH=>VWW signitures
!======cuts on observed leptons pt's and rapidity and MET
!======the transverse mass of the VH system is calculated and put in a common
!======block for the histograming routine
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'first.f'
      include 'leptcuts.f'
      include 'nwz.f'
      include 'mpicommon.f'
      real(dp) :: p(mxpart,4)
      real(dp) :: ptlepmin,ylmax,metmin,mVHtrans
      real(dp) :: metvec(4),mymet
      real(dp) :: lepvec(4)
      real(dp) m3lsq,pt3lsq,etmiss,pt3lmiss,et3l
      real(dp) :: pt,etarap,pl1,pl2,pl3,pl4
      integer j
      integer lepid1,lepid2,lepid3,lepid4
      integer metid1,metid2,metid3
      common/observables_VHWW/mVHtrans,mymet,pl1,pl2,pl3,pl4
!$omp threadprivate(/observables_VHWW/)

      gencuts_VHWW=.false.
!======lepton cuts from input
      ptlepmin=leptpt
      ylmax=leptrap
      metmin=misspt
      metvec(:)=zip
      lepvec(:)=zip
      pl1=-1._dp
      pl2=-1._dp
      pl3=-1._dp
      pl4=-1._dp

      mVHtrans=zip
      if(first) then
         first=.false.
         if (rank == 0) then
         write(6,*) '****************** VH=>WW cuts ********************'
         write(6,99) '*        pt(lepton)      >   ',ptlepmin,
     &                ' GeV            *'
         write(6,99) '*        |eta(lepton)|      <   ',ylmax,
     &     '            *'
         write(6,99) '*        MET      >   ',metmin,
     &                ' GeV            *'
         write(6,*) '***************************************************'
         endif
      endif

      if((kcase.eq.kWH__WW).or.(kcase.eq.kWH1jet)) then
!====== W+ H
         if(nwz==1) then
            lepid1=4
            metid1=3
         elseif(nwz==-1) then
            lepid1=3
            metid1=4
         endif
         lepid2=6
         lepid3=7
         do j=1,4
            metvec(j)=p(metid1,j)+p(5,j)+p(8,j)
            lepvec(j)=p(lepid1,j)+p(lepid2,j)+p(lepid3,j)
         enddo
!======check leptons
         if((pt(lepid1,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid1,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         if((pt(lepid2,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid2,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         if((pt(lepid3,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid3,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         mymet=sqrt(metvec(1)**2+metvec(2)**2)
         if(mymet.lt.metmin) then
            gencuts_VHWW=.true.
            return
         endif

         pl1=pt(lepid1,p)
         pl2=pt(lepid2,p)
         pl3=pt(lepid3,p)

      elseif((kcase.eq.kZH__WW).or.(kcase.eq.kZH1jet)) then
         lepid1=3
         lepid2=4
         lepid3=6
         lepid4=7
         do j=1,4
            metvec(j)=p(5,j)+p(8,j)
            lepvec(j)=p(lepid1,j)+p(lepid2,j)+p(lepid3,j)+p(lepid4,j)
         enddo
!======check leptons
         if((pt(lepid1,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid1,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         if((pt(lepid2,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid2,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         if((pt(lepid3,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid3,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         if((pt(lepid4,p).lt.ptlepmin)
     &        .or.(abs(etarap(lepid4,p)).gt.ylmax)) then
            gencuts_VHWW=.true.
            return
         endif
         mymet=sqrt(metvec(1)**2+metvec(2)**2)
         if(mymet.lt.metmin) then
            gencuts_VHWW=.true.
            return
         endif


         pl1=pt(lepid1,p)
         pl2=pt(lepid2,p)
         pl3=pt(lepid3,p)
         pl4=pt(lepid4,p)

      endif


!=======calculate Transverse mass
      m3lsq=lepvec(4)**2-lepvec(3)**2-lepvec(2)**2-lepvec(1)**2
      pt3lsq=lepvec(1)**2+lepvec(2)**2

      et3l=sqrt(pt3lsq+m3lsq)
      etmiss=sqrt(metvec(2)**2+metvec(1)**2)

      pt3lmiss=(lepvec(1)+metvec(1))**2+(lepvec(2)+metvec(2))**2

      mVHtrans=sqrt((et3l+etmiss)**2-pt3lmiss)

      return
 99   format(1x,a29,f6.2,a17)

      end


