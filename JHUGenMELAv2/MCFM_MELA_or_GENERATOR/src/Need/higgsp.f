      subroutine higgsp(bbbr,tautaubr,gamgambr,zgambr,wwbr,zzbr)
C--returns Higgs branching ratios as calculated by
C  interpolating the Spira tables br.sm1 br.sm2
C  Other branching ratios could be added.
      implicit none
      include 'masses.f'
      !include 'first.f'
      logical first_higgsp
      data first_higgsp/.true./
      save first_higgsp
!$omp threadprivate(first_higgsp)
      integer npt
      parameter(npt=1000)
      integer j,nlo
      character*79 string
      double precision bbbr,gamgambr,tautaubr,zgambr,wwbr,zzbr,htemp,
     . width(npt)
      double precision xmh(npt),brbb(npt),brtautau(npt),brss,brcc,brmumu
     . ,brtt,brgg,brgamgam(npt),brzgam(npt),brww(npt),brzz(npt)
      save brbb,brww,brzz,width

      if (first_higgsp) then
      first_higgsp=.false.
      open(unit=47,file='br.sm1',status='old',err=44)
      read(47,*,err=75) string
 75   read(47,*,err=76) string
 76   continue
      do j=1,npt
      read(47,*,err=77) xmh(j),brbb(j),brtautau(j),brmumu,brss,brcc,brtt
      enddo
 77   continue
      close(unit=47)

      open(unit=48,file='br.sm2',status='old',err=45)
      read(48,*,err=85) string
 85   read(48,*,err=86) string
 86   continue
      do j=1,npt
      read(48,*) xmh(j),
     &           brgg,brgamgam(j),brzgam(j),brww(j),brzz(j),width(j)
      enddo
      close(unit=48)
      endif

      if (hmass .lt. 1d0) then
      htemp=1d0
      nlo=1
      elseif (hmass .gt. 999d0) then
      htemp=999d0
      nlo=999
      else
      htemp=hmass
      nlo=int(htemp)
      endif
      bbbr=brbb(nlo)+(htemp-nlo)*(brbb(nlo+1)-brbb(nlo)) 
      tautaubr=brtautau(nlo)+(htemp-nlo)*(brtautau(nlo+1)-brtautau(nlo))
      gamgambr=brgamgam(nlo)+(htemp-nlo)*(brgamgam(nlo+1)-brgamgam(nlo))
      zgambr=brzgam(nlo)+(htemp-nlo)*(brzgam(nlo+1)-brzgam(nlo))
      wwbr=brww(nlo)+(htemp-nlo)*(brww(nlo+1)-brww(nlo))
      zzbr=brzz(nlo)+(htemp-nlo)*(brzz(nlo+1)-brzz(nlo))
      hwidth=width(nlo)+(htemp-nlo)*(width(nlo+1)-width(nlo))
      return

 44   continue
      write(6,*) 'Error opening br1.sm1'
      stop

 45   continue
      write(6,*) 'Error opening br2.sm2'
      stop
      return
      end
