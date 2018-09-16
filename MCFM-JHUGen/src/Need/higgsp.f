      subroutine higgsp(bbbr,tautaubr,gamgambr,zgambr,wwbr,zzbr)
      implicit none
      include 'types.f'
C--returns Higgs branching ratios as calculated by
C  interpolating the Spira tables br.sm1 br.sm2
C  Other branching ratios could be added.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'first.f'
      integer:: npt
      parameter(npt=1000)
      integer:: j,nlo
      character*79 string
      real(dp):: bbbr,gamgambr,tautaubr,zgambr,wwbr,zzbr,htemp,
     & width(npt)
      real(dp):: xmh(npt),brbb(npt),brtautau(npt),brss,brcc,brmumu
     & ,brtt,brgg,brgamgam(npt),brzgam(npt),brww(npt),brzz(npt)
      save brbb,brww,brzz,width

      if (first) then
      first=.false.
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

      if (hmass < one) then
      htemp=one
      nlo=1
      elseif (hmass > 999._dp) then
      htemp=999._dp
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
