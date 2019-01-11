      subroutine initialize
!******************************************************************
!     sets up masses and coupling constants for Helas
!******************************************************************
      implicit none
!
! Constants
!
      double precision  sw2
c      parameter        (sw2 = 0.2312d0)

!   masses and widths of fermions

      double precision tmass,      bmass,     cmass
      parameter       (tmass=175d0,bmass=0d0, cmass=0d0)
      double precision smass,      umass,     dmass
      parameter       (smass=0d0,  umass=0d0, dmass=0d0)
      double precision twidth,     bwidth,    cwidth
      parameter       (twidth=0d0, bwidth=0d0,cwidth=0d0)
c      parameter       (twidth=1.5525210d0, bwidth=0d0,cwidth=0d0)
      double precision swidth,     uwidth,    dwidth
      parameter       (swidth=0d0, uwidth=0d0,dwidth=0d0)
      double precision emass,      mumass,    taumass
      parameter       (emass=0d0,  mumass=0d0,taumass=0d0)
      double precision ewidth,     muwidth,    tauwidth
      parameter       (ewidth=0d0, muwidth=0d0,tauwidth=0d0)

!   masses and widths of bosons

      double precision wmass,         zmass
      parameter       (wmass=80.419d0,    zmass=91.188d0)
      double precision wwidth,        zwidth
c      parameter       (wwidth=2.06d0, zwidth=2.49d0)
      parameter       (wwidth=0d0, zwidth=0d0)
      double precision amass,         awidth
      parameter       (amass=0d0,     awidth=0d0)
      double precision hmass,         hwidth
      parameter       (hmass=100d0,   hwidth=2.4251446D-10)
!
! Local Variables
!
      integer i
!
! Global Variables
!
      double precision  gw, gwwa, gwwz
      common /coup1/    gw, gwwa, gwwz
      double precision  gal(2),gau(2),gad(2),gwf(2)
      common /coup2a/   gal,   gau,   gad,   gwf
      double precision  gzn(2),gzl(2),gzu(2),gzd(2),g1(2)
      common /coup2b/   gzn,   gzl,   gzu,   gzd,   g1
      double precision  gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      common /coup3/    gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      complex*16        gchf(2,12)
      common /coup4/    gchf
      double precision  wmass1,wwidth1,zmass1,zwidth1
      common /vmass1/   wmass1,wwidth1,zmass1,zwidth1
      double precision  amass1,awidth1,hmass1,hwidth1
      common /vmass2/   amass1,awidth1,hmass1,hwidth1
      double precision  fmass(12), fwidth(12)
      common /fermions/ fmass,     fwidth
      double precision  gg(2), g
      common /coupqcd/  gg,    g

!-----------
! Begin Code
!-----------
c      write(6,*) 'umass',umass
c      write(6,*) 'bmass',bmass
c      write(6,*) 'cmass',cmass
c      write(6,*) 'tmass',tmass
c      write(6,*) 'emass',emass
c      write(6,*) 'gw',gw

c      sw2=1d0-wmass**2/zmass**2
      sw2=0.2312d0 

      fmass(1) = emass
      fmass(2) = 0d0
      fmass(3) = umass
      fmass(4) = dmass
      fmass(5) = mumass
      fmass(6) = 0d0
      fmass(7) = cmass
      fmass(8) = smass
      fmass(9) = taumass
      fmass(10)= 0d0
      fmass(11)= tmass
      fmass(12)= bmass

      fwidth(1) = ewidth
      fwidth(2) = 0d0
      fwidth(3) = uwidth
      fwidth(4) = dwidth
      fwidth(5) = muwidth
      fwidth(6) = 0d0
      fwidth(7) = cwidth
      fwidth(8) = swidth
      fwidth(9) = tauwidth
      fwidth(10)= 0d0
      fwidth(11)= twidth
      fwidth(12)= bwidth

      wmass1=wmass
      zmass1=zmass
      amass1=amass
      hmass1=hmass
      wwidth1=wwidth
      zwidth1=zwidth
      awidth1=awidth
      hwidth1=hwidth

!  Calls to Helas routines to set couplings

      call coup1x(sw2,gw,gwwa,gwwz)
      call coup2x(sw2,gal,gau,gad,gwf,gzn,gzl,gzu,gzd,g1)
      call coup3x(sw2,zmass,hmass,gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh)
      do i=1,12
         call coup4x(sw2,zmass,fmass(i),gchf(1,i))
      enddo

!  QCD couplings

      g = 1d0
      gg(1)=-g
      gg(2)=-g

      end
