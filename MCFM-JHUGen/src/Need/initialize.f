      subroutine initialize
      implicit none
      include 'types.f'
!******************************************************************
!     sets up masses and coupling constants for Helas
!******************************************************************
      
!
! Constants
!
      real(dp)::  sw2
c      parameter        (sw2 = 0.2312_dp)

!   masses and widths of fermions

      real(dp):: tmass,      bmass,     cmass
      parameter       (tmass=175._dp,bmass=0._dp, cmass=0._dp)
      real(dp):: smass,      umass,     dmass
      parameter       (smass=0._dp,  umass=0._dp, dmass=0._dp)
      real(dp):: twidth,     bwidth,    cwidth
      parameter       (twidth=0._dp, bwidth=0._dp,cwidth=0._dp)
c      parameter       (twidth=1.5525210_dp, bwidth=0._dp,cwidth=0._dp)
      real(dp):: swidth,     uwidth,    dwidth
      parameter       (swidth=0._dp, uwidth=0._dp,dwidth=0._dp)
      real(dp):: emass,      mumass,    taumass
      parameter       (emass=0._dp,  mumass=0._dp,taumass=0._dp)
      real(dp):: ewidth,     muwidth,    tauwidth
      parameter       (ewidth=0._dp, muwidth=0._dp,tauwidth=0._dp)

!   masses and widths of bosons

      real(dp):: wmass,         zmass
<<<<<<< .mine
      parameter       (wmass=80.398d0,    zmass=91.1876d0)
=======
      parameter       (wmass=80.419_dp,    zmass=91.188_dp)
>>>>>>> .r332
      real(dp):: wwidth,        zwidth
<<<<<<< .mine
c      parameter       (wwidth=2.06d0, zwidth=2.49d0)
      parameter       (wwidth=2.1054d0, zwidth=2.4952d0)
=======
c      parameter       (wwidth=2.06_dp, zwidth=2.49_dp)
      parameter       (wwidth=0._dp, zwidth=0._dp)
>>>>>>> .r332
      real(dp):: amass,         awidth
      parameter       (amass=0._dp,     awidth=0._dp)
      real(dp):: hmass,         hwidth
<<<<<<< .mine
      parameter       (hmass=125d0,   hwidth=2.4251446D-10)
=======
      parameter       (hmass=100._dp,   hwidth=2.4251446e-10_dp)
>>>>>>> .r332
!
! Local Variables
!
      integer:: i
!
! Global Variables
!
      real(dp)::  gw, gwwa, gwwz
      common /coup1/    gw, gwwa, gwwz
      real(dp)::  gal(2),gau(2),gad(2),gwf(2)
      common /coup2a/   gal,   gau,   gad,   gwf
      real(dp)::  gzn(2),gzl(2),gzu(2),gzd(2),g1(2)
      common /coup2b/   gzn,   gzl,   gzu,   gzd,   g1
      real(dp)::  gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      common /coup3/    gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      complex*16        gchf(2,12)
      common /coup4/    gchf
      real(dp)::  wmass1,wwidth1,zmass1,zwidth1
      common /vmass1/   wmass1,wwidth1,zmass1,zwidth1
      real(dp)::  amass1,awidth1,hmass1,hwidth1
      common /vmass2/   amass1,awidth1,hmass1,hwidth1
      real(dp)::  fmass(12), fwidth(12)
      common /fermions/ fmass,     fwidth
      real(dp)::  gg(2), g
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

c      sw2=1._dp-wmass**2/zmass**2
      sw2=0.2312_dp 

      fmass(1) = emass
      fmass(2) = 0._dp
      fmass(3) = umass
      fmass(4) = dmass
      fmass(5) = mumass
      fmass(6) = 0._dp
      fmass(7) = cmass
      fmass(8) = smass
      fmass(9) = taumass
      fmass(10)= 0._dp
      fmass(11)= tmass
      fmass(12)= bmass

      fwidth(1) = ewidth
      fwidth(2) = 0._dp
      fwidth(3) = uwidth
      fwidth(4) = dwidth
      fwidth(5) = muwidth
      fwidth(6) = 0._dp
      fwidth(7) = cwidth
      fwidth(8) = swidth
      fwidth(9) = tauwidth
      fwidth(10)= 0._dp
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

      g = 1._dp
      gg(1)=-g
      gg(2)=-g

      end
