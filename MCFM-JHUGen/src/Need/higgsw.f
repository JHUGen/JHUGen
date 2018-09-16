      subroutine higgsw(br)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      real(dp):: wff,mfsq,br
c-----approximate form for the width of the standard model higgs
c-----valid for low masses
      wff(mfsq)=sqrt(2._dp)/8._dp/pi*gf*hmass*mfsq
     & *(1._dp-4._dp*mfsq/hmass**2)**1.5_dp

      hwidth=3._dp*(wff(mbsq)+wff(mcsq))+wff(mtausq)
      write(6,*) 'hmass,hwidth',hmass,hwidth
      write(6,*) 'mtausq,mcsq,mbsq',mtausq,mcsq,mbsq
      write(6,*)
      br=3._dp*wff(mbsq)/hwidth
      return
      end
