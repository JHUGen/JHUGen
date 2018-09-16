      subroutine scaleset_mVpmH(p,mu0)
      implicit none
      include 'types.f'
!==== this is the scale setting routine for
!===  mV+mH, its not really a dynmaic scale, but its
!==== more convienent than carrying around junk inthe input.dAT
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'breit.f' 
      include 'cplx.h'
      real(dp):: p(mxpart,4),mu0

      if((kcase==kWHbbar) .or.
     &   (kcase==kWH1jet) .or.
     &   (kcase==kWHgaga) .or.
     &   (kcase==kWH__WW) .or.
     &   (kcase==kZH__WW) .or.
     &   (kcase==kZHgaga) .or.
     &   (kcase==kZHbbar) .or.
     &   (kcase==kZH1jet)) then 
         mu0=mass2+mass3
      else
        write(6,*)'dynamicscale mV+mH not supported for this process.'
        stop
      endif
      
      return
      end
      
