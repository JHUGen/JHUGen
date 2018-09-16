      subroutine mcfm_getunweighted()
      implicit none
      include 'types.f'
c--- Routine to obtain a given number of unweighted events

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'eventbuffer.f'
      include 'maxwt.f'
      include 'vegas_common.f'
      include 'iterat.f'

      real(dp):: integ,integ_err

      write(6,*) 'Trying to unweight events ...'
      unweight = .true.
      nprn=-1 ! suppress VEGAS output

c--- Call integration loop again, this time unweighting :
 10   call mcfm_vegas(1,1,ncall2,.false.,integ,integ_err)

      write(6,*) 'Events generated so far: ',numstored
      if (numstored < nevtrequested) goto 10

      return
      end

