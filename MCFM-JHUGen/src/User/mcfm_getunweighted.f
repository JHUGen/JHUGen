      subroutine mcfm_getunweighted()
c--- Routine to obtain a given number of unweighted events
      implicit none
      include 'constants.f'
      include 'eventbuffer.f'
      include 'maxwt.f'
      include 'vegas_common.f'
      double precision integ,integ_err
      integer itmx1,ncall1,itmx2,ncall2
      common/iterat/itmx1,ncall1,itmx2,ncall2

      write(6,*) 'Trying to unweight events ...'
      unweight = .true.
      nprn=-1 ! suppress VEGAS output

c--- Call integration loop again, this time unweighting :
 10   call mcfm_vegas(1,1,ncall2,.false.,integ,integ_err)
 
      write(6,*) 'Events generated so far: ',numstored
      if (numstored .lt. nevtrequested) goto 10

      return
      end
      
