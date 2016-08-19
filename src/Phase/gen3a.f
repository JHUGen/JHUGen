      subroutine gen3a(r,p3,wt,*)
C---modified phase space generator, generating 2-2 and then 3 from 2
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'

      double precision p3(mxpart,4),p2(mxpart,4),r(mxdim),wt,wt2
c----although wt2 is generated we will not use it
      call gen2(r,p2,wt2,*999)      
c      write(6,*) 'wt2 in gen3a',wt2

c----this generates the full weight from both branchings
      call gen3from2(p2,r(5),r(6),r(7),p3,wt,*999)      
c      write(6,*) 
c      pause
      if (debug) then
      write(6,*) 'wt in gen3a',wt
      write(6,*) 'end of gen3a'
      endif

      return
 999  wt=0d0
      return 1
      end
