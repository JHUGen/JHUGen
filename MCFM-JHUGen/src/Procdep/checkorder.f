      subroutine checkorder(order)
      implicit none
      include 'types.f'
c--- checks the value of nproc and part against a list of processes that
c--- are calculable only at LO. If the calculation is not possible,
c--- writes an error message and aborts
      
      include 'frag.f'
      include 'kpart.f'
      include 'nproc.f'
      character*1 order

c--- special cases where there is no LO calculation
      if ((kpart.ne.kreal) .and. (order == 'R')) then 
        write(6,*)
        write(6,*)'This process can only be calculated with part=real,'
        write(6,*)'it is a subset of a NLO calculation only.'
      stop
      endif
      
c--- if we're calculating LO only, there's no problem      
      if (kpart==klord) return
      
c--- otherwise, we must be performing a NLO calculation, and this list of
c--- process numbers can't be calculated beyond LO 
      if (order == 'L') then
        write(6,*)
        write(6,*)'This process cannot be calculated beyond LO - please'
        write(6,*)'check the values of nproc and part then try again'
        stop
      endif

c--- check that fragmentation is not turned on for a non-fragmentation process
      if ((frag) .and. (order .ne. 'F')) then
        write(6,*)
        write(6,*) 'This process does not include photon fragmentation.'
        write(6,*) 'Please set frag=.false. in the input file.'
        stop
      endif
       
      return
      end
  
