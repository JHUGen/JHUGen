      subroutine smalls(s,npart,*)
c    cut if radiated parton too close
      implicit none
      include 'constants.f'
      include 'cutoff.f'
      include 'process.f'
      integer npart
      double precision s(mxpart,mxpart)

      if (s(1,2) .lt. cutoff) return 1

      if (npart .eq. 2) then
      if ( 
     .      (-s(1,4) .lt. cutoff)
     . .or. (-s(2,4) .lt. cutoff)
     . .or. (+s(3,4) .lt. cutoff)
     . .or. (-s(1,3) .lt. cutoff)
     . .or. (-s(2,3) .lt. cutoff)
     . ) return 1
      
      elseif (npart .eq. 3) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,4) .lt. cutoff)
     . .or. (-s(2,4) .lt. cutoff)
     . .or. (+s(4,5) .lt. cutoff)
     . .or. (+s(3,4) .lt. cutoff)
     . .or. (+s(3,5) .lt. cutoff)
     . .or. (-s(1,3) .lt. cutoff)
     . .or. (-s(2,3) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 4) then
      if ( 
     .      (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,4) .lt. cutoff)
     . .or. (-s(2,4) .lt. cutoff)
     . .or. (-s(1,3) .lt. cutoff)
     . .or. (-s(2,3) .lt. cutoff)
     . .or. (+s(3,4) .lt. cutoff)
     . .or. (+s(3,5) .lt. cutoff)
     . .or. (+s(3,6) .lt. cutoff)
     . .or. (+s(4,5) .lt. cutoff)
     . .or. (+s(4,6) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . ) return 1
c      if (  (case .eq. 'qq_tbg') .or. (case .eq. 'qqtbgg')
c     . .or. (case .eq. 'epem3j') .or. (case .eq. 'W_tndk')
c     . .or. (case .eq. 'Z_tjet')) then
c      if ( 
c     .      (+s(3,4) .lt. cutoff)
c     . .or. (+s(3,5) .lt. cutoff)
c     . .or. (+s(3,6) .lt. cutoff)
c     . .or. (+s(4,5) .lt. cutoff)
c     . .or. (+s(4,6) .lt. cutoff)
c     . ) return 1
c      endif
     
      elseif (npart .eq. 5) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . .or. (-s(1,4) .lt. cutoff)
     . .or. (-s(2,4) .lt. cutoff)
     . .or. (-s(1,3) .lt. cutoff)
     . .or. (-s(2,3) .lt. cutoff)
     . .or. (+s(3,4) .lt. cutoff)
     . .or. (+s(3,5) .lt. cutoff)
     . .or. (+s(3,6) .lt. cutoff)
     . .or. (+s(3,7) .lt. cutoff)
     . .or. (+s(4,5) .lt. cutoff)
     . .or. (+s(4,6) .lt. cutoff)
     . .or. (+s(4,7) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . .or. (+s(5,7) .lt. cutoff)
     . .or. (+s(6,7) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 6) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . .or. (-s(1,8) .lt. cutoff)
     . .or. (-s(2,8) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . .or. (+s(5,7) .lt. cutoff)
     . .or. (+s(5,8) .lt. cutoff)
     . .or. (+s(6,7) .lt. cutoff)
     . .or. (+s(6,8) .lt. cutoff)
     . .or. (+s(7,8) .lt. cutoff)
     . ) return 1
      
      elseif (npart .eq. 7) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . .or. (-s(1,8) .lt. cutoff)
     . .or. (-s(2,8) .lt. cutoff)
     . .or. (-s(1,9) .lt. cutoff)
     . .or. (-s(2,9) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . .or. (+s(5,7) .lt. cutoff)
     . .or. (+s(5,8) .lt. cutoff)
     . .or. (+s(5,9) .lt. cutoff)
     . .or. (+s(6,7) .lt. cutoff)
     . .or. (+s(6,8) .lt. cutoff)
     . .or. (+s(6,9) .lt. cutoff)
     . .or. (+s(7,8) .lt. cutoff)
     . .or. (+s(7,9) .lt. cutoff)
     . .or. (+s(8,9) .lt. cutoff)
     . ) return 1
      
      endif      
      
      return
      end
