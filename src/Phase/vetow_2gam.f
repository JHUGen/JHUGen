      logical function vetow_2gam(p)
      implicit none
c--- returns TRUE if the momenta passed into the function
c--- should be vetoed according to the current value of "ipsgen"
      include 'constants.f'
      include 'ipsgen.f'
      include 'masses.f'
      include 'part.f'
      double precision p(mxpart,4),dot,s345,s346,s3456,xwid
c--- note: parameter "xwid" controls the number of widths away
c--- from the peak to generate according to a BW
      parameter(xwid=5d0)
      
      vetow_2gam=.false.
      
c--- veto for fragmentation contribution
      if (part .eq. 'frag') then
        if (ipsgen .eq. 1) then
          s345=2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5))
          if (abs(sqrt(s345)-wmass) .lt. xwid*wwidth) then
            vetow_2gam=.true.
          endif
        endif
        return
      endif

c--- veto for everything else     
      if (ipsgen .eq. 1) then
        s345=2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5))
        if (abs(sqrt(s345)-wmass) .lt. xwid*wwidth) then
          vetow_2gam=.true.
          return
        endif
      endif
      if ((ipsgen .eq. 1) .or. (ipsgen .eq. 3)) then
        s346=2d0*(dot(p,3,4)+dot(p,3,6)+dot(p,4,6))
        if (abs(sqrt(s346)-wmass) .lt. xwid*wwidth) then
          vetow_2gam=.true.
          return
        endif
      endif
      if ((ipsgen .eq. 1) .or. (ipsgen .eq. 3) .or. (ipsgen .eq. 4))then
        s3456=2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &            +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
        if (abs(sqrt(s3456)-wmass) .lt. xwid*wwidth) then
          vetow_2gam=.true.
          return
        endif
      endif
      
      return
      end
      
