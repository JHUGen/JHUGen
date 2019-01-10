      subroutine interpolate_hto(mh,hwidth)
c--- This routine reads in the file hto_output.dat, prepared
c--- by the HTO program of G. Passarino, and uses linear interpolation
c--- to obtain the Higgs width parameter CPGH for any input Higgs
c--- mass (mh) between 50 GeV and 1500 GeV
      implicit none
      logical first
      double precision mh,hwidth,m1,g1,m2,g2,eps
      data first/.true./
      save first
      parameter(eps=1d-8)
      
      if ((mh .lt. 50d0) .or. (mh .gt. 1500d0)) then
        write(6,*) 'Higgs mass outside HTO interpolation range;'
        write(6,*) 'requires  50 GeV < mH < 1500 GeV'  
        stop
      endif

      if (first) then
        first=.false.
        write(6,*)
        write(6,*) 'Higgs width parameter obtained by interpolating:'
        write(6,*) '   '
        write(6,*) '*'
        write(6,*) ' *'
        write(6,*) '  * --- cpHTO v 1.1 (May 2012) by Giampiero '
        write(6,*) ' *'
        write(6,*) '*'      
      endif
      
      open(unit=12,file='hto_output.dat',status='unknown')
      read(12,*)
      read(12,*)
   77 continue
      read(12,99) m1,g1
      if (mh-m1 .gt. 0.5d0+eps) goto 77
      read(12,99) m2,g2
      close(12)
      
      hwidth=(g2-g1)*(mh-m1)/(m2-m1)+g1

      return

   99 format(f10.1,f21.12)
      end
      
