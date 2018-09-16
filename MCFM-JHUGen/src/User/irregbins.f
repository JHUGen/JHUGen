c--- Set of routines for producing histograms with irregular bins
      subroutine initirregbins(n,binedges)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'histo.f'
      include 'irregbins_incl.f'
      integer:: n,j
      real(dp):: binedges(30)
      logical:: first
      data first/.true./
      save first

      if (first) then
        nirreg=0
        first=.false.
      endif

      nirreg=nirreg+1

      if (nirreg > 10) then
        write(6,*) 'Maximum no. of irregular histograms exceeded!'
        stop
      endif

      irregbin(n)=.true.
      irregbin_ptr(n)=nirreg

      do j=1,30
      irregbinedges(nirreg,j)=binedges(j)
      enddo

      return
      end



      subroutine getirregbins(n)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'histo.f'
      include 'irregbins_incl.f'
      integer:: n,ib,j,k,nib
      real(dp):: yy(30),ee(30),del(30),ytot,ylo,yhi
c--- added these variables to scale plots at intermediate steps
      logical:: scaleplots
      real(dp):: scalefac
      common/scaleplots/scalefac,scaleplots

c--- n is the histogram number in the usual output
c--- ib is the index into irregbinedges
      ib=irregbin_ptr(n)

c--- determine number of irregular bins (nib)
      nib=0
      do j=1,30
      if (irregbinedges(ib,j) > 1.e-10_dp) nib=j
      enddo
      nib=nib-1 ! the above over-counts by 1

c--- j labels irregular bins
      do j=1,nib
      yy(j)=zip
      ee(j)=zip
      del(j)=irregbinedges(ib,j+1)-irregbinedges(ib,j)

c--- k labels regular histogram bins
      do k=1,NBIN(n)
        if ( (XHIS(n,k) >= irregbinedges(ib,j)) .and.
     &       (XHIS(n,k) <= irregbinedges(ib,j+1)) ) then
          yy(j)=yy(j)+HIST(n,k)*HDEL(n)/del(j)
          ee(j)=ee(j)+(HIST(2*maxhisto+n,k)*HDEL(n)/del(j))**2
        endif
      enddo
      ee(j)=sqrt(max(ee(j),zip))

      enddo

c--- under- and over-flow
      ylo=zip
      yhi=zip
      do k=1,NBIN(n)
        if (XHIS(n,k) < irregbinedges(ib,1)) then
          ylo=ylo+HIST(n,k)*HDEL(n)
        endif
        if (XHIS(n,k) > irregbinedges(ib,nib+1)) then
          yhi=yhi+HIST(n,k)*HDEL(n)
        endif
      enddo

c--- rescaling at intermediate stages
      if (scaleplots) then
        yy(:)=yy(:)*scalefac
        ylo=ylo*scalefac
        yhi=yhi*scalefac
      endif

      write(6,*)
      write(6,*) 'Histogram: ',trim(TITLE(n))
      write(6,*) '(results in fb/GeV, with integration error)'
      ytot=zip
      do j=1,nib
      write(6,99) irregbinedges(ib,j),irregbinedges(ib,j+1),yy(j),ee(j)
      ytot=ytot+yy(j)*del(j)
      enddo

      write(6,*) ' Cross sections in fb'
      write(6,*) '  ---> total of all bins:',ytot
      write(6,*) '  ---> underflow        :',ylo
      write(6,*) '  ---> overflow         :',yhi
      write(6,*) '  ---> total u/ & o/flow:',ytot+ylo+yhi

      return

   99 format('[',f9.3,' , ',f9.3,']',g15.6,'  +/- ',g15.6)

      end

