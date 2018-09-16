      subroutine bookupeqbins(string,binsize,xlow,xhigh)
      implicit none
      include 'types.f'
      
      character *(*) string
      real * 8 binsize,xlow,xhigh
      include 'pwhg_bookhist-new.h'
      real * 8 xx(maxbins+1)
      integer:: k
      xx(1)=xlow
      do k=2,maxbins+1
         xx(k)=xx(k-1)+binsize
         if(xx(k)-(xhigh-binsize/1e4)>0) goto 10
      enddo
      write(*,*) 'bookupeqbins: too many bins in hist ',string
      call exit(-1)
 10   continue
      if((xx(k)-xhigh)/binsize>1e-4) then
         write(*,*) 'upper limit incompatible with bin size'
         write(*,*) 'replacing ',xhigh,' with ',xx(k)
         write(*,*) ' in histogram ',string
      endif
      call bookup(string,k-1,xx)
      end


      subroutine bookup(string,n,x)
      implicit none
      include 'types.f'
c Books up a histogram characterized by the tag string <string>,
c with n bins. The array x(n+1) are the bins endpoints,
c x(i) is the low extreme of bin i.
      
      character *(*) string
      integer:: n
      real * 8 x(n+1)
      include 'pwhg_bookhist-new.h'
      integer:: j,k
      integer:: indexhist
      if(n>maxbins) then
         write(*,*) ' maximum number of bins=',maxbins
         write(*,*) ' asked for ',n
         call exit(-1)
      endif
c indexhist(string) returns the histogram index if a histogram
c with tag string was already booked, otherwise it books a new histogram,
c and returns minus the value of its index
      j=-indexhist(string)
      if(j<0) then
         write(*,*) 'Histogram ',string,' already booked'
         call exit(-1)
      endif
      do k=1,n+1
         xhistarr(k,j)=x(k)
      enddo
c y and err values go from 0 to n+1, 0 being the underflow and n+1
c the overflow.
      do k=0,n+1
         yhistarr(k,j)=0
         yhistarr1(k,j)=0
         errhistarr1(k,j)=0
         yhistarr2(k,j)=0
         errhistarr2(k,j)=0
         nhits(k,j)=0
      enddo
      nbins(j)=n
      ient1(j)=0
      end
      

      function indexhist(string)
      
      character * (*) string
      include 'pwhg_bookhist-new.h'
      integer:: indexhist
      integer:: j
      do j=1,jhist
         if(stringhist(j)==string) then
            indexhist=j
            return
         endif
      enddo
      if(jhist==nhist) then
         write(*,*) ' no more rooms for histograms'
         write(*,*) ' Histogram "',string,'" cannot be booked'
         call exit(-1)
      endif
      jhist=jhist+1
      stringhist(jhist)=string
      if(stringhist(jhist).ne.string) then
         write(*,*) ' Histogram string "',string,'" too long'
         call exit(-1)
      endif
c the negative sign indicates a new histogram
      indexhist=-jhist
      end

      subroutine filld(string,xval,weight)
      implicit none
      include 'types.f'
      
      character *(*) string
      real * 8 xval,weight
      include 'pwhg_bookhist-new.h'
      integer:: j,k,indexhist
      j=indexhist(string)
      if(j<0) then
         write(*,*) ' histogram "',string,'" was not booked'
         call exit(-1)
      endif
c underflow
      if(xval<xhistarr(1,j)) then
         yhistarr(0,j)=yhistarr(0,j)+weight
         nhits(0,j)=nhits(0,j)+1
         return
      else
         do k=1,nbins(j)
            if(xval<xhistarr(k+1,j)) then
               yhistarr(k,j)=yhistarr(k,j)+weight/
     1              (xhistarr(k+1,j)-xhistarr(k,j))
               nhits(k,j)=nhits(k,j)+1
               return
            endif
         enddo
      endif
c overflow
      yhistarr(nbins(j)+1,j)=yhistarr(nbins(j)+1,j)+weight
      end


      subroutine inihists
      implicit none
      include 'types.f'
      
      include 'pwhg_bookhist-new.h'
      jhist=0
      end

      subroutine pwhgtopout
      implicit none
      include 'types.f'
      
      include 'pwhg_bookhist-new.h'
      integer:: k,j,iun,l
      character * 50 outfile
      parameter (iun=99)
      do j=1,jhist

         write(iun,'(a,a,a,i3)')'# ',trim(stringhist(j)),'  index ',j-1
         do k=1,nbins(j)
            write(iun,'(4(1x,e14.8))') xhistarr(k,j), xhistarr(k+1,j),
     1           yhistarr2(k,j),errhistarr2(k,j)
         enddo
         write(iun,*)
         write(iun,*)
      enddo
      end


      subroutine pwhgaccumup
      implicit none
      include 'types.f'
c values histogrammed so far are transferred to array yhistarr1,
c and the square of the values are transferred to array errhistarr1.
c yhistarr is zeroed. The index ient1 is increased by one unit.
      
      include 'pwhg_bookhist-new.h'
      integer:: j,k
      do j=1,jhist
         do k=0,nbins(j)+1
            yhistarr1(k,j)=yhistarr1(k,j)+yhistarr(k,j)
            errhistarr1(k,j)=errhistarr1(k,j)+yhistarr(k,j)**2
            yhistarr(k,j)=0
         enddo
         ient1(j)=ient1(j)+1
      enddo
      end

      subroutine pwhgsetout
      implicit none
      include 'types.f'
c provides a snapshot of the current result of the
c analysis, leaving the yhistarr1 and errhistarr1 unchanged.
      
      include 'pwhg_bookhist-new.h'
      integer:: j,k
      real *8 xxx,sum,sumsq
      do j=1,jhist
         xxx=1._dp/ient1(j)
         do k=0,nbins(j)+1
            sum=yhistarr1(k,j)
            sumsq=errhistarr1(k,j)
            yhistarr2(k,j)=xxx*sum
            errhistarr2(k,j)=sqrt(xxx**2*abs(sumsq-sum**2*xxx))
         enddo
      enddo
      end

      subroutine pwhgaddout
      implicit none
      include 'types.f'
c accumulates the results obtained so far in yhistarr2 and errhistarr2.
c It zeroes yhistarr1 and errhistarr1. To be used if we compute
c a cross section with several contributions.
      
      include 'pwhg_bookhist-new.h'
      integer:: j,k
      real *8 xxx,sum,sumsq
      do j=1,jhist
         xxx=1._dp/ient1(j)
         do k=0,nbins(j)+1
            sum=yhistarr1(k,j)
            sumsq=errhistarr1(k,j)
            yhistarr2(k,j)=yhistarr2(k,j)+xxx*sum
            errhistarr2(k,j)=sqrt(errhistarr2(k,j)**2+
     1           xxx**2*abs(sumsq-sum**2*xxx))
         enddo
      enddo
      do j=1,jhist
         do k=0,nbins(j)+1
            yhistarr1(k,j)=0
            errhistarr1(k,j)=0
         enddo
         ient1(j)=0
      enddo
      end
