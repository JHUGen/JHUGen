      subroutine hisini()
      implicit none
      integer iplot,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer binnum(maxhis)
      double precision binmin(maxhis),binmax(maxhis) 
      double precision hbin(maxhis,maxbin),herr(maxhis,maxbin)
      integer hnum(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      iplot=0
      hbin(:,:)=0d0
      herr(:,:)=0d0
      hnum(:,:)=0

      return
      end

      subroutine hisres()
      implicit none
      integer i,iplot,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer binnum(maxhis)
      integer hnum(maxhis,maxbin)
      double precision binmin(maxhis),binmax(maxhis) 
      double precision hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      do i=1,iplot
         hbin(i,:)=0d0
         herr(i,:)=0d0
         hnum(i,:)=0d0
      enddo

      return
      end

      subroutine hisset(nbin,bmin,bmax,info)
      implicit none
! this is not in the parallel region and can be as complex as you like
      integer nbin,info
      integer i,nhist
      integer iplot,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer binnum(maxhis)
      integer hnum(maxhis,maxbin)
      double precision bmin,bmax
      double precision binmin(maxhis),binmax(maxhis) 
      double precision hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      iplot=iplot+1
      if (iplot.le.maxhis) then
         binnum(iplot)=min(nbin,maxbin)
         binmin(iplot)=bmin
         binmax(iplot)=bmax
      endif

      return
      end

      subroutine hisout(iter)
      implicit none
! this is not in the parallel region and can be as complex as you like
!      include 'histo.f'
!      include 'mcfmplotinfo.f'
      integer i,j,iplot,iter,maxhis,maxbin,num
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer binnum(maxhis)
      integer hnum(maxhis,maxbin)
      double precision val,bin,err,sd
      double precision binmin(maxhis),binmax(maxhis) 
      double precision hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      do i=1,iplot
         write(*,*) 'histogram',i
         do j=1,binnum(i)
            val=binmin(i)+j*(binmax(i)-binmin(i))/binnum(i)
            bin=hbin(i,j)
            err=herr(i,j)
            num=hnum(i,j)
            sd=0d0
            if (num.gt.0) sd=sqrt(err-bin**2/num)
            write(*,*) val,bin/iter,sd/iter,num
         enddo
      enddo
      
      return
      end

      subroutine hisbin(p,wt)
      implicit none
! this is in the parallel region
!      include 'histo.f'
!      include 'mcfmplotinfo.f'
      include 'mxpart.f'
      integer i,iplot,binval,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer binnum(maxhis)
      integer hnum(maxhis,maxbin)
      double precision wt,wt2,val,hisobs
      double precision p(mxpart,4)
      double precision binmin(maxhis),binmax(maxhis) 
      double precision hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      wt2=wt**2
      do i=1,iplot
         val=hisobs(i,p)
         val=(val-binmin(i))/(binmax(i)-binmin(i))
         binval=1+int((binnum(i)-1)*val)
         if (binval.le.binnum(i)) then
!$omp critical(FastHis)
            hbin(i,binval)=hbin(i,binval)+wt
            herr(i,binval)=herr(i,binval)+wt2
            hnum(i,binval)=hnum(i,binval)+1
!$omp end critical(FastHis)
         endif
      enddo
      
      return
      end

      double precision function hisobs(ihist,p)
      implicit none
! Returns the value of the observable for histogram ihist.
! This routine lives in the parallel region
      include 'mxpart.f'
      include 'mcfmplotinfo.f'
      include 'process.f'
      include 'nproc.f'
      integer ihist,ilomomenta
      double precision p(mxpart,4)
      
      hisobs=0d0
      if (nproc.eq.1) then
         if (ihist.eq.1) then
            hisobs=sqrt(p(3,1)**2+p(3,2)**2)
         endif
      endif

      return
      end
