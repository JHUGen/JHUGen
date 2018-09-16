      subroutine hisini()
      implicit none
      include 'types.f'
      
      integer:: iplot,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer:: binnum(maxhis)
      real(dp):: binmin(maxhis),binmax(maxhis) 
      real(dp):: hbin(maxhis,maxbin),herr(maxhis,maxbin)
      integer:: hnum(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      iplot=0
      hbin(:,:)=0._dp
      herr(:,:)=0._dp
      hnum(:,:)=0

      return
      end

      subroutine hisres()
      implicit none
      include 'types.f'
      
      integer:: i,iplot,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer:: binnum(maxhis)
      integer:: hnum(maxhis,maxbin)
      real(dp):: binmin(maxhis),binmax(maxhis) 
      real(dp):: hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      do i=1,iplot
         hbin(i,:)=0._dp
         herr(i,:)=0._dp
         hnum(i,:)=0._dp
      enddo

      return
      end

      subroutine hisset(nbin,bmin,bmax,info)
      implicit none
      include 'types.f'
      
! this is not in the parallel region and can be as complex as you like
      integer:: nbin,info
!      integer:: i,nhist
      integer:: iplot,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer:: binnum(maxhis)
      integer:: hnum(maxhis,maxbin)
      real(dp):: bmin,bmax
      real(dp):: binmin(maxhis),binmax(maxhis) 
      real(dp):: hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      iplot=iplot+1
      if (iplot<=maxhis) then
         binnum(iplot)=min(nbin,maxbin)
         binmin(iplot)=bmin
         binmax(iplot)=bmax
      endif

      return
      end

      subroutine hisout(iter)
      implicit none
      include 'types.f'
      
! this is not in the parallel region and can be as complex as you like
!      include 'histo.f'
!      include 'mcfmplotinfo.f'
      integer:: i,j,iplot,iter,maxhis,maxbin,num
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer:: binnum(maxhis)
      integer:: hnum(maxhis,maxbin)
      real(dp):: val,bin,err,sd
      real(dp):: binmin(maxhis),binmax(maxhis) 
      real(dp):: hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      do i=1,iplot
         write(*,*) 'histogram',i
         do j=1,binnum(i)
            val=binmin(i)+j*(binmax(i)-binmin(i))/binnum(i)
            bin=hbin(i,j)
            err=herr(i,j)
            num=hnum(i,j)
            sd=0._dp
            if (num>0) sd=sqrt(err-bin**2/num)
            write(*,*) val,bin/iter,sd/iter,num
         enddo
      enddo
      
      return
      end

      subroutine hisbin(p,wt)
      implicit none
      include 'types.f'
      
! this is in the parallel region
!      include 'histo.f'
!      include 'mcfmplotinfo.f'
      include 'mxpart.f'
      integer:: i,iplot,binval,maxhis,maxbin
      parameter(maxhis=40)
      parameter(maxbin=100)
      integer:: binnum(maxhis)
      integer:: hnum(maxhis,maxbin)
      real(dp):: wt,wt2,val,hisobs
      real(dp):: p(mxpart,4)
      real(dp):: binmin(maxhis),binmax(maxhis) 
      real(dp):: hbin(maxhis,maxbin),herr(maxhis,maxbin)
      common/hisfas/hbin,herr,hnum,binnum,binmin,binmax,iplot

      wt2=wt**2
      do i=1,iplot
         val=hisobs(i,p)
         val=(val-binmin(i))/(binmax(i)-binmin(i))
         binval=1+int((binnum(i)-1)*val)
         if (binval<=binnum(i)) then
!$omp critical(FastHis)
            hbin(i,binval)=hbin(i,binval)+wt
            herr(i,binval)=herr(i,binval)+wt2
            hnum(i,binval)=hnum(i,binval)+1
!$omp end critical(FastHis)
         endif
      enddo
      
      return
      end

      function hisobs(ihist,p)
      implicit none
      include 'types.f'
      real(dp):: hisobs
      
! Returns the value of the observable for histogram ihist.
! This routine lives in the parallel region
      include 'mxpart.f'
      include 'mcfmplotinfo.f'
      include 'kprocess.f'
      include 'nproc.f'
      integer:: ihist
      real(dp):: p(mxpart,4)
      
      hisobs=0._dp
      if (nproc==1) then
         if (ihist==1) then
            hisobs=sqrt(p(3,1)**2+p(3,2)**2)
         endif
      endif

      return
      end
