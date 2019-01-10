      program mcfm
      implicit none
      include "omp_lib.h"

      integer t,threadmin,threadmax,trials,trial
      integer mth
      parameter (mth=240)
      double precision rone,r0,r(mth),er0,er(mth)
      double precision tone,t0,t1,time(mth)
      double precision trialmin(mth),trialmax(mth),trialavg(mth)
      double precision speedmin(mth),speedmax(mth),speedavg(mth)

      tone=1e0
      threadmax=omp_get_max_threads()
!set threadmin=threadmax to run at maximum number of threads only
      threadmin=1
      threadmin=threadmax
      if (threadmax.gt.mth) then
         write(*,*) 'increase parameter mth to at least ',threadmax
         stop
      endif
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) 'maximum number of threads available: ',threadmax
      
      do t=threadmin,threadmax
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) '>>>>>>>>',t,' threads <<<<<<<<'
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         write(*,*) 
         call omp_set_num_threads(t)
         t0=omp_get_wtime()
         call mcfmsub(r0,er0)
         t1=omp_get_wtime()
         if (t.eq.1) then 
            tone=t1-t0
            rone=r0
         endif
         time(t)=t1-t0
         r(t)=r0
         er(t)=er0
      enddo
c      write(*,*)
c      write(*,*)
c      write(*,*)
c      write(*,*)
      write(*,*)
      write(*,*) 'timing results:'
      do t=threadmin,threadmax
         if (threadmin .ne. threadmax) then
            write(*,1) t,time(t),tone/time(t),r(t)
         else
            write(*,2) t,time(t),r(t)
         endif
      enddo
      
      if (threadmin .ne. threadmax) then
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) "set 'threadmin=threadmax' in 'src/Need/mcfm_omp.f'"
      write(*,*) "for normal running with max number of threads only" 
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
      endif

 1    format('threads: ',i3,' time: ',f8.2,'  ratio: ',f8.2,
     .     ' x-section: ',e18.12,' fb')
 2    format('threads: ',i3,' time: ',f8.2,
     .     '    x-section: ',e18.12,' fb')

      
      
      end
