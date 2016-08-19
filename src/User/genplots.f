************************************************************************
*                                                                      *
*                        SINGLE MOMENTUM ip                            *
*                                                                      *
************************************************************************
      subroutine genplot1(p,ip,tag,wt,wt2,n)
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to the momentum p(ip)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      implicit none
      include 'constants.f'
      integer ip,n,zero
      double precision p(mxpart,4),wt,wt2,var,pt,yrap
      character*4 tag,titley
      character*5 titlept

      zero=ichar('0')

c--- particle pt 
      titlept='pt('//char(zero+ip)//')'    
      if (tag .eq. 'book') then
        var=0d0
      else
        var=pt(ip,p)
      endif
      call bookplot(n,tag,titlept,var,wt,wt2,
     &              0d0,200d0,5d0,'lin')
       n=n+1
      
c--- particle rapidity 
      titley='y('//char(zero+ip)//')'    
      if (tag .eq. 'book') then
        var=0d0
      else
        var=yrap(ip,p)
      endif
      call bookplot(n,tag,titley,var,wt,wt2,
     &              -5d0,5d0,0.4d0,'lin')
       n=n+1
      
      return
      end
 
            
************************************************************************
*                                                                      *
*                     COMPOUND MOMENTUM ip + jp                        *
*                                                                      *
************************************************************************
      subroutine genplot2(p,ip,jp,tag,wt,wt2,n)
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to the momentum p(ip)+p(jp)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      implicit none
      include 'constants.f'
      integer ip,jp,n,zero
      double precision p(mxpart,4),wt,wt2,var,pttwo,yraptwo
      character*4 tag
      character*6 titley,titlem
      character*7 titlept

      zero=ichar('0')

c--- two-particle invariant mass 
      titlem='m('//char(zero+ip)//','//char(zero+jp)//')'    
      if (tag .eq. 'book') then
        var=0d0
      else
        var=dsqrt((p(ip,4)+p(jp,4))**2-(p(ip,1)+p(jp,1))**2
     &           -(p(ip,2)+p(jp,2))**2-(p(ip,3)+p(jp,3))**2)
      endif
      call bookplot(n,tag,titlem,var,wt,wt2,
     &              0d0,300d0,5d0,'lin')
       n=n+1

c--- two-particle pt 
      titlept='pt('//char(zero+ip)//','//char(zero+jp)//')'    
      if (tag .eq. 'book') then
        var=0d0
      else
        var=pttwo(ip,jp,p)
      endif
      call bookplot(n,tag,titlept,var,wt,wt2,
     &              0d0,200d0,5d0,'lin')
       n=n+1
      
c--- two-particle rapidity 
      titley='y('//char(zero+ip)//','//char(zero+jp)//')'    
      if (tag .eq. 'book') then
        var=0d0
      else
        var=yraptwo(ip,jp,p)
      endif
      call bookplot(n,tag,titley,var,wt,wt2,
     &              -5d0,5d0,0.4d0,'lin')
       n=n+1
      
      return
      end
            
            
************************************************************************
*                                                                      *
*                   COMPOUND MOMENTUM ip + jp + kp                     *
*                                                                      *
************************************************************************
      subroutine genplot3(p,ip,jp,kp,tag,wt,wt2,n)
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to the momentum p(ip)+p(jp)+p(kp)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      implicit none
      include 'constants.f'
      integer ip,jp,kp,n,zero
      double precision p(mxpart,4),wt,wt2,var,ptthree,yrapthree
      character*4 tag
      character*8 titley,titlem
      character*9 titlept

      zero=ichar('0')

c--- three-particle invariant mass 
      titlem='m('//char(zero+ip)//','//char(zero+jp)//
     &                            ','//char(zero+kp)//')'
      if (tag .eq. 'book') then
        var=0d0
      else
        var=dsqrt((p(ip,4)+p(jp,4)+p(kp,4))**2
     &           -(p(ip,1)+p(jp,1)+p(kp,1))**2
     &           -(p(ip,2)+p(jp,2)+p(kp,2))**2
     &           -(p(ip,3)+p(jp,3)+p(kp,3))**2)
      endif
      call bookplot(n,tag,titlem,var,wt,wt2,
     &              0d0,400d0,5d0,'lin')
       n=n+1
      
c--- three-particle pt 
      titlept='pt('//char(zero+ip)//','//char(zero+jp)//
     &                              ','//char(zero+kp)//')'   
      if (tag .eq. 'book') then
        var=0d0
      else
        var=ptthree(ip,jp,kp,p)
      endif
      call bookplot(n,tag,titlept,var,wt,wt2,
     &              0d0,200d0,5d0,'lin')
       n=n+1
      
c--- three-particle rapidity 
      titley='y('//char(zero+ip)//','//char(zero+jp)//
     &                            ','//char(zero+kp)//')'
      if (tag .eq. 'book') then
        var=0d0
      else
        var=yrapthree(ip,jp,kp,p)
      endif
      call bookplot(n,tag,titley,var,wt,wt2,
     &              -5d0,5d0,0.4d0,'lin')
       n=n+1
      
      return
      end
            
