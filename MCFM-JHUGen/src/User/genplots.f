************************************************************************
*                                                                      *
*                        SINGLE MOMENTUM ip                            *
*                                                                      *
************************************************************************
      subroutine genplot1(p,ip,tag,wt,wt2,n)
      implicit none
      include 'types.f'
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to the momentum p(ip)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: ip,n,izero
      real(dp):: p(mxpart,4),wt,wt2,var,pt,yrap
      integer tag
      character*4 titley
      character*5 titlept
      integer, parameter:: tagbook=1, tagplot=2

      izero=ichar('0')

c--- particle pt 
      titlept='pt('//char(izero+ip)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=pt(ip,p)
      endif
      call bookplot(n,tag,titlept,var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
       n=n+1
      
c--- particle rapidity 
      titley='y('//char(izero+ip)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=yrap(ip,p)
      endif
      call bookplot(n,tag,titley,var,wt,wt2,
     &              -5._dp,5._dp,0.4_dp,'lin')
       n=n+1
      
      return
      end
 
            
************************************************************************
*                                                                      *
*                     COMPOUND MOMENTUM ip + jp                        *
*                                                                      *
************************************************************************
      subroutine genplot2(p,ip,jp,tag,wt,wt2,n)
      implicit none
      include 'types.f'
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to the momentum p(ip)+p(jp)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: ip,jp,n,izero
      real(dp):: p(mxpart,4),wt,wt2,var,pttwo,yraptwo
      integer tag
      character*6 titley,titlem
      character*7 titlept
      integer, parameter:: tagbook=1, tagplot=2

      izero=ichar('0')

c--- two-particle invariant mass 
      titlem='m('//char(izero+ip)//','//char(izero+jp)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=sqrt((p(ip,4)+p(jp,4))**2-(p(ip,1)+p(jp,1))**2
     &           -(p(ip,2)+p(jp,2))**2-(p(ip,3)+p(jp,3))**2)
      endif
      call bookplot(n,tag,titlem,var,wt,wt2,
     &              0._dp,300._dp,5._dp,'lin')
       n=n+1

c--- two-particle pt 
      titlept='pt('//char(izero+ip)//','//char(izero+jp)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=pttwo(ip,jp,p)
      endif
      call bookplot(n,tag,titlept,var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
       n=n+1
      
c--- two-particle rapidity 
      titley='y('//char(izero+ip)//','//char(izero+jp)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=yraptwo(ip,jp,p)
      endif
      call bookplot(n,tag,titley,var,wt,wt2,
     &              -5._dp,5._dp,0.4_dp,'lin')
       n=n+1
      
      return
      end
            
            
************************************************************************
*                                                                      *
*                   COMPOUND MOMENTUM ip + jp + kp                     *
*                                                                      *
************************************************************************
      subroutine genplot3(p,ip,jp,kp,tag,wt,wt2,n)
      implicit none
      include 'types.f'
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to the momentum p(ip)+p(jp)+p(kp)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: ip,jp,kp,n,izero
      real(dp):: p(mxpart,4),wt,wt2,var,ptthree,yrapthree
      integer tag
      character*8 titley,titlem
      character*9 titlept
      integer, parameter:: tagbook=1, tagplot=2

      izero=ichar('0')

c--- three-particle invariant mass 
      titlem='m('//char(izero+ip)//','//char(izero+jp)//
     &                            ','//char(izero+kp)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=sqrt((p(ip,4)+p(jp,4)+p(kp,4))**2
     &           -(p(ip,1)+p(jp,1)+p(kp,1))**2
     &           -(p(ip,2)+p(jp,2)+p(kp,2))**2
     &           -(p(ip,3)+p(jp,3)+p(kp,3))**2)
      endif
      call bookplot(n,tag,titlem,var,wt,wt2,
     &              0._dp,400._dp,5._dp,'lin')
       n=n+1
      
c--- three-particle pt 
      titlept='pt('//char(izero+ip)//','//char(izero+jp)//
     &                              ','//char(izero+kp)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=ptthree(ip,jp,kp,p)
      endif
      call bookplot(n,tag,titlept,var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
       n=n+1
      
c--- three-particle rapidity 
      titley='y('//char(izero+ip)//','//char(izero+jp)//
     &                            ','//char(izero+kp)//')'
      if (tag == tagbook) then
        var=0._dp
      else
        var=yrapthree(ip,jp,kp,p)
      endif
      call bookplot(n,tag,titley,var,wt,wt2,
     &              -5._dp,5._dp,0.4_dp,'lin')
       n=n+1
      
      return
      end
            
