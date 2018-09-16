************************************************************************
*                                                                      *
*                        SINGLE MOMENTUM ip                            *
*                                                                      *
************************************************************************
      subroutine autoplot1(p,ip,tag,wt,wt2,n)
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
      include 'mcfmplotinfo.f'
      integer:: ip,n
      real(dp):: p(mxpart,4),wt,wt2,var,pt,yrap
      integer tag
      character*12 titley(12),titlept(12)
      character*6 cmom
      logical, save :: first(12)=.true.
      integer, parameter:: tagbook=1, tagplot=2
      save titley,titlept,cmom
!$omp threadprivate(first,titley,titlept,cmom)
      
c--- on first call for each plot, set up strings
      if (first(ip)) then
c--- convert integer:: to left-justified character string
        write(cmom,'(i6)') ip
        cmom=adjustl(cmom)
        titlept(ip)='pt('//trim(cmom)//')'
        titley(ip)='y('//trim(cmom)//')'
        first(ip)=.false.
      endif
      
c--- particle pt 
      if (tag == tagbook) then
        var=0._dp
      else
        var=pt(ip,p)
      endif
      call bookplot(n,tag,titlept(ip),var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
       n=n+1
      
c--- particle rapidity 
      if (tag == tagbook) then
        var=0._dp
      else
        var=yrap(ip,p)
      endif
      call bookplot(n,tag,titley(ip),var,wt,wt2,
     &              -2.4_dp,2.4_dp,0.24_dp,'lin')
       n=n+1
      
      return
      end
 
            
************************************************************************
*                                                                      *
*                     COMPOUND MOMENTUM ip + jp                        *
*                                                                      *
************************************************************************
      subroutine autoplot2(p,mm,ip,jp,tag,wt,wt2,n)
      implicit none
      include 'types.f'
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to momentum p(ip)+p(jp) == p(mm)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      
      include 'mcfmplotinfo.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: mm,ip,jp,n
      real(dp):: p(mxpart,4),wt,wt2,var,pttwo,yraptwo
      integer tag
      character*12 titley(10,10),titlept(10,10),titlem(10,10),
     & titleay(10,10)
      character*6 cmom
      logical, save :: first(10,10)=.true.
      integer, parameter:: tagbook=1, tagplot=2
      save titley,titlept,titlem,titleay,cmom
!$omp threadprivate(first,titley,titlept,titlem,titleay,cmom)
c--- on first call for each plot, set up strings      
      if (first(ip,jp)) then
c--- convert integer:: to left-justified character string
        write(cmom,'(i6)') mm
        if (mm == 90) write(cmom,'(i6)') 910 ! special code: 0 -> 10
        cmom=adjustl(cmom)
        titlem(ip,jp)='m('//trim(cmom)//')'    
        titlept(ip,jp)='pt('//trim(cmom)//')'    
        titley(ip,jp)='y('//trim(cmom)//')'    
        titleay(ip,jp)='|y('//trim(cmom)//')|'    
        first(ip,jp)=.false.
      endif
      
c--- two-particle invariant mass 
      if (tag == tagbook) then
        var=0._dp
      else
        var=sqrt((p(ip,4)+p(jp,4))**2-(p(ip,1)+p(jp,1))**2
     &           -(p(ip,2)+p(jp,2))**2-(p(ip,3)+p(jp,3))**2)
      endif
      call bookplot(n,tag,titlem(ip,jp),var,wt,wt2,
     &              0._dp,300._dp,5._dp,'lin')
      n=n+1

c--- two-particle pt 
      if (tag == tagbook) then
        var=0._dp
      else
        var=pttwo(ip,jp,p)
      endif
      call bookplot(n,tag,titlept(ip,jp),var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
      n=n+1
      
c--- two-particle rapidity and absolute rapidity 
      if (tag == tagbook) then
        var=0._dp
      else
        var=yraptwo(ip,jp,p)
      endif
      call bookplot(n,tag,titley(ip,jp),var,wt,wt2,
     &              -5._dp,5._dp,0.4_dp,'lin')
      n=n+1
      var=abs(var)
      call bookplot(n,tag,titleay(ip,jp),var,wt,wt2,
     &              0._dp,5._dp,0.25_dp,'lin')
      n=n+1
      
      return
      end
            
            
************************************************************************
*                                                                      *
*                   COMPOUND MOMENTUM ip + jp + kp                     *
*                                                                      *
************************************************************************
      subroutine autoplot3(p,mm,ip,jp,kp,tag,wt,wt2,n)
      implicit none
      include 'types.f'
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresponding to momentum p(ip)+p(jp)+p(kp) == p(mm)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mcfmplotinfo.f'
      integer:: mm,ip,jp,kp,n
      real(dp):: p(mxpart,4),wt,wt2,var,ptthree,yrapthree
      integer tag
      character*12 titley(10,10,10),titlept(10,10,10),titlem(10,10,10)
      character*6 cmom
      logical, save :: first(10,10,10)=.true.
      integer, parameter:: tagbook=1, tagplot=2
      save titley,titlept,titlem,cmom
!$omp threadprivate(first,titley,titlept,titlem,cmom)

c---  on first call for each plot, set up strings      
      if (first(ip,jp,kp)) then
c--- convert integer:: to left-justified character string
        write(cmom,'(i6)') mm
        cmom=adjustl(cmom)
        titlem='m('//trim(cmom)//')'
        titlept='pt('//trim(cmom)//')'   
        titley='y('//trim(cmom)//')'
        first(ip,jp,kp)=.false.
      endif
      
c--- three-particle invariant mass 
      if (tag == tagbook) then
        var=0._dp
      else
        var=sqrt((p(ip,4)+p(jp,4)+p(kp,4))**2
     &           -(p(ip,1)+p(jp,1)+p(kp,1))**2
     &           -(p(ip,2)+p(jp,2)+p(kp,2))**2
     &           -(p(ip,3)+p(jp,3)+p(kp,3))**2)
      endif
      call bookplot(n,tag,titlem(ip,jp,kp),var,wt,wt2,
     &              0._dp,400._dp,5._dp,'lin')
       n=n+1
      
c--- three-particle pt 
      if (tag == tagbook) then
        var=0._dp
      else
        var=ptthree(ip,jp,kp,p)
      endif
      call bookplot(n,tag,titlept(ip,jp,kp),var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
       n=n+1
      
c--- three-particle rapidity 
      if (tag == tagbook) then
        var=0._dp
      else
        var=yrapthree(ip,jp,kp,p)
      endif
      call bookplot(n,tag,titley(ip,jp,kp),var,wt,wt2,
     &              -5._dp,5._dp,0.4_dp,'lin')
       n=n+1
      
      return
      end
    
            
************************************************************************
*                                                                      *
*                 COMPOUND MOMENTUM ip + jp + kp +lp                   *
*                                                                      *
************************************************************************
      subroutine autoplot4(p,mm,ip,jp,kp,lp,tag,wt,wt2,n)
      implicit none
      include 'types.f'
c--- routine that calls histogramming routines to produce a basic
c--- set of plots corresp. to momentum p(ip)+p(jp)+p(kp)+p(lp) == p(mm)
c---  wt, wt2 are passed into the routine and hold information
c---   on the event weight
c---  tag is set by the calling routine, to either initialize or fill
c---  n is the histogram counter that is passed in, incremented
c---   appropriately and returned
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mcfmplotinfo.f'
      integer:: mm,ip,jp,kp,lp,n
      real(dp):: p(mxpart,4),wt,wt2,var,ptfour,yrapfour
      integer tag
      character*12 titley(10,10,10,10),titlept(10,10,10,10)
     & ,titlem(10,10,10,10)
      character*6 cmom
      logical, save :: first(10,10,10,10)=.true.
      integer, parameter:: tagbook=1, tagplot=2
      save titley,titlept,titlem,cmom
!$omp threadprivate(first,titley,titlept,titlem,cmom)

c--- on first call for each plot, set up strings      
      if (first(ip,jp,kp,lp)) then
c--- convert integer:: to left-justified character string
        write(cmom,'(i6)') mm
        cmom=adjustl(cmom)
        titlem='m('//trim(cmom)//')'
        titlept='pt('//trim(cmom)//')'   
        titley='y('//trim(cmom)//')'
        first(ip,jp,kp,lp)=.false.
      endif
      
c--- three-particle invariant mass 
      if (tag == tagbook) then
        var=0._dp
      else
        var=sqrt((p(ip,4)+p(jp,4)+p(kp,4)+p(lp,4))**2
     &           -(p(ip,1)+p(jp,1)+p(kp,1)+p(lp,1))**2
     &           -(p(ip,2)+p(jp,2)+p(kp,2)+p(lp,2))**2
     &           -(p(ip,3)+p(jp,3)+p(kp,3)+p(lp,3))**2)
      endif
      call bookplot(n,tag,titlem(ip,jp,kp,lp),var,wt,wt2,
     &              0._dp,1000._dp,10._dp,'lin')
       n=n+1
      
c--- four-particle pt 
      if (tag == tagbook) then
        var=0._dp
      else
        var=ptfour(ip,jp,kp,lp,p)
      endif
      call bookplot(n,tag,titlept(ip,jp,kp,lp),var,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
       n=n+1
      
c--- four-particle rapidity 
      if (tag == tagbook) then
        var=0._dp
      else
        var=yrapfour(ip,jp,kp,lp,p)
      endif
      call bookplot(n,tag,titley(ip,jp,kp,lp),var,wt,wt2,
     &              -5._dp,5._dp,0.4_dp,'lin')
       n=n+1
      
      return
      end
            
