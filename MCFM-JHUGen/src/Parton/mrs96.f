      subroutine mrs96(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C                                                               C
C     This is a package for the new MRS(R1,R2,R3,R4) parton     C
C     distributions. There are several important changes from   C
C     earlier MRS packages:                                     C
C       -- the q**2 range is enlarged to 1.25d0 < q**2 < 1d7,   C
C          the x range is still 1d-5 < x < 1d0                  C
C       -- the interpolation routine has been slightly modified C
C       -- the call is now to mrs96() rather than to MRSEB()    C 
C     Note that the grid files which the program reads in       C
C     (mrsr1.dat,...) are now larger and more obviously named.  C
C                                                               C
C     As before, x times the parton distribution is returned,   C
C     q is the scale in GeV, MSbar factorization is assumed,    C
C     and Lambda(MSbar,nf=4) = 241 MeV for R1 (mode=1)          C
C                            = 344 MeV for R2 (mode=2)          C
C                            = 241 MeV for R3 (mode=3)          C
C                            = 344 MeV for R4 (mode=4)          C
C                                                               C
C         The reference is:                                     C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         University of Durham preprint DTP/96/44 (1996)        C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE R1 GRID IS 0.00150            C
C         THE FIRST NUMBER IN THE R2 GRID IS 0.00125            C
C         THE FIRST NUMBER IN THE R3 GRID IS 0.00181            C
C         THE FIRST NUMBER IN THE R4 GRID IS 0.00085            C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      save xmin,xmax,qsqmin,qsqmax
!$omp threadprivate(xmin,xmax,qsqmin,qsqmax)
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) then
      print 99
      write(6,*) 'q2=',q2,qsqmin,qsqmax
      endif
      if(x.lt.xmin.or.x.gt.xmax) then
      print 98
      write(6,*) 'x=',x
      endif
         if(mode.eq.1) then
        call mrsr1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrsr2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.3) then
        call mrsr3(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.4) then
        call mrsr4(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      endif 
  99  format('  MRS96:WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  format('  MRS96:WARNING:   X  VALUE IS OUT OF RANGE   ')
      return
      end

      subroutine mrsr1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      character*72 filename,checkpath
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .        1d-4,2d-4,4d-4,6d-4,8d-4,
     .        1d-3,2d-3,4d-3,6d-3,8d-3,
     .        1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .     .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .     .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .     .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .     .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f
!$omp threadprivate(xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f)
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        filename=checkpath('Pdfdata/mrsr1.dat')
        open(unit=1,file=filename,status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
        do 25 i=1,np
  25     f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
      close(unit=1)
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin) qsq=qsqmin
      if(qsq.gt.qsqmax) qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .    +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
          fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
          g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrsr2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      character*72 filename,checkpath
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .        1d-4,2d-4,4d-4,6d-4,8d-4,
     .        1d-3,2d-3,4d-3,6d-3,8d-3,
     .        1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .     .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .     .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .     .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .     .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f
!$omp threadprivate(xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f)
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        filename=checkpath('Pdfdata/mrsr2.dat')
        open(unit=1,file=filename,status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
        do 25 i=1,np
  25     f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      close(unit=1)
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin) qsq=qsqmin
      if(qsq.gt.qsqmax) qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .    +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
          fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
          g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrsr3(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      character*72 filename,checkpath
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .        1d-4,2d-4,4d-4,6d-4,8d-4,
     .        1d-3,2d-3,4d-3,6d-3,8d-3,
     .        1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .     .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .     .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .     .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .     .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f
!$omp threadprivate(xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f)
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        filename=checkpath('Pdfdata/mrsr3.dat')
        open(unit=1,file=filename,status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
        do 25 i=1,np
  25     f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      close(unit=1)
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin) qsq=qsqmin
      if(qsq.gt.qsqmax) qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .    +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
          fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
          g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrsr4(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      character*72 filename,checkpath
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .        1d-4,2d-4,4d-4,6d-4,8d-4,
     .        1d-3,2d-3,4d-3,6d-3,8d-3,
     .        1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .     .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .     .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .     .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .     .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f
!$omp threadprivate(xx,qq,xmin,xmax,qsqmin,qsqmax,n0,init,f)
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        filename=checkpath('Pdfdata/mrsr4.dat')
        open(unit=1,file=filename,status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
        do 25 i=1,np
  25     f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      close(unit=1)
      init=1
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin) qsq=qsqmin
      if(qsq.gt.qsqmax) qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .    +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
          fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
          g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end

