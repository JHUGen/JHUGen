      subroutine mrstlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 LO parton            C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0201xxx                                  C
C                                                               C
C  There is 1 pdf set corresponding to mode = 1                 C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.220 C
C  corresponding to alpha_s(M_Z) of 0.130                       C
C  This set reads a grid whose first number is 0.02868          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      save xmin,xmax,qsqmin,qsqmax
!$omp threadprivate(xmin,xmax,qsqmin,qsqmax)
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
          if(mode.eq.1) then
        call mrst1lo(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrst1lo(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      character*72 filename,checkpath
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save xx,qq,xmin,xmax,qsqmin,qsqmax,init,f,emb2,emc2
      common /store/ xxl,qql,qqlc,qqlb,cc1,cc2,cc3,cc4,cc6,cc8,
     & ccc,ccb
!$omp threadprivate(xx,qq,xmin,xmax,qsqmin,qsqmax,init,f,emb2,emc2)
!$omp threadprivate(/store/)
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
!$omp critical(MrsRead)
        filename=checkpath('Pdfdata/lo2002.dat')
        open(unit=33,file=filename,status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      close(unit=33)
!$omp end critical(MrsRead)
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppelo1(nx,nq,xxl,qql,f1,cc1)
      call jeppelo1(nx,nq,xxl,qql,f2,cc2)
      call jeppelo1(nx,nq,xxl,qql,f3,cc3)
      call jeppelo1(nx,nq,xxl,qql,f4,cc4)
      call jeppelo1(nx,nq,xxl,qql,f6,cc6)
      call jeppelo1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppelo1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppelo1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppelo2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppelo2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppelo2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppelo2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppelo2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppelo2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppelo2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppelo2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
      
