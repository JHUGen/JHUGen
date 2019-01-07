      subroutine mt(x,scale,mode,upv,dnv,sea,str,chm,bot,glu)
      implicit double precision(a-h,o-z)
      dimension a(0:3,0:2,8),f(8)
      character*72 filename,checkpath
      data init/0/
      save init,q0,qcdl,a
!$omp threadprivate(init,q0,qcdl,a)
      if(init.ne.0) goto 10
      if(mode.eq.1)then
        filename=checkpath('Pdfdata/MTS1.DAT')
        open(unit=1,file=filename,status='old')
      elseif(mode.eq.2)then
        filename=checkpath('Pdfdata/MTE1.DAT')
        open(unit=1,file=filename,status='old')
      elseif(mode.eq.3)then
        filename=checkpath('Pdfdata/MTB1.DAT')
        open(unit=1,file=filename,status='old')
      elseif(mode.eq.4)then
        filename=checkpath('Pdfdata/MTB2.DAT')
        open(unit=1,file=filename,status='old')
      elseif(mode.eq.5)then
        filename=checkpath('Pdfdata/MTSN1.DAT')
        open(unit=1,file=filename,status='old')
      endif
c 2=uv 1=dv 3=glue 4=(ubar+dbar)/2 5=sbar 6=cbar 7=bbar 8=ttbar
!$omp critical(CteqRead)
      read(1,49)q0,qcdl
*      write(6,49)q0,qcdl
  49  format(2f6.3)
      do 20 n=0,3
      do 20 m=0,2
      read(1,50)(a(n,m,j),j=1,8)
*      write(6,50)(a(n,m,j),j=1,8)
  50  format(8f6.2)
  20  continue
      close(1)
!$omp end critical(CteqRead)
      init=1
  10  continue
      t=log(log(scale/qcdl)/log(q0/qcdl))
      t2=t*t
      omx=1d0-x
      dl=log(1d0+1d0/x)
      do 40 i=1,8
        a0=a(0,0,i)+a(0,1,i)*t+a(0,2,i)*t2
        a1=a(1,0,i)+a(1,1,i)*t+a(1,2,i)*t2
        a2=a(2,0,i)+a(2,1,i)*t+a(2,2,i)*t2
        a3=a(3,0,i)+a(3,1,i)*t+a(3,2,i)*t2
        f(i)=exp(a0)*x**a1*omx**a2*dl**a3
   40 continue
      upv=f(2)
      dnv=f(1)
      glu=f(3)
      sea=f(4)
      str=f(5)
      chm=f(6)
      bot=f(7)
      top=f(8)
      return
      end
