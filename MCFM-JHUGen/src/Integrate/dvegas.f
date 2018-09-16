         block data vegas_data
         implicit  real(dp):: (a-h,o-z)
         implicit integer:: (i-n)
         include 'vegas_common.f'
         include 'maxwt.f'
         parameter(mprod=50*mxdim)
         common/bveg2/xi(50,mxdim),si,si2,swgt,schi,ndo,it
c
         data ncall/10000/,itmx/15/,nprn/1000/,acc/-1d0/,
     1   xl/mxdim*0d0/
     2   xu/mxdim*1d0/
c
         data XI/mprod*1d0/
c
         data wtmax/0d0/
c
         end

C
C
C        NCALL IS THE NUMBER OF CALLS TO VEGAS.
C        NPRN >  0 VEGAS PRINTS THE RESULTS OF EACH ITERATION.
C        NPRN <= 0 VEGAS PRINTS NOTHING.
C        XL(I) IS LOWER INTEGRATION LIMIT ON I TH AXIS.
C        XU(I) IS UPPER INTEGRATION LIMIT ON I THE AXIS.
c
         subroutine vegas(fxn,avgi,sd,chi2a)
      implicit none
      include 'types.f'
c
c        routine performs n dim Monte Carlo Integration
c        written by G. P. Lepage
c
         implicit real(dp) (a-h,o-z)
         implicit integer:: (i-n)
         include 'vegas_common.f' 
         include 'gridinfo.f'
         include 'maxwt.f'
         parameter(mprod=50*mxdim)
         common/bveg2/xi(50,mxdim),si,si2,swgt,schi,ndo,it
         COMMON/ranno/idum
         dimension d(50,mxdim),di(50,mxdim),xin(50),r(50),
     1   dx(mxdim),dt(mxdim),x(mxdim),kg(mxdim),ia(mxdim)
         data ndmx/50/,alph/1.5d0/,one/1d0/,mds/1/

         if(ndim > mxdim) then
         write(6,*) 'ndim',ndim
         write(6,*) 'mxdim',mxdim
         write(6,*) 'ndim > mxdim'
         stop
         endif

         ndo=1
         do 1 j=1,ndim
 1       xi(1,j)=one
c
         entry vegas1(fxn,avgi,sd,chi2a)
c        initialises  cumulative  variables but not grid
         it=0
         si=0d0
         si2=si
         swgt=si
         schi=si
c
         entry vegas2(fxn,avgi,sd,chi2a)
c        no initialisation
         nd=ndmx
         ng=1
         if(mds==0)go to 2
         ng=int((real(ncall)/2d0)**(1d0/real(ndim)))
         mds=1
         if((2*ng-ndmx)<0)go to 2
         mds=-1
         npg=ng/ndmx+1
         nd=ng/npg
         ng=npg*nd
 2       k=ng**ndim
         npg=ncall/k
         if(npg<2)npg=2
         calls=real(npg*k)
         dxg=one/ng
         dv2g=(calls*dxg**ndim)**2/real(npg)/real(npg)/real(npg-one)
         xnd=real(nd)
         ndm=nd-1
         dxg=dxg*xnd
         xjac=one/calls
         do 3 j=1,ndim
         dx(j)=xu(j)-xl(j)
 3       xjac=xjac*dx(j)
c
c    rebin preserving bin density
c

c--- read-in grid if necessary
         if (readin) then
           open(unit=11,file=ingridfile//'.grid',status='unknown')
        write(6,*)'****************************************************'
        write(6,*)'* Reading in vegas grid from ',ingridfile,'.grid *'
        write(6,*)'****************************************************'
           call flush(6)
           do j=1,ndim
             read(11,203) jj,(xi(i,j),i=1,nd)
           enddo
           close(11)
           ndo=nd
           readin=.false.
         endif

         if(nd==ndo)go to 8
         rc=ndo/xnd
         do 7 j=1,ndim
         k=0
         xn=0d0
         dr=xn
         i=k
 4       k=k+1
         dr=dr+one
         xo=xn
         xn=xi(k,j)
 5       if(rc>dr)go to 4
         i=i+1
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr
         if(i<ndm)go to 5
         do 6 i=1,ndm
 6       xi(i,j)=xin(i)
 7       xi(nd,j)=one
         ndo=nd
c
 8       if(nprn>=0)write(6,200)ndim,calls,it,itmx,acc
     1   ,mds,nd,(xl(j),xu(j),j=1,ndim)
         call flush(6)
c
         entry vegas3(fxn,avgi,sd,chi2a)
c         main integration loop
 9       it=it+1
         ti=0d0
         tsi=ti
         do 10 j=1,ndim
         kg(j)=1
         do 10 i=1,nd
         d(i,j)=ti
 10      di(i,j)=ti
c
 11      fb=0d0
         f2b=fb
         k=0
 12      k=k+1
         wgt=xjac
         do 15 j=1,ndim
         xn=(real(kg(j))-ran1(idum))*dxg+one
         ia(j)=int(xn)
         if(ia(j)>1)go to 13
         xo=xi(ia(j),j)
         rc=(xn-real(ia(j)))*xo
         go to 14
13       xO=xi(ia(j),j)-xi(ia(j)-1,j)
         rc=xi(ia(j)-1,j)+(xn-real(ia(j)))*xo
 14      x(j)=xl(j)+rc*dx(j)
 15      wgt=wgt*xo*xnd
c
         f=wgt
c         write(6,FMT='(a20,2F20.16)') 'xo,xnd in dvegas: ',xo,xnd
         f=f*fxn(x,wgt)
         f2=f*f
         fb=fb+f
         f2b=f2b+f2
         do 16 j=1,ndim
         di(ia(j),j)=di(ia(j),j)+f
 16      if(mds>=0)d(ia(j),J)=d(ia(j),J)+f2
         if(k<npg) go to 12
c
888    FORMAT(1X,'F',G14.6,'F2',G14.6,'FB',G14.6,'F2B',G14.6)
         f2b= sqrt(f2b*      NPG)
         f2b=(f2b-fb)*(f2b+fb)
1661   FORMAT(1X,'F2B',G14.6,'NPG',  I10)
         ti=ti+fb
         tsi=tsi+f2b
33     FORMAT(1X,'TSI',G14.6,'F2B',G14.6)
         if(mds>=0)go to 18
         do 17 j=1,ndim
 17      d(ia(j),j)=d(ia(j),j)+f2b
 18      k=ndim
 19      kg(k)=mod(kg(k),ng)+1
         if(kg(k).ne.1)go to 11
         k=k-1
         if(k>0)go to 19
c
c final results for this iteration
c
        tsi=tsi*dv2g
        ti2=ti*ti
88     format(1x,'tsi',g14.6)
        wgt=ti2/tsi
        si=si+ti*wgt
        si2=si2+ti2
        swgt=swgt+wgt
        schi=schi+ti2*wgt
995    FORMAT(1X,'SWGT',G14.6,'SI2',G14.6)
        avgi=si/swgt
        sd=swgt*real(it)/si2
        chi2a=sd*(schi/swgt-avgi*avgi)/(real(it)-.999d0)
        sd=sqrt(one/sd)
c
        if(nprn==0)go to 21
        tsi=sqrt(tsi)
c        write(6,201)it,ti,tsi,avgi,sd,chi2a
        write(6,201)it,ti,avgi,tsi,sd,wtmax,chi2a
        call flush(6)
        if(nprn>=0)go to 21
        do 20 j=1,ndim
 20     write(6,202) j,(xi(i,j),di(i,j),d(i,j),i=1,nd)
c
c      refine grid
c
 21     do 23 j=1,ndim
        xo=d(1,j)
        xn=d(2,j)
        d(1,j)=(xo+xn)/2d0
        dt(j)=d(1,j)
        do 22 i=2,ndm
        d(i,j)=xo+xn
        xo=xn
        xn=d(i+1,j)
        d(i,j)=(d(i,j)+xn)/3d0
 22     dt(j)=dt(j)+d(i,j)
        d(nd,j)=(xn+xo)/2d0
 23     dt(j)=dt(j)+d(nd,j)
c
        do 28 j=1,ndim
        rc=0d0
        do 24 i=1,nd
        r(i)=0d0
        if(d(i,j)<=0d0)go to 24
        xo=dt(j)/d(i,j)
        r(i)=((xo-one)/xo/log(xo))**alph
 24     rc=rc+r(i)
        rc=rc/xnd
        k=0
        xn=0d0
        dr=xn
        i=k
 25     k=k+1
        dr=dr+r(k)
        xo=xn
        xn=xi(k,j)
 26     if(rc>dr)go to 25
        i=i+1
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
        if(i<ndm)go to 26
        do 27 i=1,ndm
 27     xi(i,j)=xin(i)
 28     xi(nd,j)=one
c
        if(it<itmx.and.acc*abs(avgi)<sd)go to 9

c--- write-out grid if necessary
         if (writeout) then
           open(unit=11,file=outgridfile//'.grid',status='unknown')
        write(6,*)'****************************************************'
        write(6,*)'* Writing out vegas grid to ',outgridfile,'.grid  *'
        write(6,*)'****************************************************'
           call flush(6)
           do j=1,ndim
             write(11,203) jj,(xi(i,j),i=1,nd)
           enddo
           close(11)
         endif

 200    format(/ 1X,' Input parameters for vegas:  ndim=',i3,
     1  '   ncall=',f8.0/28x,'  it=',i5,'    itmx=',i5/28x,
     2  '  acc=',g9.3/28x,'  mds=',i3,'     nd=',i4/28x,
     3  '  (xl,xu)=',(t40,'( ',g12.6,' , ',g12.6,' )'))
c 201    format(///' Integration by vegas' / ' iteration no.',i3,
c     1  ':  integral=',g14.8/21x,'std dev =',g14.8 /
c     2  ' accumulated results:   integral=',g14.8/
c     3  24x,'std dev =',g14.8 / 24x,'chi**2 per it''n =',g10.4)
 201    format(/'************* Integration by Vegas (iteration ',i3,
     &   ') **************' / '*',63x,'*'/,
     &   '*  integral  = ',g14.8,2x,
     &   ' accum. integral = ',g14.8,'*'/,
     &   '*  std. dev. = ',g14.8,2x,
     &   ' accum. std. dev = ',g14.8,'*'/,
     &   '*   max. wt. = ',g14.6,35x,'*'/,'*',63x,'*'/,
     &   '**************   chi**2/iteration = ',
     &   g10.4,'   ****************' /)     
 202    format(1X,' data for axis',i2,/,' ',6x,'x',7x,'  delt i ',
     1  2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  '
     2  ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/,
     3  (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))
 203    format(/(5z16))
        return
        end

        subroutine save(ndim)
      implicit none
      include 'types.f'
        implicit real(dp) (a-h,o-z)
        implicit integer:: (i-n)
        include 'mxdim.f'
        common/bveg2/xi(50,mxdim),si,si2,swgt,schi,ndo,it
c
c       stores vegas data   (unit 7) for later initialisation
c
        write(7,200) ndo,it,si,si2,swgt,schi,
     1       ((xi(i,j),i=1,ndo),j=1,ndim)
        return
        entry restr(ndim)
c
c     enters initialisation data for vegas
c
        read(7,200) ndo,it,si,si2,swgt,schi,
     1    ((xi(i,j),i= 1,ndo),j=1,ndim)
 200    format(2i8,4z16/(5z16))
        return
        end



