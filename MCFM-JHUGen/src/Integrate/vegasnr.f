      SUBROUTINE vegasnr(region,ndim,fxn,init,ncall,itmx,nprn,tgral,sd,
     *chi2a)
      
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mpif.h'
      include 'mxdim.f'
      include 'gridinfo.f'
      include 'maxwt.f'
      include 'TRtensorcontrol.f'
      include 'tensorinfo.f'
      include 'phasemin.f'
      include 'cutoff.f'
      include 'jetcuts.f'
      include 'breit.f'
      include 'zerowidth.f'
      include 'srdiags.f'
      include 'interference.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ipsgen.f'
      include 'facscale.f'
      include 'scale.f'
      include 'stopscales.f'
      include 'qlfirst.f'
      include 'ptilde.f'
      include 'bitflags.f'
      include 'flags.f'
      include 'lastphot.f'
      include 'b0.f'
      include 'nodecay.f'
      include 'swapxz.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'notag.f'
      include 'reset.f'
      include 'lc.f'
      include 'TRmaxindex.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'histo.f'
      include 'hbbparams.f'
      include 'mpicommon.f'
      integer:: init,itmx,itmxloop,ndim,nprn,NDMX
      integer(kind=8) ncall,ng,npg,nd,ndo,k
      real(dp):: tgral,chi2a,sd,region(2*mxdim),fxn,ALPH,TINY
c--- Note: NDMX increased to 100 (from 50) compared with versions 5.1 and
c---  earlier, to aid calculation of H+2 jets process
      PARAMETER (ALPH=1.5_dp,NDMX=100,TINY=1e-30_dp)
c--- DEBUG: no adapting
c      PARAMETER (ALPH=zip,NDMX=100,TINY=1e-30_dp)
      EXTERNAL fxn
C     USES fxn,ran2,rebin
      integer:: i,l,idum,it,j,jj,mds,ia(MXDIM),kg(MXDIM)
      real(dp):: calls,dv2g,dxg,f,f2,rc,ti,tsi,wgt,xjac,
     *xxn,xnd,xo,dt(MXDIM),dx(MXDIM),r(NDMX),x(MXDIM),
     *xin(NDMX),ran2,ran2nr,dumran
      double precision xi(NDMX,MXDIM),fb,fbr,f2b,f2br,d(NDMX,MXDIM),
     *dr(NDMX,MXDIM),di(NDMX,MXDIM),dir(NDMX,MXDIM)
      real(dp):: schi,si,swgt,p1ext(4),p2ext(4)
      real(dp):: t,ct,cfun,cfun2,cd(NDMX,MXDIM),cdi(NDMX,MXDIM)
      character*255 runname
      integer:: nlength
      logical:: bin,dorebin,dryrun,kahan
      common/bin/bin
      common/runname/runname
      common/nlength/nlength
      common/dryrun/dryrun
      common/pext/p1ext,p2ext
!$omp threadprivate(/pext/)
      parameter(kahan=.false.)
      integer ierr
      SAVE


      dorebin=.true.  
      mds=0
      if(init<=0)then
        ndo=1
        do 11 j=1,ndim
          xi(1,j)=1d0
11      continue
      endif
      if (init<=1)then
        si=zip
        swgt=zip
        schi=zip
      endif
      if (init<=2)then
        nd=NDMX
        ng=1
c--- DEBUG
c        write(6,*) 'DEBUG: Setting mds to zero'
c        mds=0
c--- DEBUG
        if(mds.ne.0)then
          ng=(ncall/two+0.25_dp)**(one/ndim)
          mds=1
          if((2*ng-NDMX)>=0)then
            mds=-1
            npg=ng/NDMX+1
            nd=ng/npg
            ng=npg*nd
          endif
        endif
        k=ng**ndim
        npg=max(ncall/k,2)
        calls=npg*k
        dxg=one/ng
        dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
        xnd=nd
        dxg=dxg*xnd
        xjac=one/calls
        do 12 j=1,ndim
          dx(j)=region(j+ndim)-region(j)
          xjac=xjac*dx(j)
12      continue

c--- read-in grid if necessary
        if (rank.eq.0) then
        if (readin) then
           if (dryrun) then
             open(unit=11,file=ingridfile,status='unknown')
           else
             open(unit=11,file=runname(1:nlength)//'_'
     &             //ingridfile,status='unknown')
           endif
        write(6,*)'****************************************************'
        write(6,*)'* Reading in vegas grid from ',runname(1:nlength)//
     &   '_'//ingridfile,' *'
        write(6,*)'****************************************************'
           call flush(6)
           do j=1,ndim
             read(11,203) jj,(xi(i,j),i=1,nd)
           enddo
           close(11)
           ndo=nd
           readin=.false.
c--- do not continue to adapt grid when using a small number of calls
           if (calls < 1d3) then
             write(6,*)
             write(6,*) '--> Small number of calls, so the grid is not',
     &        ' being adjusted after each iteration <--'
             write(6,*)
             call flush(6) 
             dorebin=.false.
           endif
        endif

        if(nd.ne.ndo)then
          do 13 i=1,nd
            r(i)=one
13        continue
          do 14 j=1,ndim
            call rebin(ndo/xnd,nd,r,xin,xi(1,j))
14        continue
          ndo=nd
        endif
        if(nprn>=0) write(6,200) ndim,calls,it,itmx,nprn,ALPH,mds,nd,
     *(j,region(j),j,region(j+ndim),j=1,ndim)
        call flush(6)
      endif
      endif

c--- set maximum number of iterations
c--- in case of target error, cap at 100 iterations
      if (itmx < 0) then
        itmxloop=100
      else
        itmxloop=itmx
      endif
      
      do 28 it=1,itmxloop
        ti=zip
        tsi=zip
        do 16 j=1,ndim
          kg(j)=1
          do 15 i=1,nd
            dr(i,j)=zip
            dir(i,j)=zip
15        continue
16      continue
10      continue
        fbr=zip
        f2br=zip
        if (kahan) then
           cfun = zip
           cfun2 = zip
           cd(:,:)=zip
           cdi(:,:)=zip
        endif
        call mpi_bcast(xi,NDMX*MXDIM,mpi_double_precision,
     .                 0,mpi_comm_world,ierr)
!$omp  parallel do
!$omp& schedule(dynamic)
!$omp& default(private)
!$omp& shared(npg,xjac,ncall,ndim,kg,dxg,xi,region,dx,xnd)
!$omp& shared(fbr,f2br,dir,dr,cfun,cfun2,cd,cdi,rank,size)
!$omp& copyin(/xmin/,/taumin/,/cutoff/,/jetcuts/,/breit/,/zerowidth/)
!$omp& copyin(/srdiags/,/vsymfact/,/qcdcouple/,/ewcouple/,/masses/)
!$omp& copyin(/interference/,/facscale/,/mcfmscale/,/stopscales/)
!$omp& copyin(/ipsgen/,/pext/,/ColC/,/qlfirst/,/ptildes/)
!$omp& copyin(/bitflags/,/flags/,/lastphot/)
!$omp& copyin(/QCDb0/,/nodecay/,/swapxz/,/heavyflav/)
!$omp& copyin(/notag/,/nflav/,/reset/,/pvmaxindex/)
!$omp& copyin(/epinv/,/epinv2/,/plotindex/,/hbbparams/)
          do 19 k=1,npg/size
            wgt=xjac
!$omp critical(VegasRandom)
c--- note: added (ndim+1) to use as random variable for PS branches 
c--- note: added (ndim+2) to use as random variable in new_pspace
            do j=1,size
               do l=1,ndim+2
                  dumran=ran2nr()
                  if (j.eq.rank+1) x(l)=dumran
!           x(l)=one-ran2() ! this way, agrees with vegas_Pomp
               enddo
            enddo
!$omp end critical(VegasRandom)
        do 17 j=1,ndim
              xxn=(kg(j)-x(j))*dxg+one
              ia(j)=max(min(int(xxn),NDMX),1)
              if(ia(j)>1)then
                xo=xi(ia(j),j)-xi(ia(j)-1,j)
                rc=xi(ia(j)-1,j)+(xxn-ia(j))*xo
              else
                xo=xi(ia(j),j)
                rc=(xxn-ia(j))*xo
              endif
              x(j)=region(j)+rc*dx(j)
              wgt=wgt*xo*xnd
 17        continue
            f=wgt*fxn(x,wgt)
            f2=f*f
!$omp critical(VegasAdd)
            if (kahan) then
               t=fbr+f
               if (abs(fbr)>=abs(f)) then
                  cfun=cfun+((fbr-t)+f)
               else
                  cfun=cfun+((f-t)+fbr)
               endif
               fbr=t
               t=f2br+f2
               if (abs(f2br).ge.abs(f2)) then
                  cfun2=cfun2+((f2br-t)+f2)
               else
                  cfun2=cfun2+((f2-t)+f2br)
               endif
               f2br=t
               do j = 1, ndim
                  i=ia(j)
                  ct=dir(i,j)
                  t=ct+f
                  if (abs(ct)>=abs(f)) then
                     cdi(i,j)=cdi(i,j)+((ct-t)+f)
                  else
                     cdi(i,j)=cdi(i,j)+((f-t)+ct)
                  endif
                  dir(i,j)=t
               enddo
               if (mds>=0) then
                  do j = 1, ndim
                     i=ia(j)
                     ct=dr(i,j)
                     t=ct+f2
                     if (abs(ct)>=abs(f2)) then
                        cd(i,j)=cd(i,j)+((ct-t)+f2)
                     else
                        cd(i,j)=cd(i,j)+((f2-t)+ct)
                     endif
                     dr(i,j)=t
                  enddo
               endif
            else
               fbr=fbr+f
               f2br=f2br+f2
               do 18 j=1,ndim
                  dir(ia(j),j)=dir(ia(j),j)+f
                  if(mds.ge.0) dr(ia(j),j)=dr(ia(j),j)+f2
 18            continue
            endif
!$omp end critical(VegasAdd)
 19    continue
!$omp end parallel do

       if (kahan) then
          fbr=fbr+cfun
          f2br=f2br+cfun2
          dir(:,:)=dir(:,:)+cdi(:,:)
          dr(:,:)=dr(:,:)+cd(:,:)
       endif
       call mpi_reduce(fbr,fb,1,mpi_double_precision,
     .                 mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(f2br,f2b,1,mpi_double_precision,
     .                 mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(dr,d,NDMX*MXDIM,mpi_double_precision,
     .                 mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(dir,di,NDMX*MXDIM,mpi_double_precision,
     .                 mpi_sum,0,mpi_comm_world,ierr)
       if (rank.eq.0) then
          f2b=sqrt(f2b*npg)
          f2b=(f2b-fb)*(f2b+fb)
          if (f2b<=zip) f2b=TINY
          ti=ti+fb
          tsi=tsi+f2b
          if(mds<0)then
            do 21 j=1,ndim
              d(ia(j),j)=d(ia(j),j)+f2b
21          continue
         endif
        do 22 k=ndim,1,-1
          kg(k)=mod(kg(k),ng)+1
          if(kg(k).ne.1) goto 10
22      continue
        tsi=tsi*dv2g
        wgt=one/tsi
        si=si+real(wgt,dp)*real(ti,dp)
        schi=schi+real(wgt,dp)*real(ti,dp)**2
        swgt=swgt+real(wgt,dp)
        tgral=si/swgt
        chi2a=max((schi-si*tgral)/(it-.99_dp),zip)
        sd=sqrt(one/swgt)
        tsi=sqrt(tsi)

        if(nprn.ge.0)then
c          write(6,201) it,ti,tsi,tgral,sd,chi2a
          write(6,201) it,ti,tgral,tsi,sd,wtmax,chi2a
          if (TRtensorcontrol > 0) then
          write(6,301) real(ibadpoint,dp)/real(itotal,dp),
     &                 real(ipolesfailed,dp)/real(itotal,dp),
     &                 real(ibadpoint+ipolesfailed,dp)/real(itotal,dp)
          endif
          call flush(6)
         if(nprn.ne.0)then
            do 23 j=1,ndim
              write(6,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
23          continue
          endif
        endif
        
c--- check for termination if running with target error
        if (itmx < 0) then
          if (sd < real(-itmx,dp)) exit
        endif
        
        if (abs(tgral) < 1e-9_dp) then
          write(6,*) '******** Integral is zero, no more iterations '
     &      //'required *********'
          write(6,*)
          tgral=zip
          sd=zip
          call flush(6)
          exit ! bail early (28) if result is zero
        endif
        do 25 j=1,ndim
          xo=d(1,j)
          xxn=d(2,j)
          d(1,j)=(xo+xxn)/two
          dt(j)=d(1,j)
          do 24 i=2,nd-1
            rc=xo+xxn
            xo=xxn
            xxn=d(i+1,j)
            d(i,j)=(rc+xxn)/three
            dt(j)=dt(j)+d(i,j)
24        continue
          d(nd,j)=(xo+xxn)/two
          dt(j)=dt(j)+d(nd,j)
25      continue
        do 27 j=1,ndim
           rc=zip
           do 26 i=1,nd
            if(d(i,j)<TINY) d(i,j)=TINY
            r(i)=((one-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**ALPH
            rc=rc+r(i)
26        continue
          if (dorebin) call rebin(rc/xnd,nd,r,xin,xi(1,j))
27      continue
c--- added to write out intermediate results
c      if ((bin) .and. (it .lt. itmx)) then
c        write(6,*) 'Writing out intermediate results for iteration',it
c        call histofin(tgral,sd,it,itmx) 
c      endif
      endif
28    continue

c--- write-out grid if necessary

      if (rank.eq.0) then
         if (writeout) then
            open(unit=11,file=runname(1:nlength)//'_'
     &           //outgridfile,status='unknown')
        write(6,*)'****************************************************'
        write(6,*)'* Writing out vegas grid to ',runname(1:nlength)//'_'
     &           //outgridfile,'  *'
        write(6,*)'****************************************************'
            call flush(6)
            do j=1,ndim
               write(11,203) j,(xi(i,j),i=1,nd)
            enddo
            close(11)
         endif
      endif
      return

200   FORMAT(/' input parameters for vegas:  ndim=',i3,'  ncall=',
     *f13.0/28x,'  it=',i5,'  itmx=',i5/28x,'  nprn=',i3,'  alph=',
     *f5.2/28x,'  mds=',i3,'   nd=',i4/(30x,'xl(',i2,')= ',g11.4,' xu(',
     *i2,')= ',g11.4))
c201   FORMAT(/' iteration no.',I3,': ','integral =',g14.7,'+/- ',g9.2/
c     *' all iterations:   integral =',g14.7,'+/- ',g9.2,' chi**2/iter',
c     * g9.2)
 201    format(/'************* Integration by Vegas (iteration ',i3,
     &   ') **************' / '*',63x,'*'/,
     &   '*  integral  = ',g14.8,2x,
     &   ' accum. integral = ',g14.8,'*'/,
     &   '*  std. dev. = ',g14.8,2x,
     &   ' accum. std. dev = ',g14.8,'*'/,
     &   '*   max. wt. = ',g14.6,35x,'*'/,'*',63x,'*'/,
     &   '**************   chi**2/iteration = ',
     &   g10.4,'   ****************' /)     
301     format('> TensorRed:   TensF',e9.2,
     &    ' / PolF',e9.2,' / TotalF',e9.2,' <'/)
202   FORMAT(/' data for axis ',I2/'    X       delta i       ',
     *'   x       delta i       ','    x       delta i       ',/(1x,
     *f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
203   FORMAT(/(5z16))
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


      subroutine rebin(rc,nd,r,xin,xi)
      implicit none
      include 'types.f'
      integer(kind=8):: nd
      real(dp):: rc,r(*),xin(*)
      double precision xi(*)
      integer:: i,k
      real(dp):: dr,xxn,xo
      k=0
      xxn=0._dp
      dr=0._dp
      do 11 i=1,nd-1
1       if(rc>dr)then
          k=k+1
          dr=dr+r(k)
          xo=xxn
          xxn=xi(k)
        goto 1
        endif
        dr=dr-rc
        xin(i)=xxn-(xxn-xo)*dr/r(k)
11    continue
      do 12 i=1,nd-1
        xi(i)=xin(i)
12    continue
      xi(nd)=1._dp
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
