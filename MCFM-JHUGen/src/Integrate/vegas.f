	subroutine Vegas_Pomp(integrand,result,absacc,relacc,ndim,ncall,
     >                        maxiter,init)
! Uses kahan summation to guarantee identical results independent of number of threads
	
	include 'types.f'
	include 'nf.f'
	include 'lc.f'
	include 'mxdim.f'
	include 'mxpart.f'
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
        include 'gridinfo.f'
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
	real(dp):: result, relacc, absacc,integrand
	integer:: ndim, ncall, maxiter, neval,init,i
	external integrand
	real(dp):: fun, sfun, sfun2
	real(dp):: sint,sint2,sweight
	real(dp):: fun2, weight
	real(dp):: r, dr, xo, xn, err,ran2
	real(dp):: Ingrid,Incall
	integer:: iter, calls, dim, grid, g, c, cmax
	integer:: ngrid,j,jj
	parameter (ngrid = 10)
	real(dp):: xi(ngrid, MXDIM), d(ngrid, MXDIM)
	real(dp):: x(mxdim), imp(ngrid), tmp(ngrid - 1)
	integer:: pos(NDIM)
	real(dp):: t,dt,cfun,cfun2,cd(ngrid,MXDIM)
        real(dp):: p1ext(4),p2ext(4)

	character*255 runname
	integer:: nlength
	logical:: bin,dryrun
	common/bin/bin
	common/runname/runname
	common/nlength/nlength
	common/dryrun/dryrun
        common/pext/p1ext,p2ext
        save xi
!$omp threadprivate(/pext/)
	bin=.false.
	if(init<=0)then
	   do i=1,ndim
	      xi(1,i)=1d0
           enddo
*       define the initial distribution of intervals
	   Ingrid=1d0/real(ngrid)
	   do dim = 1, ndim
	      do grid = 1, ngrid
	         r = real(grid)*Ingrid
	         xi(grid, dim) = r
	      enddo
	   enddo
	endif
	if (init<=1)then
	   neval = 0
	   sint=0d0
	   sweight=0d0
	   sint2=0d0
	   iter=0
	endif


c--- read-in grid if necessary
        if (readin) then
           if (dryrun) then
	      open(unit=11,file=ingridfile,status='unknown')
           else
	      open(unit=11,file=runname(1:nlength)//'_'
     &             //ingridfile,status='unknown')
           endif

	   write(6,*)'*********************************************'
	   write(6,*)'* Reading in vegas grid from ',
     &     	runname(1:nlength)//'_'//ingridfile,' *'
           write(6,*)'*********************************************'
           call flush(6)
           do j=1,ndim
             read(11,203) jj,(xi(i,j),i=1,ngrid)
           enddo
           close(11)
           readin=.false.
	endif

*       iterations loop
 1	continue
	iter = iter + 1
*       initialize iteration variables
	sfun = 0d0
	sfun2 = 0d0
	d(:,:)=0d0
	cfun = 0d0
	cfun2 = 0d0
	cd(:,:)=0d0
	Incall=1d0/real(ncall)
!$omp  parallel do
!$omp& schedule(dynamic)
!$omp& default(private)
!$omp& shared(incall,xi,ncall,ndim,sfun,sfun2,d,cfun,cfun2,cd)
!$omp& copyin(/xmin/,/taumin/,/cutoff/,/jetcuts/,/breit/,/zerowidth/)
!$omp& copyin(/srdiags/,/vsymfact/,/qcdcouple/,/ewcouple/,/masses/)
!$omp& copyin(/interference/,/facscale/,/mcfmscale/,/stopscales/)
!$omp& copyin(/ipsgen/,/pext/,/ColC/,/qlfirst/,/ptildes/)
!$omp& copyin(/bitflags/,/flags/,/lastphot/)
!$omp& copyin(/QCDb0/,/nodecay/,/swapxz/,/heavyflav/)
!$omp& copyin(/notag/,/nflav/,/reset/)
	do calls = 1, ncall
	   weight = Incall
!$omp critical
           do i=1,ndim
              x(i)=ran2()
           enddo
!	   call GetRandom(x)
!$omp end critical
	   do dim = 1, ndim+2
	      r = x(dim)*ngrid + 1
	      grid = int(r)
	      xo = 0
	      if( grid > 1 ) xo = xi(grid - 1, dim)
	      xn = xi(grid, dim) - xo
	      x(dim) = xo + (r - grid)*xn
	      pos(dim) = grid
	      weight = weight*xn*ngrid
	   enddo
*       compute the integrand
	   fun=integrand(x,weight)
	   fun = fun*weight
	   fun2 = fun**2
!$omp critical
           t=sfun+fun
	   if (abs(sfun)>=abs(fun)) then
	      cfun=cfun+((sfun-t)+fun)
	   else
	      cfun=cfun+((fun-t)+sfun)
	   endif
	   sfun=t
           t=sfun2+fun2
	   if (abs(sfun2)>=abs(fun2)) then
	      cfun2=cfun2+((sfun2-t)+fun2)
	   else
	      cfun2=cfun2+((fun2-t)+sfun2)
	   endif
	   sfun2=t
	   do dim = 1, ndim
	      i=pos(dim)
	      dt=d(i,dim)
	      t=dt+fun2
	      if (abs(dt)>=abs(fun2)) then
		 cd(i,dim)=cd(i,dim)+((dt-t)+fun2)
	      else
		 cd(i,dim)=cd(i,dim)+((fun2-t)+dt)
	      endif
	      d(i,dim)=t
	   enddo
!$omp end critical
 666	   continue
	enddo
!$omp end parallel do
	sfun=sfun+cfun
	sfun2=sfun2+cfun2
	d(:,:)=d(:,:)+cd(:,:)
	neval = neval + ncall
*       compute the integral and error values
	err = 0
	fun2 = sfun**2
	sint2 = sint2 + fun2
	r = sfun2*ncall - fun2
	if( r .ne. 0 ) then
	   weight = fun2/abs(r)*(ncall - 1)
	   sweight = sweight + weight
	   sint = sint + sfun*weight
	endif
	if( sweight == 0 ) then
	   result = 0
	else
	   r = sint/sweight
	   result = r
*       if the integrand is very close to zero, it is pointless (and costly)
*       to insist on a certain relative accuracy
	   if( abs(r) > absacc ) then
              r = sint2/(sint*r)
	      if( r > err ) then
		 err = r
	      endif
	   endif
	endif
	err = sqrt(err/iter)

        print *, "iteration ", iter, ":",result,"+/-",result*err
        if (abs(result) < 1d-9) then
          print *, "integral is zero, exiting"
          return
        endif
!	if( err <= relacc ) then
!	   print *, "iteration ", iter, ":",result,"+/-",result*err
!	   call Outhist(iter)
!	   return
!	endif

*       redefine the grid (importance sampling)
*       - smooth the f^2 value stored for each interval
	do dim = 1, ndim
	   xo = d(1, dim)
	   xn = d(2, dim)
	   d(1, dim) = .5D0*(xo + xn)
	   x(dim) = d(1, dim)
	   do grid = 2, ngrid - 1
	      r = xo + xn
	      xo = xn
	      xn = d(grid + 1, dim)
	      d(grid, dim) = (r + xn)/3D0
	      x(dim) = x(dim) + d(grid, dim)
	   enddo
	   d(ngrid, dim) = .5D0*(xo + xn)
	   x(dim) = x(dim) + d(ngrid, dim)
	enddo

*       - compute the importance function of each interval
	do dim = 1, ndim
	   r = 0
	   do grid = 1, ngrid
	      imp(grid) = 0
	      if( d(grid, dim) > 0 ) then
		 xo = x(dim)/d(grid, dim)
		 imp(grid) = ((xo - 1)/xo/log(xo))**1.5D0
	      endif
	      r = r + imp(grid)
	   enddo
	   r = r/ngrid

*       - redefine the size of each interval
	   dr = 0
	   xn = 0
	   g = 0
	   do grid = 1, ngrid - 1
	      do while( dr < r )
		 g = g + 1
		 dr = dr + imp(g)
		 xo = xn
		 xn = xi(g, dim)
	      enddo
	      dr = dr - r
	      tmp(grid) = xn - (xn - xo)*dr/imp(g)
	   enddo
	   do grid = 1, ngrid - 1
	      xi(grid, dim) = tmp(grid)
	   enddo
	   xi(ngrid, dim) = 1
	enddo
c--- added to write out intermediate results
!	if ((bin) .and. (iter < maxiter)) then
!	   write(6,*)'Writing out intermediate results for iteration',iter
!	   call histofin(tgral,sd,iter,maxiter) 
!	endif
 28	continue
c--- write-out grid if necessary
	if (writeout) then
           open(unit=11,file=runname(1:nlength)//'_'
     &           //outgridfile,status='unknown')
	   write(6,*)'***********************************************'
	   write(6,*)'* Writing out vegas grid to ',
     &               runname(1:nlength)//'_'//outgridfile,'  *'
           write(6,*)'***********************************************'
           call flush(6)
           do j=1,ndim
             write(11,203) j,(xi(i,j),i=1,ngrid)
           enddo
           close(11)
         endif

	if( iter >= maxiter ) then
!	   print *,"iterations reached set maximum of ",maxiter
!	   print *, "iteration ", iter, ":",result,"+/-",result*err
	   return
	endif

	goto 1
203     FORMAT(/(5z16))

	end








************************************************************************
*       * IniRandom sets up the random-number generator to produce at most
*       * max dims-dimensional quasi-random vectors

	subroutine IniRandom(dims)
	
	integer:: max, dims

	integer:: ndim,a,b,ni(55),n(55)
	data ni /
     1  980629335, 889272121, 422278310,1042669295, 531256381,
     2  335028099,  47160432, 788808135, 660624592, 793263632,
     3  998900570, 470796980, 327436767, 287473989, 119515078,
     4  575143087, 922274831,  21914605, 923291707, 753782759,
     5  254480986, 816423843, 931542684, 993691006, 343157264,
     6  272972469, 733687879, 468941742, 444207473, 896089285,
     7  629371118, 892845902, 163581912, 861580190,  85601059,
     8  899226806, 438711780, 921057966, 794646776, 417139730,
     9  343610085, 737162282,1024718389,  65196680, 954338580,
     1  642649958, 240238978, 722544540, 281483031,1024570269,
     2  602730138, 915220349, 651571385, 405259519, 145115737 /
	common /rngdata/ ndim,a,b,n

	a=55
	b=31
	n(:)=ni(:)
	ndim = dims
	end


************************************************************************
*       * GetRandom is a subtractive Mitchell-Moore random-number generator.
*       * The algorithm is n(i) = (n(i - 24) - n(i - 55)) mod m, implemented
*       * as a circular array with n(i + 55) = n(i) and m = 2^30 in this
*       * version.  The array n has been initialized by setting n(i) = i and
*       * running the algorithm 100,000 times.  Code by Ronald Kleiss.

	subroutine GetRandom(array)
        include 'types.f'	
 	real(dp):: array(*)
	integer:: ndim
	integer:: dim, a, b, j, m, n(55)
	common /rngdata/ ndim,a,b,n
	parameter (m = 2**30)

	do dim = 1, ndim
	   a = mod(a, 55) + 1
	   b = mod(b, 55) + 1
	   j = n(b) - n(a)
	   if( j < 0 ) j = j + m
	   n(a) = j
	   array(dim) = real(j)/m
	enddo
	end
