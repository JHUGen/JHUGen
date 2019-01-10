      subroutine ovBtensor(p1,m0sq,m1sq,FB0,FB1,FB2,B00)
C     p1 is the external momenta,
C     m1s,m2s are the squares of the internal masses
C     FC0...FB2 are the rank 0,...2 bubble functions
C     Lorentz indices are stored as linear array, thus FD2(y2(n1,n2),ep)
C     Author: R.K.Ellis (January 2013)
C     Implementing the formula of Denner and Dittmaier arXiv:hep-ph/0509141 
      implicit none
      include 'TRconstants.f'
      include 'TRonshellcutoff.f'
      include 'TRscale.f'
      include 'TRmetric.f'
      include 'TRydef.f'
      include 'TRclear.f'
      include 'ovBnames.f'
      include 'ovBsave.f'
      integer N,ep,nl,in,iP,n1,n2
      double precision p1sq,m0sq,m1sq,f1,iep,p1(4)
      double complex A0(-2:0),B0(-2:0),B1(-2:0),B00(-2:0),B11(-2:0),
     & FB0(-2:0),FB1(y1max,-2:0),FB2(y2max,-2:0),trI1,
     & xp,xm,rt,arg,arg1,pvfndd,cln,xpvfndd
      logical p1sqnonzero
      double precision fac,facnp
      double precision,save::idp1(0:2),id(0:2),idm1(0:2)
      logical,save:: scaleset=.false.
      logical,save:: first=.true.
      double precision para(Pbb)
      double precision,save::tableB(Pbb,Nbmax)      
      integer, save:: Nstore=0
      integer :: jtable,j,Ntrue
!$omp threadprivate(scaleset,first,idp1,id,idm1,tableB,Nstore)

C-----statement functions
      fac(n)=(-1d0)**n/dfloat(n+1)
      facnp(in,iP)=(-1d0)**(iP-2*in-1)
C-----end statement functions
      
      if (clear(2)) then
      clear(2)=.false.
      Nstore=0
      endif

      if (Nstore .gt. Nbmax) then
      print * 
      print *, 'ovBtensor: Nstore .gt. Nbmax'
      print *, 'Nstore,Nbmax',Nstore,Nbmax
      print *, 'Either adjust Nbmax in Bnames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      do j=1,4
      para(j)=p1(j)
      enddo
      para(5)=m0sq
      para(6)=m1sq
C if parameter set is found set pvBcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Pbb
        if (abs(para(j)-tableB(j,jtable)) .lt. 1d-8) then
          Ntrue=Ntrue+1
        else
          exit
        endif 
        enddo
        if (Ntrue .eq. Pbb) then
c--- retrieve from cache
c          write(6,*) 'Retrieving from cache: ',jtable
          do ep=-2,0
            FB0(ep)=FB0save(jtable,ep)
            B00(ep)=B00save(jtable,ep)
            do j=1,y1max
              FB1(j,ep)=FB1save(jtable,j,ep)
            enddo
            do j=1,y2max
              FB2(j,ep)=FB2save(jtable,j,ep)
            enddo
          enddo
          return
        endif
      enddo

C    if parameter set is not found we have to calculate
 20   continue
      Nstore=Nstore+1
      do j=1,Pbb
      tableB(j,Nstore)=para(j)
      enddo
c      write(6,*) 'Computing new Nstore: ',Nstore

      if (first) then
      first=.false.
      call ovarraysetup
C--idp1=1/[D+1]
      idp1(0)=0.2d0
      idp1(1)=idp1(0)*0.4d0
      idp1(2)=idp1(1)*0.4d0
C--id=1/D
      id(0)=0.25d0
      id(1)=id(0)*0.5d0
      id(2)=id(1)*0.5d0
C--idm1=1/[D-1]
      idm1(0)=1d0/3d0
      idm1(1)=idm1(0)*2d0/3d0
      idm1(2)=idm1(1)*2d0/3d0
      endif
      p1sq=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2

      if (scaleset .neqv. .true.) then
      scaleset=.true.
      if ((scale .eq. -1d12) .and. (musq .eq. -1d12)) then
      write(6,*) 'Did you forget to call setmudim?'
      write(6,*) 'Setting scale to scale=1d0'
      scale=1d0
      musq=1d0
      endif
      endif

C----self energies never contain double poles -- set to zero
      B0(-2)=czip
      B1(-2)=czip
      B11(-2)=czip
      B00(-2)=czip

C----calculate B0
c      do ep=-1,0
c      Bv(bb0+N,ep)=trI2(p1sq,m0sq,m1sq,musq,ep)
c      enddo

      if ((abs(p1sq/musq) .lt. onshellcutoff) 
     . .and. (abs(m0sq/musq) .lt. onshellcutoff) 
     . .and. (abs(m1sq/musq) .lt. onshellcutoff)) then
      
      
      do ep=-1,0
      B0(ep)=czip
      B1(ep)=czip
      B00(ep)=czip
      B11(ep)=czip
      enddo 
      goto 99

      elseif (abs(m0sq/musq) .lt. onshellcutoff) then

C---deal with special cases for m0sq=0, p1sq=m1sq, DD(4.13)
      if (abs((p1sq-m1sq)/musq) .lt. onshellcutoff) then
      nl=0
      B0(-1)=dcmplx(fac(nl))
      B0(0)=B0(-1)
     . *dcmplx(log(musq/m1sq)+2d0/dfloat(nl+1))

      nl=1
      B1(-1)=dcmplx(fac(nl))
      B1(0)=B1(-1)
     . *dcmplx(log(musq/m1sq)+2d0/dfloat(nl+1))

      nl=2
      B11(-1)=dcmplx(fac(nl))
      B11(0)=B11(-1)
     . *dcmplx(log(musq/m1sq)+2d0/dfloat(nl+1))

      else
C---deal with special cases for m0sq=0, DD(4.12)
      arg=dcmplx(1d0-m1sq/p1sq)
      arg1=dcmplx(m1sq-p1sq)
      iep=sign(1d0,p1sq)

c--- if p1sq=0 too, root is formally infinite and pvfndd should be zero
      if (abs(p1sq/musq) .gt. onshellcutoff) then
        p1sqnonzero=.true.
      else
        p1sqnonzero=.false.
        xpvfndd=czip
      endif

      nl=0
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      B0(-1)=dcmplx(fac(nl))
      B0(0)=B0(-1)
     . *(dcmplx(log(musq)+1d0/dfloat(nl+1))-cln(arg1,-1d0)-xpvfndd)
      nl=1
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      B1(-1)=dcmplx(fac(nl))
      B1(0)=B1(-1)
     . *(dcmplx(log(musq)+1d0/dfloat(nl+1))-cln(arg1,-1d0)-xpvfndd)
      nl=2
      if (p1sqnonzero) xpvfndd=pvfndd(nl,arg,iep)
      B11(-1)=dcmplx(fac(nl))
      B11(0)=B11(-1)
     . *(dcmplx(log(musq)+1d0/dfloat(nl+1))-cln(arg1,-1d0)-xpvfndd)

      endif

      elseif (abs(p1sq/musq) .lt. onshellcutoff) then
C---deal with special case, p1sq=0
      xp=dcmplx(m0sq/(m0sq-m1sq))  ! other root is formally infinite

      nl=0
      B0(-1)=dcmplx(fac(nl))
      B0(0)=B0(-1)
     . *(dcmplx(log(musq/m0sq))-pvfndd(nl,xp,1d0))
      nl=1
      B1(-1)=dcmplx(fac(nl))
      B1(0)=B1(-1)
     . *(dcmplx(log(musq/m0sq))-pvfndd(nl,xp,1d0))
      nl=2
      B11(-1)=dcmplx(fac(nl))
      B11(0)=B11(-1)
     . *(dcmplx(log(musq/m0sq))-pvfndd(nl,xp,1d0))
      
      else
C----general case, DD (4.8)
      rt=sqrt(dcmplx((m1sq-m0sq-p1sq)**2-4d0*p1sq*m0sq))
      xp=0.5d0*(-dcmplx(m1sq-m0sq-p1sq)+rt)/p1sq
      xm=0.5d0*(-dcmplx(m1sq-m0sq-p1sq)-rt)/p1sq

      nl=0
      B0(-1)=dcmplx(fac(nl))
      B0(0)=B0(-1)
     . *(dcmplx(log(musq/m0sq))-pvfndd(nl,xp,1d0)-pvfndd(nl,xm,-1d0))
      nl=1
      B1(-1)=dcmplx(fac(nl))
      B1(0)=B1(-1)
     . *(dcmplx(log(musq/m0sq))-pvfndd(nl,xp,1d0)-pvfndd(nl,xm,-1d0))
      nl=2
      B11(-1)=dcmplx(fac(nl))
      B11(0)=B11(-1)
     . *(dcmplx(log(musq/m0sq))-pvfndd(nl,xp,1d0)-pvfndd(nl,xm,-1d0))

      endif

c--- Construct tensors involving the metric (these are general)
c---  following the formulae of DD (4.4)
      f1=m1sq-m0sq-p1sq
      in=0
      iP=2
      do ep=-1,0
      A0(ep)=trI1(m1sq,musq,ep)
      B00(ep)=-0.5d0/dfloat(iP-2*in-1)*(
     . facnp(in,iP)*A0(ep)-f1*B1(ep)
     . +2d0*p1sq*B11(ep))
      enddo
   
   99 continue
      
      FB0(-2)=czip
      FB1(:,-2)=czip
      FB2(:,-2)=czip
       do ep=-1,0
      FB0(ep)=B0(ep)
      do n1=1,4
      FB1(n1,ep)=B1(ep)*p1(n1)
      do n2=n1,4
      FB2(y2(n1,n2),ep)=g(n1,n2)*B00(ep)+p1(n1)*p1(n2)*B11(ep)
      enddo
      enddo
      enddo

c--- store in cache
      do ep=-2,0
        FB0save(Nstore,ep)=FB0(ep)
        B00save(Nstore,ep)=B00(ep)
        do j=1,y1max
          FB1save(Nstore,j,ep)=FB1(j,ep)
        enddo
        do j=1,y2max
          FB2save(Nstore,j,ep)=FB2(j,ep)
        enddo
      enddo
      
      return 
      end

