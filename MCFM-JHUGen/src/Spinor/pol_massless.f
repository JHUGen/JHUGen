      subroutine pol_massless(p,i,f)
      implicit none
      include 'constants.f'
      include 'swapxz.f'
c     massless vector polarization  subroutine
c     mcfm style momentum assignment
      integer i,pol
      double complex p(4),f(4)
      double precision p0,px,py,pz,pv,ct,st,cphi,sphi
      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'pol_massless:swapxz=',swapxz
      first=.false.
      if ((i .ne. 1) .and. (i .ne. -1)) then
      write(6,*) 'pol_massless:pol out of bounds, i= ',i
      stop
      endif
      endif

      if (swapxz) then
      p0=dreal(p(4))
      px=+dreal(p(3))
      py=-dreal(p(2))
      pz=+dreal(p(1))
      else
      p0=dreal(p(4))
      px=dreal(p(1))
      py=dreal(p(2))
      pz=dreal(p(3))
      endif
      pv=dsqrt(dabs(p0**2))
      ct=pz/pv
      st=dsqrt(dabs(1d0-ct**2))

      if (st.lt.1d-8) then 
        cphi=1d0
        sphi=0d0
      else
        cphi= px/pv/st
        sphi= py/pv/st
      endif


c     the following ifstatement distinguishes between 
c     positive and negative energies
      if ( p0.gt.0d0) then  
        pol=i
      else
        pol=-i
      endif

      if (swapxz) then
      f(4)=dcmplx(0d0,0d0)
      f(3)=+dcmplx(ct*cphi/rt2,+pol*sphi/rt2)
      f(2)=-dcmplx(ct*sphi/rt2,-pol*cphi/rt2)
      f(1)=+dcmplx(-st/rt2,0d0)
      else
      f(4)=dcmplx(0d0,0d0)
      f(1)=+dcmplx(ct*cphi/rt2,+pol*sphi/rt2)
      f(2)=+dcmplx(ct*sphi/rt2,-pol*cphi/rt2)
      f(3)=+dcmplx(-st/rt2,0d0)
      endif      
      return
      end


