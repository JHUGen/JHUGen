      subroutine pol_real(p,i,f)
      implicit none
      include 'constants.f'
      include 'swapxz.f'
c     massless vector polarization  subroutine
c     mcfm style momentum assignment
c     polarization only has real component
      integer i
      double complex p(4),f(4)
      double precision p0,px,py,pz,pv,ct,st,cphi,sphi
      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'pol_real:swapxz=',swapxz
      first=.false.
      if ((i .ne. 1) .and. (i .ne. -1)) then
      write(6,*) 'pol_real:pol out of bounds, i= ',i
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


      if (i .eq. 1) then
      if (swapxz) then
          f(4)=dcmplx(0d0,0d0)
          f(3)=+dcmplx(ct*cphi,0d0)
          f(2)=-dcmplx(ct*sphi,0d0)
          f(1)=+dcmplx(-st,0d0)
      else
          f(4)=dcmplx(0d0,0d0)
          f(1)=+dcmplx(ct*cphi,0d0)
          f(2)=+dcmplx(ct*sphi,0d0)
          f(3)=+dcmplx(-st,0d0)
      endif
      elseif (i .eq. -1) then
      if (swapxz) then
          f(4)=dcmplx(0d0,0d0)
          f(3)=-dcmplx(sphi,0d0)
          f(2)=-dcmplx(cphi,0d0)
          f(1)=+dcmplx(0d0,0d0)
      else
          f(4)=dcmplx(0d0,0d0)
          f(1)=-dcmplx(sphi,0d0)
          f(2)=+dcmplx(cphi,0d0)
          f(3)=+dcmplx(0d0,0d0)
      endif
      endif
      return
      end


