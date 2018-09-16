      subroutine pol_real(p,i,f)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'swapxz.f'
c     massless vector polarization  subroutine
c     mcfm style momentum assignment
c     polarization only has real component
      integer:: i
      complex(dp):: p(4),f(4)
      real(dp):: p0,px,py,pz,pv,ct,st,cphi,sphi
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
      p0=real(p(4))
      px=+real(p(3))
      py=-real(p(2))
      pz=+real(p(1))
      else
      p0=real(p(4))
      px=real(p(1))
      py=real(p(2))
      pz=real(p(3))
      endif

      pv=sqrt(abs(p0**2))
      ct=pz/pv
      st=sqrt(abs(1._dp-ct**2))

      if (st<1.e-8_dp) then 
        cphi=1._dp
        sphi=zip
      else
        cphi= px/pv/st
        sphi= py/pv/st
      endif


      if (i == 1) then
      if (swapxz) then
          f(4)=cplx2(zip,zip)
          f(3)=+cplx2(ct*cphi,zip)
          f(2)=-cplx2(ct*sphi,zip)
          f(1)=+cplx2(-st,zip)
      else
          f(4)=cplx2(zip,zip)
          f(1)=+cplx2(ct*cphi,zip)
          f(2)=+cplx2(ct*sphi,zip)
          f(3)=+cplx2(-st,zip)
      endif
      elseif (i == -1) then
      if (swapxz) then
          f(4)=cplx2(zip,zip)
          f(3)=-cplx2(sphi,zip)
          f(2)=-cplx2(cphi,zip)
          f(1)=+cplx2(zip,zip)
      else
          f(4)=cplx2(zip,zip)
          f(1)=-cplx2(sphi,zip)
          f(2)=+cplx2(cphi,zip)
          f(3)=+cplx2(zip,zip)
      endif
      endif
      return
      end


