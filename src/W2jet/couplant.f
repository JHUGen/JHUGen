      real*8 function CoupleC(j,k,nwz)
      implicit none
      include 'constants.f'
      include 'ckm.f'
      include 'ckm1.f'
      logical decay
      integer j,k,nwz
      double complex bwf
      common/bwf/bwf
      common/decay/decay
!$omp threadprivate(/bwf/,/decay/)

      decay=.true.

      if (decay) then
         if (nwz .eq. 0) then
            if (j.eq.k) then
            CoupleC=
     &      +(glsq(j,j)*flsq+grsq(j,j)*frsq)
     &      -2d0*e(j,j)*(gl(j,j)*fl+gr(j,j)*fr)*dble(bwf)
     &      +2d0*e(j,j)**2*abs(bwf)**2
            else
            CoupleC=0d0
            endif
         else
            CoupleC=Vsq(j,k)*glsq(j,k)*flsq
         endif
      else 
         if (nwz .eq. 0) then
            if (j.eq.k) then
            CoupleC=glsq(j,j)+grsq(j,j)
            else
            CoupleC=0d0
            endif
         else
            CoupleC=Vsq(j,k)*glsq(j,k)
         endif
      endif

      return
      end

      real*8 function CoupleV(j,k,nwz)
      implicit none
      include 'constants.f'
      include 'ckm.f'
      include 'ckm1.f'
      logical decay
      integer j,k,nwz
      double complex bwf
      common/bwf/bwf
      common/decay/decay
!$omp threadprivate(/bwf/,/decay/)
      decay=.true.
      if (decay) then
          if (nwz .eq. 0) then
            if (j.eq.k) then
            CoupleV=
     &      +(glsq(j,j)*frsq+grsq(j,j)*flsq)
     &      -2d0*e(j,j)*(gl(j,j)*fr+gr(j,j)*fl)*dble(bwf)
     &      +2d0*e(j,j)**2*abs(bwf)**2
            else
            CoupleV=0d0
            endif
          else
            CoupleV=Vsq(j,k)*glsq(j,k)*flsq
          endif
      else 
          if (nwz .eq. 0) then
            if (j.eq.k) then
            CoupleV=glsq(j,j)+grsq(j,j)
            else
            CoupleV=0d0
            endif
          else
            CoupleV=Vsq(j,k)*glsq(j,k)
          endif
      endif

      return
      end

