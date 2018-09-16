      function CoupleC(j,k,nwz)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'ckm1.f'
      real(dp) CoupleC
      logical:: decay
      integer:: j,k,nwz
      complex(dp):: bwf
      common/bwf/bwf
      common/decay/decay
!$omp threadprivate(/bwf/,/decay/)

      decay=.true.

      if (decay) then
         if (nwz == 0) then
            if (j==k) then
            CoupleC=
     &      +(glsq(j,j)*flsq+grsq(j,j)*frsq)
     &      -2._dp*e(j,j)*(gl(j,j)*fl+gr(j,j)*fr)*real(bwf)
     &      +2._dp*e(j,j)**2*abs(bwf)**2
            else
            CoupleC=0._dp
            endif
         else
            CoupleC=Vsq(j,k)*glsq(j,k)*flsq
         endif
      else 
         if (nwz == 0) then
            if (j==k) then
            CoupleC=glsq(j,j)+grsq(j,j)
            else
            CoupleC=0._dp
            endif
         else
            CoupleC=Vsq(j,k)*glsq(j,k)
         endif
      endif

      return
      end

      function CoupleV(j,k,nwz)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'ckm1.f'
      real(dp):: CoupleV
      logical:: decay
      integer:: j,k,nwz
      complex(dp):: bwf
      common/bwf/bwf
      common/decay/decay
!$omp threadprivate(/bwf/,/decay/)
      decay=.true.
      if (decay) then
          if (nwz == 0) then
            if (j==k) then
            CoupleV=
     &      +(glsq(j,j)*frsq+grsq(j,j)*flsq)
     &      -2._dp*e(j,j)*(gl(j,j)*fr+gr(j,j)*fl)*real(bwf)
     &      +2._dp*e(j,j)**2*abs(bwf)**2
            else
            CoupleV=0._dp
            endif
          else
            CoupleV=Vsq(j,k)*glsq(j,k)*flsq
          endif
      else 
          if (nwz == 0) then
            if (j==k) then
            CoupleV=glsq(j,j)+grsq(j,j)
            else
            CoupleV=0._dp
            endif
          else
            CoupleV=Vsq(j,k)*glsq(j,k)
          endif
      endif

      return
      end

