      subroutine arraysort(imax,vals,sortorder)
      implicit none
      include 'types.f'
c--- Given a a 1-dimensional array of quantities, vals(1:mxpart),
c--- this routine sorts the first imax entries such that
c--- vals(sortorder(1)) > vals(sortorder(2)) >...> vals(sortorder(imax)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: imax,sortorder(mxpart),i1,i2,iswap
      real(dp):: vals(mxpart)

c-- initial ordering      
      do i1=1,imax
      sortorder(i1)=i1
      enddo

c--- catch trivial case      
      if (imax == 1) return
      
      do i1=imax-1,1,-1
      do i2=1,i1
        if (vals(sortorder(i2)) < vals(sortorder(i2+1))) then
          iswap=sortorder(i2)
          sortorder(i2)=sortorder(i2+1)
          sortorder(i2+1)=iswap
        endif
      enddo
      enddo
      
      return
      end
      
