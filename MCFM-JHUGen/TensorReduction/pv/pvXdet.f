************************************************************************
* Det computes the determinant of a matrix.
* Input:
*   A: n-by-n matrix A
*   n: dimension of A
* Output:
*   determinant of A
* Warning: A is overwritten

      Double Complex function pvXDet(A, n)
      implicit none
      integer n
        include 'pvNmax.f'
      Double Complex A(n,n)

      integer i, perm(Nmax)

      call XLUDecomp(A, n, perm)
      pvXDet = dcmplx(1d0,0d0)
      do i = 1, n
        pvXDet = pvXDet*A(i,i)
        if( perm(i) .ne. i ) pvXDet = -pvXDet
      enddo
      end

