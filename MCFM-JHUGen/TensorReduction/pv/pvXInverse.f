************************************************************************
* Inverse computes the inverse of a matrix.
* Input:
*   A: n-by-n matrix A
*   n: dimension of A
* Output:
*   A: mangled LU decomposition of A
*   Ainv: inverse of A
*   perm: permutation vector

      subroutine pvXInverse(A, Ainv, n, perm)
      implicit none
      integer n, perm(n)
      Double Complex A(n,n), Ainv(n,n)

      integer i, j

      call XLUDecomp(A, n, perm)
      do i = 1, n
        do j = 1, n
          Ainv(j,i) = dcmplx(0d0,0d0)
        enddo
        Ainv(i,i) = dcmplx(1d0,0d0)
        call XLUBackSubst(A, n, perm, Ainv(1,i))
      enddo
      end
