      double complex function pvfndd(n,x,iep)
C----Implementation of DD Eq. 4.11
      implicit none
      include 'TRonshellcutoff.f'
      integer j,n,infty
      double complex xm1,x,cln,cone
      double precision iep
      parameter(cone=(1d0,0d0),infty=16) ! number of terms in sum
      
      xm1=x-cone
      if (abs(x) .lt. 10d0) then
        if (abs(x-cone) .lt. onshellcutoff) then
          pvfndd=0d0
        else
          pvfndd=(cone-x**(n+1))*(cln(xm1,iep)-cln(x,iep))
        endif
        do j=0,n
          pvfndd=pvfndd-x**(n-j)/dfloat(j+1)
        enddo
      else
        pvfndd=cln(cone-cone/x,iep)
        do j=n+1,n+infty
          pvfndd=pvfndd+x**(n-j)/dfloat(j+1)
        enddo
      endif
           
      return
      end
