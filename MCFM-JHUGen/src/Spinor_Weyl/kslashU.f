c-----multiplication of k-slash from the left with a spinor sp
C-----and return resultant spinor f. Weyl representation
      subroutine kslashU(k,sp,f) 
      implicit none
      include 'swapxz.f'
      double complex sp(4),k(4),f(4),czip,im,kslash(4,4),E,kx,ky,kz
      integer i,j
      parameter(czip=(0d0,0d0),im=(0d0,1d0))

      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'kslashU:swapxz=',swapxz
      first=.false.
      endif

      if (swapxz) then
C----create kslash after performing the swap (x<->z),(y->-y)
      E=k(4)
      kx=+k(3)
      ky=-k(2)
      kz=+k(1)
      else
      E=k(4)
      kx=+k(1)
      ky=+k(2)
      kz=+k(3)
      endif


      kslash(1,1)=czip
      kslash(1,2)=czip
      kslash(1,3)=E+kz
      kslash(1,4)=kx-im*ky

      kslash(2,1)=czip
      kslash(2,2)=czip
      kslash(2,3)=kx+im*ky
      kslash(2,4)=E-kz

      kslash(3,1)=E-kz
      kslash(3,2)=-kx+im*ky
      kslash(3,3)=czip
      kslash(3,4)=czip

      kslash(4,1)=-kx-im*ky
      kslash(4,2)=E+kz
      kslash(4,3)=czip
      kslash(4,4)=czip

      do i=1,4
      f(i)=czip
      do j=1,4
      f(i)=f(i)+kslash(i,j)*sp(j)
      enddo
      enddo  
      return
      end



