      subroutine kslashU(zk,spinor,f) 
      implicit none
      include 'types.f'
      
c-----multiplication of k-slash from the left with a spinor spinor
C-----and return resultant spinor f, Majorana representation
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'swapxz.f'
      complex(dp):: spinor(4),zk(4),f(4)
      real(dp):: spr(4),spi(4),kslash(4,4),fr(4),fi(4),Ek,kx,ky,kz
      integer:: j
      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'kslashU:swapxz=',swapxz
      first=.false.
      endif

      if (swapxz) then
C----create kslash after performing the swap (x<->z),(y->-y)
      Ek=+real(zk(4))
      kx=+real(zk(3))
      ky=-real(zk(2))
      kz=+real(zk(1))
      else
      Ek=+real(zk(4))
      kx=+real(zk(1))
      ky=+real(zk(2))
      kz=+real(zk(3))
      endif

C---kslash in Majorana representation with factor of i removed
      kslash(1,1)=-kx
      kslash(1,2)=+kz
      kslash(1,3)=zip
      kslash(1,4)=-ky-Ek

      kslash(2,1)=+kz
      kslash(2,2)=+kx
      kslash(2,3)=+ky+Ek
      kslash(2,4)=zip

      kslash(3,1)=zip
      kslash(3,2)=+ky-Ek
      kslash(3,3)=-kx
      kslash(3,4)=+kz

      kslash(4,1)=-ky+Ek
      kslash(4,2)=zip
      kslash(4,3)=+kz
      kslash(4,4)=+kx

      spr(:)=real(spinor(:))
      spi(:)=aimag(spinor(:))
      do j=1,4
      fr(j)=kslash(j,1)*spr(1)+kslash(j,2)*spr(2)
     &     +kslash(j,3)*spr(3)+kslash(j,4)*spr(4)
      fi(j)=kslash(j,1)*spi(1)+kslash(j,2)*spi(2)
     &     +kslash(j,3)*spi(3)+kslash(j,4)*spi(4)
C---recompose complex f restoring factor of i
      f(j)=cplx2(-fi(j),+fr(j))
      enddo  
      return
      end



