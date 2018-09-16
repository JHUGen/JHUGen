      subroutine qqWZggampf(i1,i2,i3,i4,i5,i6,i7,i8,
     & b7,b8,za,zb,msq)
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,h5,h7,h8
      double precision msq
      double complex a78(2,2,2),a87(2,2,2),aq(2,2,2)

      call qqWZggamp(i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,za,zb,a78,a87)
      aq(:,:,:)=a78(:,:,:)+a87(:,:,:)

      msq=zip
      do h5=1,2
      do h7=1,2
      do h8=1,2
C23456789012345678901234567890123456789012345678901234567890123456789012
      msq=msq+gsq**2*esq**4*V*xn*(
     & +dble(a78(h5,h7,h8)*Dconjg(a78(h5,h7,h8)))
     & +dble(a87(h5,h7,h8)*Dconjg(a87(h5,h7,h8)))
     & -dble(aq(h5,h7,h8)*Dconjg(aq(h5,h7,h8)))
     &  /xn**2)
      enddo
      enddo
      enddo
      return
      end
