      subroutine qq4lggampf(i1,i2,i3,i4,i5,i6,i7,i8,
     & b7,b8,za,zb,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,
     & jdu,h1,h3,h5,h7,h8
      real(dp):: msq(2)
      complex(dp)::
     & a78(2,2,2,2,2,2),a87(2,2,2,2,2,2),aq(2,2,2,2,2,2),
     & amp78xy(2,2,2,2,2,2),amp87xy(2,2,2,2,2,2),
     & amp78yx(2,2,2,2,2,2),amp87yx(2,2,2,2,2,2),
     & tmp78(2,2,2,2,2,2),tmp87(2,2,2,2,2,2)

      call qq4lggamp(i1,i2,i3,i4,i5,i6,i7,i8,
     & b7,b8,za,zb,amp78xy,amp87xy)
      call qq4lggamp(i1,i2,i5,i6,i3,i4,i7,i8,
     & b7,b8,za,zb,amp78yx,amp87yx)
      tmp78=amp78yx
      tmp87=amp87yx
      do h3=1,2
      do h5=1,2
      amp78yx(:,:,h3,h5,:,:)=tmp78(:,:,h5,h3,:,:)
      amp87yx(:,:,h3,h5,:,:)=tmp87(:,:,h5,h3,:,:)
      enddo
      enddo
      

      a78(:,:,:,:,:,:)=amp78xy(:,:,:,:,:,:)+amp78yx(:,:,:,:,:,:)
      a87(:,:,:,:,:,:)=amp87xy(:,:,:,:,:,:)+amp87yx(:,:,:,:,:,:)
      aq(:,:,:,:,:,:)=a78(:,:,:,:,:,:)+a87(:,:,:,:,:,:)

      msq(:)=zip
      do jdu=1,2
      do h1=1,2
      do h3=1,2
      do h5=1,2
      do h7=1,2
      do h8=1,2
C23456789012345678901234567890123456789012345678901234567890123456789012
      msq(jdu)=msq(jdu)+gsq**2*esq**4*V*xn*(
     & +real(a78(jdu,h1,h3,h5,h7,h8)*conjg(a78(jdu,h1,h3,h5,h7,h8)))
     & +real(a87(jdu,h1,h3,h5,h7,h8)*conjg(a87(jdu,h1,h3,h5,h7,h8)))
     & -real(aq(jdu,h1,h3,h5,h7,h8)*conjg(aq(jdu,h1,h3,h5,h7,h8)))
     &  /xn**2)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
