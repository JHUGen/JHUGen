
      subroutine qqb_zaj_swap(p,msq)
      implicit none
      include 'types.f'
**************************************************************
*     Matrix element for                                     *
*     f(-p1)+f(-p2)-->Z(p3+p4)+gamma(p6)+g(p5)               *
*     this is basically qqb_zgam_g with p5 and p6 exchanged  *
*     Needed for Zgamgam dipoles                             *
**************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),pswap(mxpart,4)
      integer:: ii,jj
c-----swap p5 and p6
      do ii=1,4
         pswap(1,ii)=p(1,ii)
         pswap(2,ii)=p(2,ii)
         pswap(3,ii)=p(3,ii)
         pswap(4,ii)=p(4,ii)
         pswap(5,ii)=p(6,ii)
         pswap(6,ii)=p(5,ii)
      enddo
c-----initialize msq
      do ii=-nf,nf
      do jj=-nf,nf
         msq(ii,jj)=zip
      enddo
      enddo
c-----call qqb_zaj matelem
      call qqb_zaj(pswap,msq)
c-----done
      return
      end


