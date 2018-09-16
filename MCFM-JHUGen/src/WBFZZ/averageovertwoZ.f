      subroutine averageovertwoZ(fxn,p,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      integer:: kcount,jcount
      real(dp):: msq(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),
     & p(mxpart,4),pout1(mxpart,4),pout(mxpart,4),statfac
      pout(:,:)=0d0
      msq(:,:)=0d0
C---sum over 36 values;
      do jcount=1,6
      call pgen(p,jcount,3,4,pout1)
      do kcount=1,6
      call pgen(pout1,kcount,5,6,pout)
      call fxn(pout,msq1)
      msq(:,:)=msq(:,:)+msq1(:,:)
      enddo
      enddo
      statfac=0.5d0
      msq(:,:)=msq(:,:)*(zwidth**2/(4d0*esq*(le**2+re**2)))**2*statfac
      return
      end
