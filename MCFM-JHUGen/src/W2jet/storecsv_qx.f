      subroutine storecsv_qx(i,j)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the W2jet_gvec matrix elements into elements (..,i,j) of q1q2
      
      include 'mmsqv_cs.f'
      integer:: i,j,k
      real(dp):: q1q2(0:2,-1:1,-1:1)
      common/q1q2/q1q2
!$omp threadprivate(/q1q2/)
      
      do k=0,2
        q1q2(k,i,j)=mmsqv_cs(k,+1,+1)
      enddo
      
      return
      end
