      subroutine GZHqaggvsqPoles(p,sqres)
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f' 
      include 'scale.f' 
      include 'b0.f'
      real(dp)::   p(4,5),sqres(-2:0)
      real(dp)::   Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym
      real(dp)::   q(mxpart,4)
      integer:: i,j
      real(dp)::   s12,s13,s14,s23,s24,s34
      complex(dp):: lnrat 


C --  reversal of GZ notation and restore dimension of the momenta 
      do i=1,4
         do j=1,4
            q(i,j) = p(j,i)*hmass
         enddo
      enddo
      call spinoru(4,q,za,zb)

      call hqqgg(1,2,3,4,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)

      s12=s(1,2)
      s13=s(1,3)
      s14=s(1,4)
      s23=s(2,3)
      s24=s(2,4)
      s34=s(3,4)
               
c--- NB: removed factor of ason2pi which is restored
c--- in the wrapping routine
      sqres(-2) = (-3._dp*xn+one/xn)*Hqagg
      sqres(-1) = (
     &  (-3._dp*Cf+(11._dp*xn/3._dp-2._dp/3._dp*nf-1._dp/xn*real(lnrat(-s12,musq))))
     &  *Hqagg
     &  +(xn*real(lnrat(-s13,musq)+lnrat(-s34,musq)+lnrat(-s24,musq)))
     &  *Hqagg_ab 
     &  +(xn*real(lnrat(-s14,musq)+lnrat(-s34,musq)+lnrat(-s23,musq)))
     &  *Hqagg_ba 
     &  +(xn*real(-lnrat(-s12,musq)+lnrat(-s13,musq)+lnrat(-s23,musq)
     &       +lnrat(-s14,musq)+lnrat(-s24,musq)))*Hqagg_sym)


      sqres(0)  = 11._dp*Hqagg
c--- The finite part which will be calculated using the semi-numerical
c--- procedure has to be added to this piece, which comes from the 
c--- finite renormalization of the effective Lagrangian

c--- Add on renormalization factor
      sqres(-1)=sqres(-1)-4._dp*b0*Hqagg

C      write(*,*) 'Poles qagggv K: 1/eps^2', sqres(-2) 
C      write(*,*) 'Poles qagggv K: 1/eps^1', sqres(-1) 
      return
      end
