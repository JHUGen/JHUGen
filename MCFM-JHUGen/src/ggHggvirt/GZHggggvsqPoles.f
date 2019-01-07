      subroutine GZHggggvsqPoles(p,sqres)
      implicit none 
      include 'constants.f' 
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f' 
      include 'scale.f' 
      include 'b0.f'
      double precision   p(4,5),sqres(-2:0)
      double precision   Hgggg,Hgggg_1234,Hgggg_1243,Hgggg_1423
      double precision   q(mxpart,4)
      integer i,j
      double precision   s12,s13,s14,s23,s24,s34
      double complex lnrat 


C --  reversal of GZ notation and restore dimension of the momenta 
      do i=1,4
         do j=1,4
            q(i,j) = p(j,i)*hmass
         enddo
      enddo
      call spinoru(4,q,za,zb)

      call h4g(1,2,3,4,Hgggg,Hgggg_1234,Hgggg_1243,Hgggg_1423)

      s12=s(1,2)
      s13=s(1,3)
      s14=s(1,4)
      s23=s(2,3)
      s24=s(2,4)
      s34=s(3,4)
         
c--- NB: removed factor of ason2pi which is restored
c--- in the wrapping routine
      sqres(-2) = -4d0*xn*Hgggg
      sqres(-1) = xn*(
     . dble(lnrat(-s12,musq)+lnrat(-s23,musq)
     .     +lnrat(-s34,musq)+lnrat(-s14,musq))*Hgggg_1234+

     . dble(lnrat(-s12,musq)+lnrat(-s24,musq)
     .     +lnrat(-s34,musq)+lnrat(-s13,musq))*Hgggg_1243+

     . dble(lnrat(-s14,musq)+lnrat(-s24,musq)
     .     +lnrat(-s23,musq)+lnrat(-s13,musq))*Hgggg_1423)

      sqres(0)  = 11d0*Hgggg
c--- The finite part which will be calculated using the semi-numerical
c--- procedure has to be added to this piece, which comes from the 
c--- finite renormalization of the effective Lagrangian

c--- Add on renormalization factor
      sqres(-1)=sqres(-1)-4d0*b0*Hgggg

C       write(*,*) 'Poles Hggggv K: Born', Hgggg
C      write(*,*) 'Poles Hggggv K: 1/eps^2', sqres(-2) 
C      write(*,*) 'Poles Hggggv K: 1/eps^1', sqres(-1) 

      return
      end
