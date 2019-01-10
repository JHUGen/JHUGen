      subroutine checkndotp(p,n,ig)
c--- routine to check that the vector "n" contracted with the LO
c--- matrix elements, is properly defined such that n.p(ig)=0
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),n(4),test,tolerance
      integer ig
      parameter (tolerance=1d-4)
      
      test=n(4)*p(ig,4)-n(1)*p(ig,1)-n(2)*p(ig,2)-n(3)*p(ig,3)
      test=test/max(abs(n(4)),abs(n(1)),abs(n(2)),abs(n(3)))
     .         /max(abs(p(ig,4)),abs(p(ig,1)),abs(p(ig,2)),abs(p(ig,3)))
      
      if (test .gt. tolerance) then
        write(6,*) 'Tolerance for n.p check in gvec routine exceeded!'
        write(*,*) 'tolerance: ',tolerance
        write(6,*) 'test: ',test
        write(*,*) 'index of gluon: ',ig
        write(6,*) 'n: ',n(1),n(2),n(3),n(4)
        write(6,*) 'p: ',p(ig,1),p(ig,2),p(ig,3),p(ig,4)
        write(6,*)
        call writeout(p)
        call flush(6)
        stop
      endif
      
      return
      end
      
