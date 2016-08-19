      subroutine qqb_dirgam_gvec(p,n,in,msq)
C*********************************************************************** 
c     Author: R.K. Ellis                                               *
c     October, 2002.                                                   *
c     Matrix element for gamma production                              *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p4)              *
c     q(-p1)+qbar(-p2)--> gamma(p3)+ g(p4)                             *
C*********************************************************************** 
      implicit none
      include 'constants.f'
      integer j,k,in
C--in is the label of the parton dotted with n
      double precision msq(-nf:nf,-nf:nf),msqa(-nf:nf,-nf:nf),
     . p(mxpart,4),n(4),nDn

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
 
      call qqb_dirgam(p,msqa)
      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2

      call checkndotp(p,n,in)

      do j=-nf,nf

      if (in .eq. 1) then
        msq(0,j)=-0.5d0*nDn*msqa(0,j)
      elseif (in .eq. 2) then
        msq(j,0)=-0.5d0*nDn*msqa(j,0)
      elseif (in .eq. 4) then      
        msq(j,-j)=-0.5d0*nDn*msqa(j,-j)
      endif

      enddo

      return
      end
