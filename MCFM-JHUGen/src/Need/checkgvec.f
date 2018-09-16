      subroutine checkgvec(j,k,ip,p,routine_lo,routine_gvec)
c--- routine that checks the squared LO matrix elements (calculated
c--- using routine_lo) against the squared "gvec" matrix elements
c--- (calculated using routine_gvec), which are spin-summed
c--- by using particular forms for the polarization vectors ("n")
      implicit none
      include 'constants.f'
      integer j,k,ip
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),
     . msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),n(4)
      external routine_lo,routine_gvec
      
c--- calculate LO
      call routine_lo(p,msq)
      
c--- calculate gvec routine with first polarization vector
      if (ip .le. 2) then
        n(4)=0d0
        n(1)=1d0
        n(2)=0d0
        n(3)=0d0
      else
        n(1)=p(ip,2)/dsqrt(p(ip,1)**2+p(ip,2)**2)
        n(2)=-p(ip,1)/dsqrt(p(ip,1)**2+p(ip,2)**2)
        n(3)=0d0
        n(4)=0d0       
      endif 
      call routine_gvec(p,n,ip,msq1)
      
c--- calculate gvec routine with second polarization vector
      if (ip .le. 2) then
        n(4)=0d0
        n(1)=0d0
        n(2)=1d0
        n(3)=0d0
      else
        n(1)=p(ip,1)*p(ip,3)/p(ip,4)/dsqrt(p(ip,1)**2+p(ip,2)**2)
        n(2)=p(ip,2)*p(ip,3)/p(ip,4)/dsqrt(p(ip,1)**2+p(ip,2)**2)
        n(3)=-dsqrt(p(ip,1)**2+p(ip,2)**2)/p(ip,4)
        n(4)=0d0       
      endif 
      call routine_gvec(p,n,ip,msq2)
      
      write(6,*) 'initial state: ',j,k
      write(6,*) 'gvec 1',msq1(j,k)
      write(6,*) 'gvec 2',msq2(j,k)
      write(6,*) 'gvec sum',msq1(j,k)+msq2(j,k)
      write(6,*) 'l. order',msq(j,k)
      write(6,*) '   RATIO',msq(j,k)/(msq1(j,k)+msq2(j,k))
      pause
      
      return
      end
      
