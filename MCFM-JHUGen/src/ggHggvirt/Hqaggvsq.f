      function Hqaggvsq(j1,j2,j3,j4)
      implicit none
      include 'types.f'
      real(dp):: Hqaggvsq
      
      include 'epinv.f'
      include 'epinv2.f'
      include 'first_time.f'
C---  matrix element squared for H--q(j1)+a(j2)+g(j3)+g(j4)
      integer:: al,j1,j2,j3,j4
      real(dp):: q(4,5),qswap(4,5),sqres(-2:0)
      common/GZmom/q
!$omp threadprivate(/GZmom/)

      do al=1,4
      qswap(al,1)=q(al,j1)
      qswap(al,2)=q(al,j2)
      qswap(al,3)=q(al,j3)
      qswap(al,4)=q(al,j4)
      qswap(al,5)=-q(al,j1)-q(al,j2)-q(al,j3)-q(al,j4)
      enddo

C      write(*,*) 'scale, hmass',scale, hmass
      call GZHqaggvsqPoles(qswap,sqres)
c      call GZHqaggvsq(qswap,sqres,first_time,new_event,scale,hmass) 
      Hqaggvsq=epinv*epinv2*sqres(-2)+epinv*sqres(-1)+sqres(0)
      if (first_time) first_time = .false. 
      !if (new_event) new_event = .false. 
      return
      end

