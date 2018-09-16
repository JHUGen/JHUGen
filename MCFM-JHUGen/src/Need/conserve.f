      subroutine conserve(p)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),dot
      integer:: nu

      do nu=1,4
      write(6,*) nu,p(1,nu)+p(2,nu)+p(3,nu)+p(4,nu)+p(5,nu)+p(6,nu)
     & +p(7,nu)   
      enddo
      write(6,*) 'dot',1,1,dot(p,1,1),p(1,4),p(1,3),p(1,2),p(1,1)
      write(6,*) 'dot',2,2,dot(p,2,2),p(2,4),p(2,3),p(2,2),p(2,1)
      write(6,*) 'dot',3,3,dot(p,3,3),p(4,4),p(4,3),p(4,2),p(4,1)
      write(6,*) 'dot',4,4,dot(p,4,4),p(3,4),p(3,3),p(3,2),p(3,1)
      write(6,*) 'dot',5,5,dot(p,5,5),p(5,4),p(5,3),p(5,2),p(5,1)
      write(6,*) 'dot',6,6,dot(p,6,6),p(6,4),p(6,3),p(6,2),p(6,1)
      write(6,*) 'dot',7,7,dot(p,7,7),p(7,4),p(7,3),p(7,2),p(7,1)
       
      write(6,*) 
   
      pause
      return 
      end

      subroutine conserve5(p)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),dot
      integer:: nu
      write(6,*) 
      do nu=1,4
      write(6,*) 'sum',
     & p(1,nu)+p(2,nu)+p(3,nu)+p(4,nu)+p(5,nu)+p(6,nu)+p(7,nu)   
      enddo
      write(6,*) 'dot',1,1,dot(p,1,1)
      write(6,*) 'dot',2,2,dot(p,2,2)
      write(6,*) 'dot',4,4,dot(p,4,4)
      write(6,*) 'dot',3,3,dot(p,3,3)
      write(6,*) 'dot',5,5,dot(p,5,5)
      write(6,*) 'dot',6,6,dot(p,6,6)
      write(6,*) 'dot',7,7,dot(p,7,7)
      write(6,*) 

      pause
      return 
      end

      subroutine conserve8(p)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4)
      integer:: nu
      write(6,*) 
      do nu=1,4
      write(6,*) 'sum',
     & p(1,nu)+p(2,nu)+p(3,nu)+p(4,nu)+p(5,nu)+p(6,nu)+p(7,nu)+p(8,nu)   
      enddo
      pause
      return 
      end
