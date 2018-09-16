c--- suite of dummy routines that replace OneLOop routines if
c--- that library is not linked
      subroutine error_olo_notlinked
      implicit none
      include 'types.f'
      
      write(6,*) 'OneLOop library has not been linked in'
      write(6,*) 'the MCFM makefile. Please recompile with'
      write(6,*) 'appropriate flag set.'  
      stop
          
      return
      end
      

      subroutine olo_unit(i1,aa)
      implicit none
      include 'types.f'
      
      integer:: i1
      character*(*) aa
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_onshell(x1)
      implicit none
      include 'types.f'
      
      real(dp):: x1
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_a0(c0,x1,x2)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_b0(c0,x1,x2,x3,x4)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2,x3,x4
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_c0(c0,x1,x2,x3,x4,x5,x6,x7)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2,x3,x4,x5,x6,x7
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_d0(c0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      
