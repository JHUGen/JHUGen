      subroutine gencol(x,xjac,xmin,emit,r)
      implicit none
      include 'types.f'
      
c---Generate an x value and store it for later retrieval
      integer:: emit,lemit
      real(dp):: x,xjac,xmin,xl,xljac,xlmin,r
      save xl,xljac,xlmin,lemit
      x=1._dp-(1._dp-xmin)*abs(1._dp-2._dp*r)
      xjac=2._dp*(1._dp-xmin)
      xl=x
      xljac=xjac
      xlmin=xmin
      lemit=emit
      return

      entry getcol(x,xjac,xmin,emit)
c---return the same values as last time
      x=xl
      xjac=xljac
      xmin=xlmin
      emit=lemit
      end
