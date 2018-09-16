      function cdot(e1,e2)
      implicit none
      include 'types.f'
      complex(dp):: cdot
      
      complex(dp):: e1(4),e2(4)
      cdot=e1(4)*e2(4)-e1(1)*e2(1)-e1(2)*e2(2)-e1(3)*e2(3)
      return
      end
