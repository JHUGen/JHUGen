      function cdotpr(J1,J2)
      implicit none
      include 'types.f'
      complex(dp):: cdotpr
      
      complex(dp):: J1(4),J2(4)
      cdotpr=J1(4)*J2(4)-J1(1)*J2(1)-J1(2)*J2(2)-J1(3)*J2(3)
      return
      end
