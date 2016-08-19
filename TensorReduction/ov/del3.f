      double precision function del3(p1,p2,p3,k1,k2,k3)
      implicit none
      double precision ovdot,p1(4),p2(4),p3(4),k1(4),k2(4),k3(4),del2
      del3=ovdot(p1,k1)*del2(p2,p3,k2,k3)
     &    -ovdot(p1,k2)*del2(p2,p3,k1,k3)
     &    +ovdot(p1,k3)*del2(p2,p3,k1,k2)
      return
      end
      
