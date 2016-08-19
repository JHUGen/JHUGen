      double precision function del4(p1,p2,p3,p4,k1,k2,k3,k4)
      implicit none
      double precision ovdot,del3,p1(4),p2(4),p3(4),p4(4),
     &                            k1(4),k2(4),k3(4),k4(4)
      del4=ovdot(p1,k1)*del3(p2,p3,p4,k2,k3,k4)
     &    -ovdot(p1,k2)*del3(p2,p3,p4,k1,k3,k4)
     &    +ovdot(p1,k3)*del3(p2,p3,p4,k1,k2,k4)
     &    -ovdot(p1,k4)*del3(p2,p3,p4,k1,k2,k3)
      return
      end
      
