      function spstrng0(Ub,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng0

      complex(dp)::Ub(4),V(4)
      spstrng0=Ub(1)*V(1)+Ub(2)*V(2)+Ub(3)*V(3)+Ub(4)*V(4)
      return
      end

      function spstrng1(Ub,zk1,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng1

      complex(dp)::Ub(4),V(4),zk1(4),temp(4),spstrng0
C -   contract Ubm*Gamma(zk1)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng1=spstrng0(temp,V)
      return
      end

      function spstrng2(Ub,zk1,zk2,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng2

      complex(dp)::Ub(4),V(4),zk1(4),zk2(4),temp(4),spstrng1
C -   contract Ubm*Gamma(zk1)*Gamma(zk2)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng2=spstrng1(temp,zk2,V)
      return
      end


      function spstrng3(Ub,zk1,zk2,zk3,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng3
      complex(dp)::Ub(4),V(4),zk1(4),zk2(4),zk3(4),temp(4),spstrng2
C -   contract Ubm*Gamma(zk1)*Gamma(zk2)*Gamma(zk3)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng3=spstrng2(temp,zk2,zk3,V)
      return
      end

      function spstrng4(Ub,zk1,zk2,zk3,zk4,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng4

      complex(dp)::Ub(4),V(4),zk1(4),zk2(4),zk3(4),zk4(4),temp(4),
     & spstrng3
C -   contract Ubm*Gamma(zk1)*Gamma(zk2)*Gamma(zk3)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng4=spstrng3(temp,zk2,zk3,zk4,V)
      return
      end

      function spstrng5(Ub,zk1,zk2,zk3,zk4,zk5,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng5

      complex(dp)::Ub(4),V(4),zk1(4),zk2(4),zk3(4),zk4(4),zk5(4),
     & temp(4),spstrng4
C -   contract Ubm*Gamma(zk1)*Gamma(zk2)*Gamma(zk3)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng5=spstrng4(temp,zk2,zk3,zk4,zk5,V)
      return
      end

      function spstrng6(Ub,zk1,zk2,zk3,zk4,zk5,zk6,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng6

      complex(dp)::Ub(4),V(4),zk1(4),zk2(4),zk3(4),zk4(4),zk5(4),
     & zk6(4),temp(4),spstrng5
C -   contract Ubm*Gamma(zk1)*Gamma(zk2)*Gamma(zk3)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng6=spstrng5(temp,zk2,zk3,zk4,zk5,zk6,V)
      return
      end


      function spstrng7(Ub,zk1,zk2,zk3,zk4,zk5,zk6,zk7,V)
      implicit none
      include 'types.f'
      complex(dp)::spstrng7

      complex(dp)::Ub(4),V(4),zk1(4),zk2(4),zk3(4),zk4(4),zk5(4),
     & zk6(4),zk7(4),temp(4),spstrng6
C -   contract Ubm*Gamma(zk1)*Gamma(zk2)*Gamma(zk3)*V
      call Ubkslash(Ub,zk1,temp)
      spstrng7=spstrng6(temp,zk2,zk3,zk4,zk5,zk6,zk7,V)
      return
      end


