      function fmt(s12,s34,s56)
      implicit none
      include 'types.f'
      complex(dp):: fmt
      include 'masses.f'
      real(dp):: s12,s34,s56
      fmt=((1._dp+(2._dp*s34+s12+s56)/15._dp/mt**2)/(24._dp*mt**2))
      return
      end

