      function hzgamwidth(mh)
      implicit none
      include 'types.f'
      real(dp):: hzgamwidth
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      real(dp):: mh,mhsq
      complex(dp):: f0DDHK
      
      mhsq=mh**2
      
      hzgamwidth=esq*Gf**2*wmass**2*xw/256._dp/pi**5
     & *mh**3*(1._dp-zmass**2/mhsq)**3*abs(f0DDHK(mhsq,zmass**2))**2
      
      return
      end
      
