      double precision function hzgamwidth(mh)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      double precision mh,mhsq
      double complex f0DDHK
      
      mhsq=mh**2
      
      hzgamwidth=esq*Gf**2*wmass**2*xw/256d0/pi**5
     & *mh**3*(1d0-zmass**2/mhsq)**3*abs(f0DDHK(mhsq,zmass**2))**2
      
      return
      end
      
