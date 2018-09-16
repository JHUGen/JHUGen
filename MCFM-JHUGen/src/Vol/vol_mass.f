      function vol_mass(decaymass,s)
      implicit none
      include 'types.f'
      real(dp):: vol_mass
      
      real(dp):: decaymass,b0,bp,bm,s,vol,ro,ddilog
      
      ro=4._dp*decaymass**2/s
      b0=sqrt(1._dp-ro)
      bp=0.5_dp*(1._dp+b0)      
      bm=0.5_dp*(1._dp-b0)      

c      temp=(4._dp*b0*(27._dp - 26._dp*b0**2 + 3._dp*b0**4)
c     & + 3._dp*(-1._dp + b0**2)**2*log((-1._dp + b0)**2)**2 + 
c     & 6._dp*(-1._dp + b0**2)**2*log(1._dp-b0)*(1._dp-4._dp*log(1._dp + b0))+ 
c     & 6._dp*(-1._dp + b0**2)*log(1._dp + b0)*
c     & (9._dp - b0**4 + (-1._dp + b0**2)*log(16._dp)
c     & - 2._dp*(-1._dp + b0**2)*log(1._dp + b0)) 
c     & + 3._dp*(-1._dp + b0**2)*log((-1._dp + b0)**2)*
c     & (-8._dp + b0**4 + log(16._dp) - b0**2*(1._dp + log(16._dp)) + 
c     & 4._dp*(-1._dp + b0**2)*log(1._dp + b0)) + 
c     & 24._dp*(-1._dp + b0**2)**2*Ddilog((1._dp - b0)/2._dp) - 
c     & 24._dp*(-1._dp + b0**2)**2*Ddilog((1._dp + b0)/2._dp))/192._dp
c     
c--- normalize to massless case (1/12)
c      temp=temp/(1._dp/12._dp)
 
      vol_mass=b0*(1._dp+5._dp*ro+0.75_dp*ro**2)
     & +1.5_dp*ro**2*(Ddilog(bm)-Ddilog(bp))
     & -ro*log(bp/bm)*(24._dp+6._dp*ro-3._dp*ro**2+6._dp*ro*log(bp*bm))/8._dp
      vol_mass=vol_mass*vol(s,4)
      
      return
      end
      
