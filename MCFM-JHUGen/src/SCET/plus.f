      function plus(pz,fx,L0,p1,fx1,L01,z,jaco)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: plus
      real(dp), intent(in)::pz,p1,fx,fx1,z,L0,L01,jaco
      if(abs(one-z) < 1.e-5_dp) then
      plus=L01*p1*fx1
      else
      plus=(pz*fx/z-p1*fx1)*L0*jaco+L01*p1*fx1
      endif
      return
      end
