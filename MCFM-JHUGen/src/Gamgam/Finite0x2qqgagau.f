      function Finite0x2qqgagau(s,t)
!     Results taken from hep-ph/0201274
! \bibitem{Anastasiou:2002zn} 
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!  %``Two loop QED and QCD corrections to massless fermion boson scattering,''
!  Nucl.\ Phys.\ B {\bf 629}, 255 (2002)
!  [hep-ph/0201274].
      implicit none
      include 'types.f'
      real(dp)::Finite0x2qqgagau
      include 'constants.f'
      include 'nf.f'
      include 'scale.f'
      real(dp)::s,t,u,x,y,z,Lx,Ly,Lu,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x,
     & Li4,Li3,Li2,
     & AGTYAu,xAu,AGTYBu,xBu,AGTYD2u,xD2u,AGTYE3u,xE3u
      real(dp),parameter::sumq=11._dp/9._dp
      u=-s-t
      x=-t/s
      y=-u/s
      z=-y/x

      Lx=log(x)
      Ly=log(y)
      Lu=log(-u/musq)

      Li4x=Li4(x)
      Li4y=Li4(y)
      Li4z=Li4(z)

      Li3x=Li3(x)
      Li3y=Li3(y)

      Li2x=Li2(x)
      xAu=AGTYAu(s,t,u,Lx,Ly,Li2x,Li3x,Li4x,Li4y,Li3y,Li4z)
      xBu=AGTYBu(s,t,u,Lx,Ly,Lu,Li2x,Li3x,Li4x,Li3y,Li4y,Li4z)
      xD2u=AGTYD2u(s,t,u,Lx,Ly,Lu,Li2x,Li3x,Li4x,Li3y,Li4y,Li4z)
      xE3u=AGTYE3u(s,t,Lx,Ly,Lu)

      Finite0x2qqgagau=two*xn
     & *(sumq*TR*CF*xAu+CF**2*xBu+CF*CA*xD2u+nf*CF*xE3u)
      return
      end





