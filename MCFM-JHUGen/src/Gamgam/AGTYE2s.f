      function AGTYE2s(s,t,u,Lx,Ly,Ls,Li2x,Li2y,Li3x,Li3y)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.6
      include 'types.f'
      real(dp):: AGTYE2s,E2sx
      real(dp):: s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li2y,Li3y
      AGTYE2s=E2sx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li3y)
     &       +E2sx(s,u,t,Ly,Lx,Ls,Li2y,Li3y,Li3x)

      return
      end


      function E2sx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li3y)
      implicit none
      include 'types.f'
      real(dp):: E2sx
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li3y


      E2sx= (-4._dp/3._dp*Li3x-4._dp/
     & 3*Li3y+ (-4._dp/3._dp*Ly+4._dp/3._dp*Lx
     & )*Li2x
     & -11._dp/ 18*Lx**3+ (-1._dp/ 2*Ly+5._dp/
     & 6-Ls )*Lx**2
     & + (-5._dp/3._dp*Ly**2+ (10._dp/ 9-2*Ls
     & )*Ly-1._dp/ 9*pisq+37._dp/ 9-31._dp/ 6*Ls
     & )*Lx
     & -41._dp/ 18*Ly**2+ (5._dp/ 9*pisq+43._dp/ 9-31._dp/
     & 6*Ls )*Ly+65._dp/ 81+19._dp/
     & 18*pisq*Ls-11._dp/3._dp*Ls**2
     & -13._dp/ 9*zeta3+206._dp/
     & 27*Ls-275._dp/ 108*pisq )*(t/u)
     & +(-Lx**2+ (2*Ly+2._dp)*Lx-Ly**2-2*Ly-pisq)*(t/s)**2
     & + (14._dp/ 9*Lx**2+ (2._dp/3._dp*Ly-1._dp/
     & 3*Ly**2+38._dp/ 9 )*Lx-11._dp/
     & 9*Ly**3-2*Ly**2*Ls
     & -17._dp/ 18*pisq+ (-8._dp/
     & 9*pisq-2*Ls )*Ly )
! + \t \leftrightarrow u \
      return
      end
