      function AGTYF2s(t,u,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.11
      include 'types.f'
      real(dp):: AGTYF2s,F2sx,t,u,Ls

      AGTYF2s=F2sx(t,u,Ls)
     &       +F2sx(u,t,Ls)

      return
      end


      function F2sx(t,u,Ls)
      implicit none
      include 'types.f'
      real(dp):: F2sx
      include 'zeta.f'
      real(dp):: t,u,Ls

      F2sx=(((-160*Ls)/27. + (32*Ls**2)/9. - (92*pisq)/27.)*t)/u

!     + \t \leftrightarrow u \
      return
      end
