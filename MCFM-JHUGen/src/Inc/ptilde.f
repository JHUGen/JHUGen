      integer:: maxd,ndmax
!----maxd=The maximum possible number of dipoles
!----ndmax=The maximum number of dipoles for the problem at hand
      parameter (maxd=40)
      real(dp):: ptilde(0:maxd,mxpart,4)
      real(dp):: ptildejet(0:maxd,mxpart,4)
      common/ptildes/ptilde,ptildejet,ndmax
!$omp threadprivate(/ptildes/)
