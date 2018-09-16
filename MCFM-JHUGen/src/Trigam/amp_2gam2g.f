      function amp_2gam2g_mppppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_mppppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      -       +         +         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=2, m=2 and labels permuted appropriately
c--- negative helicity gluon i=1
      amp_2gam2g_mppppm=      
     & za(j5,j1)*za(j6,j1)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_2gam2g_pmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_pmpppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +       -         +         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=2, m=2 and labels permuted appropriately
c--- negative helicity gluon i=2
      amp_2gam2g_pmpppm=      
     & za(j5,j2)*za(j6,j2)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_2gam2g_ppmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_ppmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +       +         -         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=2, m=2 and labels permuted appropriately
c--- negative helicity photon i=3
      amp_2gam2g_ppmppm=      
     & za(j5,j3)*za(j6,j3)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_2gam2g_pppmpm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_pppmpm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +       +         +         -
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=2, m=2 and labels permuted appropriately
c--- negative helicity photon i=4
      amp_2gam2g_pppmpm=      
     & za(j5,j4)*za(j6,j4)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_2gam2g_mmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_mmpppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      -       -         +         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
      
c--- Equation (4.1)
c--- note: I have interpreted, e.g. <5|1+2|6> as [5|1+2|6> in order
c--- to get the correct <| and |] required by the helicities
      amp_2gam2g_mmpppm=      
     &+(zb(j5,j1)*za(j1,j6)+zb(j5,j2)*za(j2,j6))*za(j5,j6)*t(j1,j2,j6)
     & /(zb(j5,j2)*zb(j1,j6)*zb(j2,j1)*za(j6,j3)
     &  *za(j6,j4)*za(j5,j3)*za(j5,j4))
     &-za(j1,j6)*zb(j4,j5)*(zb(j3,j1)*za(j1,j2)+zb(j3,j6)*za(j6,j2))
     & /(zb(j1,j6)*za(j4,j5)*zb(j5,j2)*za(j6,j3)*t(j1,j3,j6))
     &-za(j1,j6)*zb(j3,j5)*(zb(j4,j1)*za(j1,j2)+zb(j4,j6)*za(j6,j2))
     & /(zb(j1,j6)*za(j3,j5)*zb(j5,j2)*za(j6,j4)*t(j1,j4,j6))
      return
      end

      
      function amp_2gam2g_pmmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_pmmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +       -         -         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
      
c--- Equation (4.2)
c--- note: I have interpreted, e.g. <4|3+6|2> as [4|3+6|2> in order
c--- to get the correct <| and |] required by the helicities
      amp_2gam2g_pmmppm=      
     &+zb(j4,j5)*za(j2,j6)**2*(zb(j1,j2)*za(j2,j3)+zb(j1,j6)*za(j6,j3))
     & /(za(j4,j5)*zb(j5,j3)*za(j1,j6)*s(j1,j2)*t(j1,j2,j6))
     &+za(j3,j6)*zb(j5,j4)*(zb(j1,j3)*za(j3,j2)+zb(j1,j6)*za(j6,j2))
     & /(zb(j3,j6)*za(j5,j4)*zb(j5,j2)*za(j1,j6)*t(j1,j3,j6))
     &+za(j3,j6)*zb(j1,j5)**2*(zb(j4,j3)*za(j3,j2)+zb(j4,j6)*za(j6,j2))
     & /(zb(j3,j6)*zb(j2,j5)*za(j6,j4)*s(j1,j2)*t(j1,j2,j5))
     &+zb(j5,j1)*za(j6,j2)*(zb(j5,j1)*za(j1,j6)+zb(j5,j4)*za(j4,j6))
     & *t(j1,j4,j5)/(zb(j3,j6)*zb(j5,j2)*zb(j5,j3)
     &              *za(j6,j1)*za(j6,j4)*za(j4,j5)*s(j1,j2))
     &+(zb(j1,j3)*za(j3,j2)+zb(j1,j6)*za(j6,j2))
     & *(zb(j5,j1)*za(j1,j6)*(zb(j5,j2)*za(j2,j6)+zb(j5,j3)*za(j3,j6))
     &  +zb(j5,j4)*za(j4,j6)*zb(j5,j2)*za(j2,j6))
     & /(zb(j3,j6)*zb(j5,j2)*zb(j5,j3)*za(j6,j1)
     &  *za(j6,j4)*za(j4,j5)*s(j1,j2))

      return
      end
    
      
      function amp_2gam2g_mpmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2g_mpmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      -       +         -         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
      
c--- Equation (4.3)
c--- note: I have interpreted, e.g. <4|2+5|3> as [4|2+5|3> in order
c--- to get the correct <| and |] required by the helicities
      amp_2gam2g_mpmppm=      
     &+za(j1,j6)*zb(j2,j5)*(zb(j4,j2)*za(j2,j3)+zb(j4,j5)*za(j5,j3))
     & /(zb(j1,j6)*za(j2,j5)*zb(j5,j3)*za(j6,j4)*t(j1,j4,j6))
     &+za(j1,j6)*zb(j4,j5)*zb(j6,j2)
     & *(zb(j2,j4)*za(j4,j3)+zb(j2,j5)*za(j5,j3))
     & /(zb(j1,j6)*za(j4,j5)*zb(j5,j3)*s(j1,j2)*t(j1,j2,j6))
     &+za(j3,j6)*zb(j2,j5)*za(j5,j1)
     & *(zb(j4,j2)*za(j2,j1)+zb(j4,j5)*za(j5,j1))
     & /(zb(j3,j6)*za(j2,j5)*za(j6,j4)*s(j1,j2)*t(j3,j4,j6))
     &+za(j1,j6)*zb(j2,j5)*(zb(j2,j4)*za(j4,j3)+zb(j2,j5)*za(j5,j3))
     & /(zb(j1,j6)*zb(j3,j5)*za(j5,j4)*za(j4,j6)*s(j1,j2))
     &+t(j1,j3,j6)*(za(j5,j1)*zb(j2,j6)*zb(j5,j3)*za(j3,j6)
     &             +zb(j1,j2)*zb(j6,j5)*za(j1,j6)*za(j5,j1))
     & /(zb(j1,j6)*zb(j5,j3)*zb(j3,j6)*za(j2,j5)
     &  *za(j5,j4)*za(j4,j6)*s(j1,j2))

      return
      end
      
