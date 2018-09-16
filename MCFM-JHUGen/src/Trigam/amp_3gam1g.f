      function amp_3gam1g_mppppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_3gam1g_mppppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + gam(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      -         +         +         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=1, m=3 and labels permuted appropriately
c--- negative helicity gluon i=1
      amp_3gam1g_mppppm=      
     & za(j5,j1)*za(j6,j1)**3/(za(j5,j6)*za(j6,j1)*za(j1,j5))
     & *za(j6,j5)/(za(j6,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
           
      return
      end

      
      function amp_3gam1g_pmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_3gam1g_pmpppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + gam(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +         -         +         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=1, m=3 and labels permuted appropriately
c--- negative helicity photon i=2
      amp_3gam1g_pmpppm=      
     & za(j5,j2)*za(j6,j2)**3/(za(j5,j6)*za(j6,j1)*za(j1,j5))
     & *za(j6,j5)/(za(j6,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
           
      return
      end

      
      function amp_3gam1g_ppmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_3gam1g_ppmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + gam(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +         +         -         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=1, m=3 and labels permuted appropriately
c--- negative helicity photon i=3
      amp_3gam1g_ppmppm=      
     & za(j5,j3)*za(j6,j3)**3/(za(j5,j6)*za(j6,j1)*za(j1,j5))
     & *za(j6,j5)/(za(j6,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
           
      return
      end

      
      function amp_3gam1g_pppmpm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_3gam1g_pppmpm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + gam(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +         +         +         -
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=1, m=3 and labels permuted appropriately
c--- negative helicity photon i=4
      amp_3gam1g_pppmpm=      
     & za(j5,j4)*za(j6,j4)**3/(za(j5,j6)*za(j6,j1)*za(j1,j5))
     & *za(j6,j5)/(za(j6,j2)*za(j2,j5))
     & *za(j6,j5)/(za(j6,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
           
      return
      end

      
      function amp_3gam1g_mmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_3gam1g_mmpppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + gam(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      -         -         +         +
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
      integer:: j1,j2,j3,j4,j5,j6,k1(2),k2(2),k3(2),k4(2),i12,i34
      real(dp):: t
      
c--- Equation (5.1)
      amp_3gam1g_mmpppm=      
     & s(j5,j6)*t(j6,j1,j2)
     & *(za(j6,j3)*zb(j3,j5)+za(j6,j4)*zb(j4,j5))
     & /(za(j5,j3)*za(j3,j6)*za(j5,j4)*za(j4,j6)
     &  *zb(j5,j1)*zb(j1,j6)*zb(j5,j2)*zb(j2,j6))

c--- Initialize arrays for permutations
      k1(1)=j1
      k2(1)=j2
      k1(2)=j2
      k2(2)=j1
      k3(1)=j3
      k4(1)=j4
      k3(2)=j4
      k4(2)=j3
      
      do i12=1,2
      do i34=1,2
c--- There appear to be errors in this equation in the paper
c--- (even in the journal version) that I have corrected here
      amp_3gam1g_mmpppm=amp_3gam1g_mmpppm
     & +zb(j5,k3(i34))*za(j6,k2(i12))
     & *(zb(k4(i34),j5)*za(j5,k1(i12))
     &  +zb(k4(i34),k3(i34))*za(k3(i34),k1(i12)))
     & /(za(j5,k3(i34))*zb(j6,k2(i12))*za(j6,k4(i34))*zb(j5,k1(i12)))
     & /t(j6,k4(i34),k2(i12))      
      enddo
      enddo
      
      return
      end
      
