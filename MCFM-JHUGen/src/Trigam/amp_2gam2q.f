      function amp_2gam2q_mpmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2q_mpmppp
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  -       +        -         +        +         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34

c--- Equation (2.15) with m=2, r=0 (n=6) and f(+,+) from Eq. (2.13)
      amp_2gam2q_mpmppp=      
     & (Q12*za(j2,j1)/(za(j2,j5)*za(j5,j1))
     & +Q34*za(j4,j3)/(za(j4,j5)*za(j5,j3)))
     &*(Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
     &/(za(j1,j2)*za(j3,j4))
     &*(-za(j1,j3)**2)
          
      return
      end

      
      function amp_2gam2q_mppmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2q_mppmpp
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  -       +        +         -        +         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34

c--- Equation (2.15) with m=2, r=0 (n=6) and f(+,-) from Eq. (2.13)
      amp_2gam2q_mppmpp=      
     & (Q12*za(j2,j1)/(za(j2,j5)*za(j5,j1))
     & +Q34*za(j4,j3)/(za(j4,j5)*za(j5,j3)))
     &*(Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
     &/(za(j1,j2)*za(j3,j4))
     &*(+za(j1,j4)**2)
          
      return
      end

      
      function amp_2gam2q_pmmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2q_pmmppp
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  +       -        -         +        +         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34

c--- Equation (2.15) with m=2, r=0 (n=6) and f(-,+) from Eq. (2.13)
      amp_2gam2q_pmmppp=      
     & (Q12*za(j2,j1)/(za(j2,j5)*za(j5,j1))
     & +Q34*za(j4,j3)/(za(j4,j5)*za(j5,j3)))
     &*(Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
     &/(za(j1,j2)*za(j3,j4))
     &*(+za(j2,j3)**2)
          
      return
      end

      
      function amp_2gam2q_pmpmpp(j1,j2,j3,j4,j5,j6,za,zb,
     & Q12,Q34)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2q_pmpmpp
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  +       -        +         -        +         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34

c--- Equation (2.15) with m=2, r=0 (n=6) and f(-,-) from Eq. (2.13)
      amp_2gam2q_pmpmpp=      
     & (Q12*za(j2,j1)/(za(j2,j5)*za(j5,j1))
     & +Q34*za(j4,j3)/(za(j4,j5)*za(j5,j3)))
     &*(Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
     &/(za(j1,j2)*za(j3,j4))
     &*(-za(j2,j4)**2)
          
      return
      end

      
      function amp_2gam2q_pmpmmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2q_pmpmmp
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  +       -        +         -        -         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      complex(dp):: g1_2gam2q,g2_2gam2q

c--- Equation (4.5), first line
      amp_2gam2q_pmpmmp=      
     & +Q12**2*g1_2gam2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34**2*g1_2gam2q(j3,j4,j1,j2,j5,j6,za,zb)
     & +Q12*Q34*g2_2gam2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q12*Q34*g2_2gam2q(j3,j4,j1,j2,j5,j6,za,zb)
         
      return
      end

      
      function amp_2gam2q_pmmpmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam2q_pmmpmp
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  +       -        -         +        -         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      complex(dp):: g1_2gam2q,g2_2gam2q

c--- Equation (4.5), second line
      amp_2gam2q_pmmpmp=      
     & +Q12**2*g1_2gam2q(j1,j2,j4,j3,j5,j6,za,zb)
     & +Q34**2*g1_2gam2q(j4,j3,j1,j2,j5,j6,za,zb)
     & -Q12*Q34*g2_2gam2q(j1,j2,j4,j3,j5,j6,za,zb)
     & -Q12*Q34*g2_2gam2q(j4,j3,j1,j2,j5,j6,za,zb)
         
      return
      end

      
      function g1_2gam2q(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: g1_2gam2q
c--- -i*g1(j1,j2,j3,j4,j5,j6): auxiliary function for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
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

c--- Equation (3.9)
      g1_2gam2q=
     &  (zb(j1,j3)*za(j3,j4)+zb(j1,j6)*za(j6,j4))
     &  *(zb(j3,j1)*za(j1,j2)+zb(j3,j6)*za(j6,j2))
     &  /(za(j1,j6)*za(j6,j2)*zb(j1,j5)*zb(j5,j2)*s(j3,j4))
     & +zb(j1,j6)*za(j2,j4)*(zb(j3,j1)*za(j1,j5)+zb(j3,j6)*za(j6,j5))
     &  /(za(j1,j6)*zb(j1,j5)*s(j3,j4)*t(j1,j5,j6))
     & +za(j5,j2)*zb(j3,j1)*(zb(j6,j2)*za(j2,j4)+zb(j6,j5)*za(j5,j4))
     &  /(zb(j5,j2)*za(j6,j2)*s(j3,j4)*t(j1,j3,j4))
     
      return
      end
            
            
      function g2_2gam2q(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: g2_2gam2q
c--- -i*g2(j1,j2,j3,j4,j5,j6): auxiliary function for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
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

c--- Equation (3.10)
      g2_2gam2q=
     & -(zb(j3,j1)*za(j1,j2)+zb(j3,j6)*za(j6,j2))**2
     &  /(za(j1,j6)*za(j6,j2)*zb(j3,j5)*zb(j5,j4)*t(j1,j2,j6))
     
      return
      end
      
      
     
      
