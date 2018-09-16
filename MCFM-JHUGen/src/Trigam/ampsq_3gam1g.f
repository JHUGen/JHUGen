      function ampsq_3gam1g(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      real(dp):: ampsq_3gam1g
c--- Matrix element squared for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + gam(j2) + gam(j3) + gam(j4)
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer:: j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,h5
      complex(dp):: amp(2,2,2,2,2),amp_3gam1g_mppppm,
     & amp_3gam1g_pmpppm,
     & amp_3gam1g_ppmppm,amp_3gam1g_pppmpm,amp_3gam1g_mmpppm
      
c--- ordering of labels in amp is as follows:
c---  (gluon j1, photon j2, photon j3, photon j4, antiquark j5)
c--- for consistency with function names (and quark j6=3-j5)

c--- trivial amplitude
      amp(2,2,2,2,2)=czip
c--- basic MHV amplitudes
      amp(1,2,2,2,2)=amp_3gam1g_mppppm(j1,j2,j3,j4,j5,j6,za,zb)      
      amp(2,1,2,2,2)=amp_3gam1g_pmpppm(j1,j2,j3,j4,j5,j6,za,zb)      
      amp(2,2,1,2,2)=amp_3gam1g_ppmppm(j1,j2,j3,j4,j5,j6,za,zb)      
      amp(2,2,2,1,2)=amp_3gam1g_pppmpm(j1,j2,j3,j4,j5,j6,za,zb)      

c--- non-MHV amplitude and permutations
      amp(1,1,2,2,2)=amp_3gam1g_mmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(1,2,1,2,2)=amp_3gam1g_mmpppm(j1,j3,j2,j4,j5,j6,za,zb)
      amp(1,2,2,1,2)=amp_3gam1g_mmpppm(j1,j4,j2,j3,j5,j6,za,zb)            

c--- parity and charge conjugation
      amp(1,1,1,1,2)=czip
      
      amp(2,1,1,1,2)=-amp_3gam1g_pppmpm(j4,j3,j2,j1,j6,j5,za,zb)      
      amp(1,2,1,1,2)=-amp_3gam1g_ppmppm(j4,j3,j2,j1,j6,j5,za,zb)      
      amp(1,1,2,1,2)=-amp_3gam1g_pmpppm(j4,j3,j2,j1,j6,j5,za,zb)      
      amp(1,1,1,2,2)=-amp_3gam1g_mppppm(j4,j3,j2,j1,j6,j5,za,zb)      
      
      amp(2,2,1,1,2)=-amp_3gam1g_mmpppm(j4,j3,j2,j1,j5,j6,za,zb)      
      amp(2,1,2,1,2)=-amp_3gam1g_mmpppm(j2,j4,j1,j3,j5,j6,za,zb)
      amp(2,1,1,2,2)=-amp_3gam1g_mmpppm(j2,j3,j1,j4,j5,j6,za,zb)  
            
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      amp(h1,h2,h3,h4,1)=amp(3-h1,3-h2,3-h3,3-h4,2)
      enddo
      enddo
      enddo
      enddo

c--- note: obvious redundancy in this routine, but might be
c--- worth checking relations for use in virtual
      
      ampsq_3gam1g=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
c      write(6,*) h1*10000+h2*1000+h3*100+h4*10+h5,
c     & abs(amp(h1,h2,h3,h4,h5))**2
c     & *8._dp*esq**3*gsq*xn*Cf/6._dp*aveqq*Q(2)**6*8._dp
      ampsq_3gam1g=ampsq_3gam1g+abs(amp(h1,h2,h3,h4,h5))**2
      enddo
      enddo
      enddo
      enddo
      enddo
c      pause

      return
      end
      
