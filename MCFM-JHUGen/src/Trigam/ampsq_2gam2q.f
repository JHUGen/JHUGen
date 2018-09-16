      subroutine ampsq_2gam2q(j1,j2,j3,j4,j5,j6,za,zb,
     & ampsq,ampsqid)
      implicit none
      include 'types.f'
c--- Returns matrix element squared for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c---
c--- returns arrays of matrix elements for non-identical quarks (ampsq)
c--- and for identical quarks (ampsqid), indexed by charges of
c--- (j1,j2) and (j3,j4) quark lines
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
      integer:: j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,h5,j12,j34
      complex(dp):: amp(2,2,2,2,2,2),amp_swap(2,2,2,2,2,2),
     & amp_2gam2q_mpmppp,amp_2gam2q_mppmpp,
     & amp_2gam2q_pmmppp,amp_2gam2q_pmpmpp,
     & amp_2gam2q_pmpmmp,amp_2gam2q_pmmpmp
      real(dp):: Q12,Q34,ampsq(2,2),ampsqid(2,2)
              
c--- ordering of labels in amp is as follows:
c---  (charge of quark line 12, charge of quark line 34,
c---   helicity of antiquark j1, antiquark j3, photon j5, photon j6;
c---   quark helicity j2=3-j1, quark j4=3-j3)

c--- loop over possible charges of quark lines
      do j12=1,2
      do j34=1,2

      if (j12 == 1) then
        Q12=-1._dp/3._dp
      else
        Q12=+2._dp/3._dp
      endif

      if (j34 == 1) then
        Q34=-1._dp/3._dp
      else
        Q34=+2._dp/3._dp
      endif

c--- basic amplitudes
      amp(j12,j34,1,1,2,2)=
     & amp_2gam2q_mpmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      amp(j12,j34,1,2,2,2)=
     & amp_2gam2q_mppmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      amp(j12,j34,2,1,2,2)=
     & amp_2gam2q_pmmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      amp(j12,j34,2,2,2,2)=
     & amp_2gam2q_pmpmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)

      amp(j12,j34,2,2,1,2)=
     & amp_2gam2q_pmpmmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)
      amp(j12,j34,2,1,1,2)=
     & amp_2gam2q_pmmpmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34)

c--- simple symmetry
      amp(j12,j34,2,2,2,1)=
     & amp_2gam2q_pmpmmp(j1,j2,j3,j4,j6,j5,za,zb,Q12,Q34)
      amp(j12,j34,2,1,2,1)=
     & amp_2gam2q_pmmpmp(j1,j2,j3,j4,j6,j5,za,zb,Q12,Q34)
      
c--- remaining amplitudes by complex conjugation
      do h1=1,2
      do h2=1,2
      amp(j12,j34,h1,h2,1,1)=conjg(amp(j12,j34,3-h1,3-h2,2,2))
      enddo
      enddo
      amp(j12,j34,1,1,2,1)=conjg(amp(j12,j34,2,2,1,2))    
      amp(j12,j34,1,2,2,1)=conjg(amp(j12,j34,2,1,1,2))    
      amp(j12,j34,1,1,1,2)=conjg(amp(j12,j34,2,2,2,1))    
      amp(j12,j34,1,2,1,2)=conjg(amp(j12,j34,2,1,2,1))    

c--- additional amplitudes needed for identical quarks 
c--- (obtained by interchanging j2 and j4)     
      if (j12. eq. j34) then
c--- basic amplitudes
      amp_swap(j12,j34,1,1,2,2)=
     & amp_2gam2q_mpmppp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34)
      amp_swap(j12,j34,1,2,2,2)=
     & amp_2gam2q_mppmpp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34)
      amp_swap(j12,j34,2,1,2,2)=
     & amp_2gam2q_pmmppp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34)
      amp_swap(j12,j34,2,2,2,2)=
     & amp_2gam2q_pmpmpp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34)

      amp_swap(j12,j34,2,2,1,2)=
     & amp_2gam2q_pmpmmp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34)
      amp_swap(j12,j34,2,1,1,2)=
     & amp_2gam2q_pmmpmp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34)

c--- simple symmetry
      amp_swap(j12,j34,2,2,2,1)=
     & amp_2gam2q_pmpmmp(j1,j4,j3,j2,j6,j5,za,zb,Q12,Q34)
      amp_swap(j12,j34,2,1,2,1)=
     & amp_2gam2q_pmmpmp(j1,j4,j3,j2,j6,j5,za,zb,Q12,Q34)
      
c--- remaining amplitudes by complex conjugation
      do h1=1,2
      do h2=1,2
      amp_swap(j12,j34,h1,h2,1,1)
     & =conjg(amp_swap(j12,j34,3-h1,3-h2,2,2))
      enddo
      enddo
      amp_swap(j12,j34,1,1,2,1)=conjg(amp_swap(j12,j34,2,2,1,2))    
      amp_swap(j12,j34,1,2,2,1)=conjg(amp_swap(j12,j34,2,1,1,2))    
      amp_swap(j12,j34,1,1,1,2)=conjg(amp_swap(j12,j34,2,2,2,1))    
      amp_swap(j12,j34,1,2,1,2)=conjg(amp_swap(j12,j34,2,1,2,1))    

      endif
      
      enddo
      enddo
 
c--- square amplitudes     
      do j12=1,2
      do j34=1,2
      
      ampsq(j12,j34)=0._dp
      ampsqid(j12,j34)=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      ampsq(j12,j34)=ampsq(j12,j34)+abs(amp(j12,j34,h1,h2,h3,h4))**2
      if (j12 == j34) then
      ampsqid(j12,j34)=ampsqid(j12,j34)
     & +abs(amp(j12,j34,h1,h2,h3,h4))**2
     & +abs(amp_swap(j12,j34,h1,h2,h3,h4))**2
        if (h1 == h2) then
        ampsqid(j12,j34)=ampsqid(j12,j34)
     &   +2._dp/xn*real(amp(j12,j34,h1,h2,h3,h4)
     &                *conjg(amp_swap(j12,j34,h1,h2,h3,h4)))
        endif
      endif
      enddo
      enddo
      enddo
      enddo
      
      enddo
      enddo
      
      return
      end
      
