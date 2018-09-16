      subroutine qqcoeff(Qsq,musq,coeff)
      implicit none
!   Wilson Coefficients for qq->V in units of 
!   as/4/pi
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp),intent(in)::Qsq,musq
      complex(dp)::L,lnrat
      complex(dp),intent(out)::coeff(2)
      
      L=lnrat(-Qsq,musq)
      
      coeff(1)=CF*(-(L**2-3*L+8._dp-zeta2))
      coeff(2)=CF**2*(1/2._dp*(L**2-3*L+8._dp-zeta2)**2
     & +(3/2._dp-12*zeta2+24*zeta3)*L
     & -1/8._dp+29._dp*zeta2-30*zeta3-44/5._dp*zeta2**2)
     & +CF*nf*(-2/9._dp*L**3+19/9._dp*L**2-(209/27._dp+4/3._dp*zeta2)*L
     & +4085/324._dp+23/9._dp*zeta2+2/9._dp*zeta3)
     & +CF*CA*(11/9._dp*L**3+(2*zeta2-233/18._dp)*L**2
     & +(2545/54._dp+22/3._dp*zeta2-26*zeta3)*L
     & -51157/648._dp-zeta2*(337/18._dp-44/5._dp*zeta2)+313/9._dp*zeta3)
     
      return
      end

