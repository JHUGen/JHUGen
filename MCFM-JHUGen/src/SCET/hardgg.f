      subroutine hardgg(Mhsq,musq,coeff)
      implicit none
!   Wilson Coefficients for gg->H in units of as/4/pi
!   Taken from Ahrens,Becher,Neubert,Yang, 0809.4283v3
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp),intent(in)::Mhsq,musq
      real(dp),intent(out)::coeff(2)
      real(dp)::Ct(2),logmtsq
      complex(dp)::L,lnrat
      complex(dp)::Cs(2)

      L=lnrat(-Mhsq,musq)

!----Eq.(17) of 0809.4283v3
      Cs(1)=CA*(-L**2+zeta2)
      Cs(2)=
     & +CA**2*(1/2._dp*L**4+11/9._dp*L**3+(-67/9._dp+zeta2)*L**2
     & +L*(80/27._dp-22/3._dp*zeta2-2*zeta3)
     & +5105/162._dp+67/6._dp*zeta2+1/2._dp*zeta2**2-143/9._dp*zeta3)
     & +TR*nf*CF*(4*L-67/3._dp+16*zeta3)
     & +TR*nf*CA*(
     & -4/9._dp*L**3+20/9._dp*L**2+(104/27._dp+8/3._dp*zeta2)*L
     & -1832/81._dp-10/3._dp*zeta2-92/9._dp*zeta3)

!----Eq.(12) of 0809.4283v3
      logmtsq=log(mt**2/musq)
      Ct(1)=5*CA-3*CF
      Ct(2)=27/2._dp*CF**2+(11*logmtsq-100/3._dp)*CF*CA
     & -(7*logmtsq-1063/36._dp)*CA**2
     & -4/3._dp*CF*TR-5/6._dp*CA*TR
     & -(8*logmtsq+5._dp)*CF*TR*nf-47/9._dp*CA*TR*nf

! coeff is returned in units of as/2/pi
      coeff(1)=Ct(1)+real(Cs(1),dp)
      coeff(2)=(Ct(1)**2+2*Ct(2)+2*real(Cs(2),dp)
     & +real(Cs(1)*conjg(Cs(1)),dp)
     & +four*real(Cs(1),dp)*Ct(1))/four

      return
      end

