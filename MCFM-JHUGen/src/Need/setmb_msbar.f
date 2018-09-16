      subroutine setmb_msbar
      implicit none
      include 'types.f'
c--- set up the value of the running b-mass in the MS-bar scheme,
c--- evaluated at the pole mass (mb)
c--- expressions are taken from 
c--- J.~Fleischer, F.~Jegerlehner, O.~V.~Tarasov and O.~L.~Veretin,
c--- %``Two-loop {QCD} corrections of the massive fermion propagator,''
c--- Nucl.\ Phys.\ B {\bf 539}, 671 (1999)
c--- [Erratum-ibid.\ B {\bf 571}, 511 (2000)]
c--- [arXiv:hep-ph/9803493].
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'couple.f'
      include 'kpart.f'
      real(dp):: alphas,c1,c2,a,zeta3,i31 
      parameter(zeta3=1.20205690315959428539_dp)

      if (mb_msbar < 0._dp) then
        c1=4._dp*Cf
        c2=0._dp
c--- calculate the MS-bar mass from the pole mass
        if (kpart==klord) then
          a=alphas(mb,amz,1)/(4._dp*pi)
        else
          a=alphas(mb,amz,2)/(4._dp*pi)
          i31=3._dp/2._dp*zeta3-6._dp*pisqo6*log(2._dp)
          c2=Cf*xn*(1111._dp/24._dp-8._dp*pisqo6-4._dp*i31)
     &      -Cf*half*real(nf-1,dp)*(71._dp/6._dp+8._dp*pisqo6)
     &      +Cf**2*(121._dp/8._dp+30._dp*pisqo6+8._dp*i31)
     &      -12._dp*Cf*half*(1._dp-2._dp*pisqo6)
        endif
        mb_msbar=mb/(1._dp+c1*a+c2*a**2)
      endif

c--- For comparison, these were the choices made in previous publications:
c--- mb(mb)=4.20 is our choice for the H+b paper
c--- mb(mb)=4.25 is the value used in the Les Houches write-up
 
      write(6,99) mb_msbar

      return

 99   format(/,
     &       ' ************* Running b-mass at pole mass **********'/, 
     &       ' *                                                  *'/, 
     &       ' *                mb_MSbar(mb)  = ',f8.4,'          *'/,
     &       ' ****************************************************')

      end
      
