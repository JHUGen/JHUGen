      subroutine setmb_msbar
c--- set up the value of the running b-mass in the MS-bar scheme,
c--- evaluated at the pole mass (mb)
c--- expressions are taken from 
c--- J.~Fleischer, F.~Jegerlehner, O.~V.~Tarasov and O.~L.~Veretin,
c--- %``Two-loop {QCD} corrections of the massive fermion propagator,''
c--- Nucl.\ Phys.\ B {\bf 539}, 671 (1999)
c--- [Erratum-ibid.\ B {\bf 571}, 511 (2000)]
c--- [arXiv:hep-ph/9803493].
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'couple.f'
      include 'part.f'
      double precision alphas,c1,c2,a,zeta3,i31 
      parameter(zeta3=1.20205690315959428539d0)

      if (mb_msbar .lt. 0d0) then
        c1=4d0*Cf
        c2=0d0
c--- calculate the MS-bar mass from the pole mass
        if (part .eq. 'lord') then
          a=alphas(mb,amz,1)/(4d0*pi)
        else
          a=alphas(mb,amz,2)/(4d0*pi)
          i31=3d0/2d0*zeta3-6d0*pisqo6*dlog(2d0)
          c2=Cf*xn*(1111d0/24d0-8d0*pisqo6-4d0*i31)
     .      -Cf*half*dfloat(nf-1)*(71d0/6d0+8d0*pisqo6)
     .      +Cf**2*(121d0/8d0+30d0*pisqo6+8d0*i31)
     .      -12d0*Cf*half*(1d0-2d0*pisqo6)
        endif
        mb_msbar=mb/(1d0+c1*a+c2*a**2)
      endif

c--- For comparison, these were the choices made in previous publications:
c--- mb(mb)=4.20 is our choice for the H+b paper
c--- mb(mb)=4.25 is the value used in the Les Houches write-up
 
      write(6,99) mb_msbar

      return

 99   format(/,
     .       ' ************* Running b-mass at pole mass **********'/, 
     .       ' *                                                  *'/, 
     .       ' *                mb_MSbar(mb)  = ',f8.4,'          *'/,
     .       ' ****************************************************')

      end
      
