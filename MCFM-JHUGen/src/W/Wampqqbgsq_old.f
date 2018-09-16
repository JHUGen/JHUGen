      function Wampqqbgsq(p1,p2,p3,p4,p5,region,za,zb)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'hpls.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     index is helicity of outgoing gluon line
      real(dp):: Wampqqbgsq,u,v
      integer:: p1,p2,p3,p4,p5,region,NFV
      complex(dp):: al,be,ga,de,alC,beC,gaC,deC
      complex(dp):: amp(2),ampqqbgll
      real(dp):: alR,alI,
     & Alpha_2a0re,Alpha_2a0im,
     & Alpha_3a0re,Alpha_3a0im,
     & Alpha_4a0re,Alpha_4a0im
      real(dp):: beR,beI,
     & Beta_2a0re,Beta_2a0im,
     & Beta_3a0re,Beta_3a0im,
     & Beta_4a0re,Beta_4a0im
      real(dp):: gaR,gaI,
     & Gamma_2a0re,Gamma_2a0im,
     & Gamma_3a0re,Gamma_3a0im,
     & Gamma_4a0re,Gamma_4a0im

      NFV=0

      u=-s(p1,p3)/s(p1,p2)
      v=+s(p4,p5)/s(p1,p2)
      
      if (region == 2) then
        u=-s(p1,p3)/s(p1,p2)
        v=+s(p4,p5)/s(p1,p2)
      elseif (region == 3) then
        u=-s(p2,p3)/s(p1,p3)
        v=+s(p4,p5)/s(p1,p3)
      elseif (region == 4) then
        u=-s(p1,p3)/s(p2,p3)
        v=+s(p4,p5)/s(p2,p3)
      endif
      
      call tdhpl(u,v,4,G1,G2,G3,G4,H1,H2,H3,H4)

c      write(6,*) 'region',region
c      write(6,*) 's(p1,p2)',s(p1,p2)
c      write(6,*) 's(p1,p3)',s(p1,p3)
c      write(6,*) 's(p2,p3)',s(p2,p3)
c      write(6,*) 'u,v',u,v
c      pause
      
      if (region == 2) then
      alR=Alpha_2a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      alI=Alpha_2a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      beR=Beta_2a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      beI=Beta_2a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      gaR=Gamma_2a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      gaI=Gamma_2a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)

      elseif (region == 3) then
      alR=Alpha_3a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      alI=Alpha_3a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      beR=Beta_3a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      beI=Beta_3a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      gaR=Gamma_3a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      gaI=Gamma_3a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)

      elseif (region == 4) then

      alR=Alpha_4a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      alI=Alpha_4a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      beR=Beta_4a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      beI=Beta_4a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      gaR=Gamma_4a0re(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
      gaI=Gamma_4a0im(u,v,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)

      endif

      al=cmplx(alR,alI,dp)
      be=cmplx(beR,beI,dp)
      ga=cmplx(gaR,gaI,dp)
      de=s(p1,p2)*(al-be-ga)/(2._dp*(s(p1,p2)+s(p2,p3)+s(p3,p1)))

      alC=conjg(al)
      beC=conjg(be)
      gaC=conjg(ga)
      deC=conjg(de)

!RH gluon
      amp(2)=-ampqqbgll(p2,p1,p3,p4,p5,al, be, ga, de ,za,zb)
!LH gluon
      amp(1)=+ampqqbgll(p1,p2,p3,p5,p4,alC,beC,gaC,deC,zb,za)

      Wampqqbgsq=
     & +real(amp(1)*conjg(amp(1)),dp)+real(amp(2)*conjg(amp(2)),dp)

      return
      end


