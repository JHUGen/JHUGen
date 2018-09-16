      subroutine zgampsbis(order,p1,p2,p3,p4,p5,region,za,zb,amps2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'hpls.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     1st index is helicity quark line
!     2nd index is helicity lepton line
!     helicity of gluon is summe.e-_dpover
!
!     Order of calculation is passed in via integer "order"
!     Region is passed in via integer "region"
!
      integer:: p1,p2,p3,p4,p5,hq,hl,k,NFV,region,order
      complex(dp):: al,be,ga,de,alC,beC,gaC,deC
      complex(dp):: amps(0:2,2,2,2),ampqqbgll
      real(dp):: uu,vv,amps2(2,2)
      real(dp):: alR(0:2),alI(0:2),
     & Alpha_2a1re,Alpha_2a1im,
     & Alpha_3a1re,Alpha_3a1im,
     & Alpha_4a1re,Alpha_4a1im
      real(dp):: beR(0:2),beI(0:2),
     & Beta_2a1re,Beta_2a1im,
     & Beta_3a1re,Beta_3a1im,
     & Beta_4a1re,Beta_4a1im
      real(dp):: gaR(0:2),gaI(0:2),
     & Gamma_2a1re,Gamma_2a1im,
     & Gamma_3a1re,Gamma_3a1im,
     & Gamma_4a1re,Gamma_4a1im
      real(dp):: 
     & Alpha_2a0re,Alpha_2a0im,
     & Alpha_3a0re,Alpha_3a0im,
     & Alpha_4a0re,Alpha_4a0im
      real(dp):: 
     & Beta_2a0re,Beta_2a0im,
     & Beta_3a0re,Beta_3a0im,
     & Beta_4a0re,Beta_4a0im
      real(dp):: 
     & Gamma_2a0re,Gamma_2a0im,
     & Gamma_3a0re,Gamma_3a0im,
     & Gamma_4a0re,Gamma_4a0im
      real(dp):: 
     & Alpha_2a2re,Alpha_2a2im,
     & Alpha_3a2re,Alpha_3a2im,
     & Alpha_4a2re,Alpha_4a2im
      real(dp):: 
     & Beta_2a2re,Beta_2a2im,
     & Beta_3a2re,Beta_3a2im,
     & Beta_4a2re,Beta_4a2im
      real(dp):: 
     & Gamma_2a2re,Gamma_2a2im,
     & Gamma_3a2re,Gamma_3a2im,
     & Gamma_4a2re,Gamma_4a2im

      NFV=0


      if (region == 2) then
        uu=-s(p1,p3)/s(p1,p2)
        vv=+s(p4,p5)/s(p1,p2)
      elseif (region == 3) then
        uu=-s(p2,p3)/s(p1,p3)
        vv=+s(p4,p5)/s(p1,p3)
      elseif (region == 4) then
        uu=-s(p1,p3)/s(p2,p3)
        vv=+s(p4,p5)/s(p2,p3)
      endif
      
      if (order > 0) then
c--- fill arrays for 2DHPLs
        call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)
      else
c--- fill with dummy values
        G1=zip
        G2=zip
        G3=zip
        G4=zip
        H1=zip
        H2=zip
        H3=zip
        H4=zip
      endif

c--- Coefficients for Region 3
      if (region == 2) then
        alR(0)=Alpha_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        alI(0)=Alpha_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        beR(0)=Beta_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        beI(0)=Beta_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        gaR(0)=Gamma_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        gaI(0)=Gamma_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)

        if (order > 0) then
          alR(1)=Alpha_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          alI(1)=Alpha_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beR(1)=Beta_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beI(1)=Beta_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaR(1)=Gamma_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaI(1)=Gamma_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        endif
      
        if (order > 1) then
          alR(2)=Alpha_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          alI(2)=Alpha_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beR(2)=Beta_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beI(2)=Beta_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaR(2)=Gamma_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaI(2)=Gamma_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        endif
      
c--- Coefficients for Region 2
      elseif (region == 3) then
        alR(0)=Alpha_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        alI(0)=Alpha_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        beR(0)=Beta_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        beI(0)=Beta_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        gaR(0)=Gamma_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        gaI(0)=Gamma_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)

        if (order > 0) then
          alR(1)=Alpha_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          alI(1)=Alpha_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beR(1)=Beta_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beI(1)=Beta_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaR(1)=Gamma_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaI(1)=Gamma_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        endif

        if (order > 1) then
          alR(2)=Alpha_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          alI(2)=Alpha_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beR(2)=Beta_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beI(2)=Beta_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaR(2)=Gamma_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaI(2)=Gamma_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        endif

c--- Coefficients for Region 4
      elseif (region == 4) then
        alR(0)=Alpha_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        alI(0)=Alpha_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        beR(0)=Beta_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        beI(0)=Beta_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        gaR(0)=Gamma_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        gaI(0)=Gamma_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)

        if (order > 0) then
          alR(1)=Alpha_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          alI(1)=Alpha_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beR(1)=Beta_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          beI(1)=Beta_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaR(1)=Gamma_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
          gaI(1)=Gamma_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        endif

        if (order > 1) then
          alR(2)=Alpha_4a2re(uu,vv,H2,H2,H3,H4,G2,G2,G3,G4,nf,NFV)
          alI(2)=Alpha_4a2im(uu,vv,H2,H2,H3,H4,G2,G2,G3,G4,nf,NFV)
          beR(2)=Beta_4a2re(uu,vv,H2,H2,H3,H4,G2,G2,G3,G4,nf,NFV)
          beI(2)=Beta_4a2im(uu,vv,H2,H2,H3,H4,G2,G2,G3,G4,nf,NFV)
          gaR(2)=Gamma_4a2re(uu,vv,H2,H2,H3,H4,G2,G2,G3,G4,nf,NFV)
          gaI(2)=Gamma_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,nf,NFV)
        endif

      endif

      do k=0,2
      al=cmplx(alR(k),alI(k),dp)
      be=cmplx(beR(k),beI(k),dp)
      ga=cmplx(gaR(k),gaI(k),dp)
      de=s(p1,p2)*(al-be-ga)/(2._dp*(s(p1,p2)+s(p2,p3)+s(p3,p1)))

      alC=conjg(al)
      beC=conjg(be)
      gaC=conjg(ga)
      deC=conjg(de)

!     1st index is order of amplitude
!     2nd index is helicity quark line
!     3rd index is helicity lepton line
!     4th index is helicity of outgoing gluon line
      amps(k,1,1,2)=-ampqqbgll(p2,p1,p3,p4,p5,al, be, ga, de ,za,zb)
      amps(k,2,2,1)=-ampqqbgll(p2,p1,p3,p4,p5,alC,beC,gaC,deC,zb,za)
      amps(k,1,2,2)=-ampqqbgll(p2,p1,p3,p5,p4,al, be, ga, de ,za,zb)
      amps(k,2,1,1)=-ampqqbgll(p2,p1,p3,p5,p4,alC,beC,gaC,deC,zb,za)

      amps(k,1,1,1)=+ampqqbgll(p1,p2,p3,p5,p4,alC,beC,gaC,deC,zb,za)
      amps(k,2,2,2)=+ampqqbgll(p1,p2,p3,p5,p4,al, be, ga, de ,za,zb)
      amps(k,1,2,1)=+ampqqbgll(p1,p2,p3,p4,p5,alC,beC,gaC,deC,zb,za)
      amps(k,2,1,2)=+ampqqbgll(p1,p2,p3,p4,p5,al, be, ga, de ,za,zb)
      enddo
 
      do hq=1,2
      do hl=1,2
        if (order == 0) then
          amps2(hq,hl)=abs(amps(0,hq,hl,1))**2+abs(amps(0,hq,hl,2))**2
        elseif (order == 1) then
          amps2(hq,hl)=
     &     +real(amps(0,hq,hl,1)*conjg(amps(1,hq,hl,1)),dp)
     &     +real(amps(0,hq,hl,2)*conjg(amps(1,hq,hl,2)),dp)
        else
          write(6,*) 'zgampsbis not implemented for order = ',order
          stop
        endif
      enddo
      enddo
      
      return
      end


