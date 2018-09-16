      subroutine gagaAGTYassemble1loop(ss,tt,uu,M1finM0,BoldC0,M0sq,M1rensq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f' 
      real(dp), intent(in) :: ss,tt,uu
      complex(dp), intent(in) :: M1finM0,BoldC0,M0sq
      real(dp), intent(out) :: M1rensq
      real(dp) :: AGTYG1s,BigX,BigY

c--- AGTY function for |M1fin|^2 pieces
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      M1rensq=AGTYG1s(tt,uu,BigX,BigY)
      
      M1rensq=M1rensq+two*real(BoldC0*conjg(M1finM0),dp)
     & +real(BoldC0*conjg(BoldC0)*M0sq,dp)

c--- restore overall color factor
      M1rensq=M1rensq*CF**2
      
      return
      end
      
      
      
