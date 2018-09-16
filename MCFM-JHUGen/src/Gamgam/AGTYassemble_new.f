
      subroutine AGTYassemble_new(ss,tt,uu,F1x1,F2x0,Remain)
!======C.Williams August 2015
!======routine to reproduce eq. 4.6 of AGTY, without overall
!=====factor of N_c, which lives with the LO in our units (factor of 2 is included) 

!=====Apie is the piece which goes like sum over quark charges
!=====Qpie is the piece which goes like the LO quark charge

!=====note that neither Apie nor Qpie is dressed with EM charges in this
!===== routine (for consistency) 
      
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f' 
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'scet_const.f'
      real(dp), intent(in) :: ss,tt,uu
      real(dp):: ddilog,Li3,Li4

      real(dp), intent(out) :: F2x0,F1x1,Remain
!=======basis functions from AGTY
      real(dp) :: AGTYAs,AGTYBs,AGTYCs,AGTYD1s,AGTYD2s,AGTYE1s
      real(dp) :: AGTYE2s,AGTYE3s,AGTYE4s,AGTYF1s,AGTYF2s,AGTYG1s

      real(dp) :: Bigx,Bigy,Bigs,BigT,BigU,x,y,z
      real(dp) Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li2z,Li3z,Li4z,Li4zinv
      
      
!----- define various pieces 
      BigX=log(-tt/ss)
      BigY=log(-uu/ss) 
      BigS=log(ss/musq) 
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)
      
      x=-tt/ss
      y=-uu/ss 
      z=-uu/tt

      Li2x=ddilog(x)
      Li3x=Li3(x)
      Li4x=Li4(x)

      Li2y=ddilog(y)
      Li3y=Li3(y)
      Li4y=Li4(y)

      Li2z=ddilog(z)
      Li3z=Li3(z)
      Li4z=Li4(z)

      Li4zinv=Li4(1._dp/z)


!======== now build the functions

!====== this bit is special and doesnt scale with the LO quark charge. 
      Remain=two*cf*tr*AGTYAs(ss,tt,uu,BigX,BigY,Li2x,Li2y,Li3x,Li3y,
     & Li4x,Li4y,Li4z,Li4zinv)

!======= B piecs (CF**2)
      F2x0=two*Cf**2*AGTYBs(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)

!======= D2 piecs (CF*CA)
      F2x0=F2x0
     & +two*CF*CA*AGTYD2s(ss,tt,uu,BigX,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)

!======= E3 pieces (NF*CF)
      F2x0=F2x0+two*nf*Cf*AGTYE3s(tt,uu,BigX,Bigy,Bigs)

      F1x1=CF**2*AGTYG1s(tt,uu,BigX,BigY)

      
      return
      end 
