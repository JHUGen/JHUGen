      function HAQggvsqanal(j1,j2,j3,j4)
      implicit none
      include 'types.f'
      real(dp):: HAQggvsqanal
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'scale.f' 
c      include 'masses.f' 
c      include 'deltar.f'
C---  matrix element squared for 0 --> H + a(j1)+q(j2)+g(j3)+g(j4)
c---  implemented according to arXiv:0906.0008, Eq. (2.23)
      integer:: j1,j2,j3,j4,h1,h2,h3
      complex(dp):: A0ab(2,2,2),A0ba(2,2,2),
     & A41ab(2,2,2),A41ba(2,2,2),A43ab(2,2,2),A43ba(2,2,2)
      real(dp):: temp,ren,H4prenorm
      
      call   Amplo_AQgg(j1,j2,j3,j4,A0ab,A0ba)
      call Ampvirt_AQgg(j1,j2,j3,j4,A41ab,A41ba,A43ab,A43ba)

c--- get renormalization factor
      ren=H4prenorm()

      temp=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      A41ab(h1,h2,h3)=A41ab(h1,h2,h3)+ren*A0ab(h1,h2,h3)
      A41ba(h1,h2,h3)=A41ba(h1,h2,h3)+ren*A0ba(h1,h2,h3)
c--- Note: A43 receives no renormalization

      temp=temp+real(conjg(A0ab(h1,h2,h3))*
     & (V*A41ab(h1,h2,h3)-A41ba(h1,h2,h3)+A43ab(h1,h2,h3)) 
     &              +conjg(A0ba(h1,h2,h3))*
     & (V*A41ba(h1,h2,h3)-A41ab(h1,h2,h3)+A43ba(h1,h2,h3))) 
      enddo
      enddo
      enddo

c--- Note: additional factor of 1/4 here due to difference between
c--- definition of overall Hgg coupling C^2 (sentence following (2.1))
c--- and our factor, "Asq" in gg_hgg_v.f
      HAQggvsqanal=V*temp/4._dp
            
      return
      end


