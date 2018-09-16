!==== This file contains the box coefficients for the Q(i2) pieces of the following amplitude 

!===== (i1)^- (i2)^+ (i3)^- + (i4)^+ + gamma(i5) + gluon(i6) 

! where i3 and i4 are the leptons from a W decay. 

!=== C. Williams Nov 2013 

!================================================================================================      

!======= Q2 h5=1 h6=1 LEADING COLOR AMPLITUDES (order matches KC)
!=====  s(34)

      function Bcoeff_q2lcpp_s34(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Bcoeff_q2lcpp_s34
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2,t,ans 
!===== statement functions 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      Bcoeff_q2lcpp_s34=(-za(i1,i2)*zb(i4,i5)*
     -    (za(i1,i3)/(-s(i3,i4) + t(i1,i2,i6)) + 
     -      (za(i1,i5)*za(i4,i3)*zb(i5,i4))/
     -       (2.*(-s(i3,i4) + t(i1,i2,i6))**2)))/
     -  (za(i1,i6)*za(i2,i5)*za(i2,i6))


      return 
      end 

      function Bcoeff_q2lcpp_s16(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Bcoeff_q2lcpp_s16
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2,t,ans 
!===== statement functions 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      Bcoeff_q2lcpp_s16=   (-za(i3,i2)*zab2(i2,i3,i4,i5)*zb(i2,i6)*
     -    ((2*za(i1,i3))/(-s(i1,i6) + t(i1,i2,i6)) - 
     -      (za(i1,i6)*za(i2,i3)*zb(i6,i2))/
     -       (2.*(-s(i1,i6) + t(i1,i2,i6))**2)))/
     -  (za(i2,i5)*za(i2,i6)*za(i3,i4)*zab2(i5,i3,i4,i5))
      return 
      end 


!============ SUBLEADING COLORS + +  

      function Bcoeff_q2slcpp_s34
     &     (i1,i2,i3,i4,i5,i6,za,zb)
       
            include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2,t,ans 
!===== statement functions 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      Bcoeff_q2slcpp_s34=((za(i1,i3)/zab2(i5,i3,i4,i5) - 
     -       (za(i3,i4)*za(i5,i1)*zb(i4,i5))/
     -        (2.*zab2(i5,i3,i4,i5)**2))*zb(i5,i4))/
     -   (za(i2,i6)*za(i5,i6)) + 
     -  (za(i1,i5)**2*(za(i6,i3)/zab2(i5,i3,i4,i5) - 
     -       (za(i3,i4)*za(i5,i6)*zb(i4,i5))/
     -        (2.*zab2(i5,i3,i4,i5)**2))*zb(i5,i4))/
     -   (za(i1,i6)*za(i2,i5)*za(i5,i6)**2) + 
     -  (((za(i1,i5)*za(i3,i6))/
     -        (za(i5,i6)*zab2(i6,i3,i4,i6)) + 
     -       (za(i3,i4)*za(i6,i1)*zb(i4,i6))/
     -        zab2(i6,i3,i4,i6)**2)*zb(i6,i4))/
     -   (za(i2,i5)*za(i5,i6))

      return 
      end 

      function Bcoeff_q2slcpp_s16
     &     (i1,i2,i3,i4,i5,i6,za,zb)
       
            include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2,t,ans 
!===== statement functions 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

       Bcoeff_q2slcpp_s16=  (za(i3,i2)*(-(s(i2,i6)*za(i1,i6)*za(i2,i3))/
     -       (2.*zab2(i2,i1,i6,i2)**2) + 
     -      (za(i1,i3)*za(i2,i6) + za(i1,i2)*za(i3,i6))/
     -       zab2(i2,i1,i6,i2))*zab2(i2,i3,i4,i5)*zb(i2,i6))/
     -  (za(i2,i5)*za(i2,i6)**2*za(i3,i4)*zab2(i5,i3,i4,i5))

      return 
      end
      
      
      
      function Bcoeff_q2slcpp_s125
     &     (i1,i2,i3,i4,i5,i6,za,zb)
       
            include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2,t,ans 
!===== statement functions 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      Bcoeff_q2slcpp_s125= -((za(i1,i3)*za(i3,i5))/
     -     (za(i2,i5)*za(i3,i4)*za(i5,i6)**2)) + 
     -  (za(i1,i2)*za(i3,i5)**2*zb(i5,i2))/
     -   (za(i2,i5)*za(i3,i4)*za(i5,i6)**2*zab2(i5,i1,i2,i5))
     -   - (((za(i1,i5)*za(i3,i6))/
     -        (za(i5,i6)*zab2(i6,i3,i4,i6)) + 
     -       (za(i3,i4)*za(i6,i1)*zb(i4,i6))/
     -        zab2(i6,i3,i4,i6)**2)*zb(i6,i4))/
     -   (za(i2,i5)*za(i5,i6))

      return 
      end


      function Bcoeff_q2slcpp_s12
     &     (i1,i2,i3,i4,i5,i6,za,zb)
       
            include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2,t,ans 
!===== statement functions 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

        Bcoeff_q2slcpp_s12=-((za(i1,i2)*za(i3,i5)**2*zb(i5,i2))/
     -     (za(i2,i5)*za(i3,i4)*za(i5,i6)**2*
     -       zab2(i5,i1,i2,i5))) + 
     -  (za(i2,i1)*za(i3,i6)*zb(i6,i2)*
     -     (((za(i2,i3) - (za(i2,i6)*za(i3,i5))/za(i5,i6))*
     -          zab2(i6,i3,i4,i5))/
     -        (za(i2,i6)*zab2(i6,i1,i2,i6)) + 
     -       (t(i3,i4,i5)*za(i3,i6)*zb(i6,i5))/
     -        zab2(i6,i1,i2,i6)**2))/
     -   (za(i2,i6)*za(i3,i4)*za(i5,i6)*zab2(i5,i3,i4,i5))

        return 
        end
