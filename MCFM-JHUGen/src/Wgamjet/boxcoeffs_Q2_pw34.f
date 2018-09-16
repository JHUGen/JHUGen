!==== This file contains the box coefficients for the Q(i2) pieces of the following amplitude 

!===== (i1)^- (i2)^+ (i3)^- + (i4)^+ + gamma(i5) + gluon(i6) 

! where i3 and i4 are the leptons from a W decay. 

!=== C. Williams Nov 2013 

!================================================================================================      

!======= Q2 h5=1 h6=1 LEADING COLOR AMPLITUDES (order matches KC)
!===== two mass easy 
      function Dcoeff2me_q2lcpp_s26s34s126s134
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      ans=-(((-(s(i2,i6)*s(i3,i4)) + 
     -        t(i1,i2,i6)*t(i1,i3,i4))*
     -      za(i1,i2)*za(i1,i3)**2)/
     -    (2._dp*za(i1,i5)*za(i1,i6)*za(i2,i5)*
     -      za(i2,i6)*za(i3,i4)))

      
      Dcoeff2me_q2lcpp_s26s34s126s134=ans
      return 
      end

!==== one mass (s126)
      function Dcoeff1m_q2lcpp_s126
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      Dcoeff1m_q2lcpp_s126=-((s(i1,i6)*s(i2,i6)*za(i1,i3)**2*
     -      zab2(i2,i3,i4,i5))/
     -    (2._dp*za(i1,i6)*za(i2,i5)*za(i2,i6)*
     -      za(i3,i4)*zab2(i5,i3,i4,i5)))
      return 
      end


!==== one mass (s134)
      function Dcoeff1m_q2lcpp_s134
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      Dcoeff1m_q2lcpp_s134= -((s(i2,i5)*s(i2,i6)*za(i1,i3)**2)/
     -    (2._dp*za(i1,i6)*za(i2,i5)*za(i3,i4)*
     -      za(i5,i6)))
      return 
      end


!================== SUBLEADING COLOR ============================ 

!============ two mass s25 s134 

      function Dcoeff2me_q2slcpp_s25s34s125s134
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      ans= ((-(s(i2,i5)*s(i3,i4)) + t(i1,i2,i5)*t(i1,i3,i4))*
     -    za(i1,i3)**2)/
     -  (2.*za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6))
      
      Dcoeff2me_q2slcpp_s25s34s125s134=ans
      return 
      end


!======= two mass s26 s34 
      function Dcoeff2me_q2slcpp_s26s34s126s134
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      ans= -((-(s(i2,i6)*s(i3,i4)) + t(i1,i2,i6)*t(i1,i3,i4))*
     -     za(i1,i3)**2)/
     -  (2.*za(i1,i5)*za(i2,i6)*za(i3,i4)*za(i5,i6))
      
      Dcoeff2me_q2slcpp_s26s34s126s134=ans
      return 
      end


!===== two mass s12 s34 
      
      function Dcoeff2me_q2slcpp_s12s34s126s125
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      ans= ((-(s(i1,i2)*s(i3,i4)) + t(i1,i2,i5)*t(i1,i2,i6))*
     -    za(i1,i5)**2*za(i3,i6)**2)/
     -  (2.*za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6)**3)
      
      Dcoeff2me_q2slcpp_s12s34s126s125=ans
      return 
      end


!=====one mass s126 s12 s26
      function Dcoeff1m_q2slcpp_s126s12s26
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      Dcoeff1m_q2slcpp_s126s12s26=(za(i1,i3)**2*zab2(i6,i1,i2,i5))/
     -  (2.*za(i1,i6)*za(i2,i6)*za(i3,i4)*za(i5,i6)*
     -    zab2(i5,i3,i4,i5))*s(i1,i2)*s(i2,i6)

      return 
      end


!=====one mass s126 s16 s12
      function Dcoeff1m_q2slcpp_s126s16s12
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      Dcoeff1m_q2slcpp_s126s16s12=
     &     -(za(i1,i2)**2*za(i3,i6)**2*zab2(i2,i1,i6,i5))/
     -  (2.*za(i1,i6)*za(i2,i5)*za(i2,i6)**3*za(i3,i4)*
     -    zab2(i5,i3,i4,i5))*s(i1,i2)*s(i1,i6)

      return 
      end


!=====one mass s125 s12 s25
      function Dcoeff1m_q2slcpp_s125s12s25
     & (i1,i2,i3,i4,i5,i6,za,zb) 
       
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
      
      Dcoeff1m_q2slcpp_s125s12s25=(s(i1,i2)*s(i2,i5)*za(i1,i3)**2)/
     -  (2.*za(i1,i6)*za(i2,i5)*za(i3,i4)*za(i5,i6))
      return 
      end



