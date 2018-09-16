      subroutine real_aaajj_fill(j1,j2,j3,j4,j5,j6,j7,za,zb,real_aaajj)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      complex(dp)::real_aaajj(2,2,2,2,2,2), aaajj_aMHV, aaajj_j6MHV
      complex(dp)::aaajj_j7MHV, aaajj_aaNMHV, aaajj_aj6NMHV
      complex(dp)::aaajj_aj7NMHV, aaajj_jjNMHV

      real_aaajj(:,:,:,:,:,:)=czip

! real_aaajj is dimension 6 array 
! with entries real_aaajj(hel(j1),hel(j3),...)
! because helicity of quark line is conserved
      real_aaajj(2,1,2,2,2,2)=aaajj_aMHV(j1,j2,j3,j4,j5,j6,j7,za,zb)
      real_aaajj(2,2,1,2,2,2)=aaajj_aMHV(j1,j2,j4,j3,j5,j6,j7,za,zb)
      real_aaajj(2,2,2,1,2,2)=aaajj_aMHV(j1,j2,j5,j4,j3,j6,j7,za,zb)
      real_aaajj(2,1,1,2,2,2)=aaajj_aaNMHV(j1,j2,j3,j4,j5,j6,j7,za,zb)
      real_aaajj(2,1,2,1,2,2)=aaajj_aaNMHV(j1,j2,j3,j5,j4,j6,j7,za,zb)
      real_aaajj(2,2,1,1,2,2)=aaajj_aaNMHV(j1,j2,j5,j4,j3,j6,j7,za,zb)
      real_aaajj(2,1,1,1,2,2)=-aaajj_jjNMHV(j2,j1,j3,j4,j5,j7,j6,zb,za)

      real_aaajj(2,2,2,2,2,1)=aaajj_j7MHV(j1,j2,j3,j4,j5,j6,j7,za,zb)
      real_aaajj(2,1,2,2,2,1)=aaajj_aj7NMHV(j1,j2,j3,j4,j5,j6,j7,za,zb)
      real_aaajj(2,2,1,2,2,1)=aaajj_aj7NMHV(j1,j2,j4,j3,j5,j6,j7,za,zb)
      real_aaajj(2,2,2,1,2,1)=aaajj_aj7NMHV(j1,j2,j5,j4,j3,j6,j7,za,zb)
      real_aaajj(2,1,1,2,2,1)=-aaajj_aj7NMHV(j2,j1,j5,j4,j3,j7,j6,zb,za)
      real_aaajj(2,1,2,1,2,1)=-aaajj_aj7NMHV(j2,j1,j4,j5,j3,j7,j6,zb,za)
      real_aaajj(2,2,1,1,2,1)=-aaajj_aj7NMHV(j2,j1,j3,j4,j5,j7,j6,zb,za)
      real_aaajj(2,1,1,1,2,1)=-aaajj_j7MHV(j2,j1,j3,j4,j5,j7,j6,zb,za)

      real_aaajj(2,2,2,2,1,2)=aaajj_j6MHV(j1,j2,j3,j4,j5,j6,j7,za,zb)
      real_aaajj(2,1,2,2,1,2)=aaajj_aj6NMHV(j1,j2,j3,j4,j5,j6,j7,za,zb)
      real_aaajj(2,2,1,2,1,2)=aaajj_aj6NMHV(j1,j2,j4,j3,j5,j6,j7,za,zb)
      real_aaajj(2,2,2,1,1,2)=aaajj_aj6NMHV(j1,j2,j5,j4,j3,j6,j7,za,zb)
      real_aaajj(2,1,1,2,1,2)=-aaajj_aj6NMHV(j2,j1,j5,j4,j3,j7,j6,zb,za)
      real_aaajj(2,1,2,1,1,2)=-aaajj_aj6NMHV(j2,j1,j4,j5,j3,j7,j6,zb,za)
      real_aaajj(2,2,1,1,1,2)=-aaajj_aj6NMHV(j2,j1,j3,j4,j5,j7,j6,zb,za)
      real_aaajj(2,1,1,1,1,2)=-aaajj_j6MHV(j2,j1,j3,j4,j5,j7,j6,zb,za)

      real_aaajj(2,2,2,2,1,1)=aaajj_jjNMHV(j1,j2,j5,j4,j3,j6,j7,za,zb)
      real_aaajj(2,1,2,2,1,1)=-aaajj_aaNMHV(j2,j1,j5,j4,j3,j7,j6,zb,za)
      real_aaajj(2,2,1,2,1,1)=-aaajj_aaNMHV(j2,j1,j5,j3,j4,j7,j6,zb,za)
      real_aaajj(2,2,2,1,1,1)=-aaajj_aaNMHV(j2,j1,j3,j4,j5,j7,j6,zb,za)
      real_aaajj(2,1,1,2,1,1)=-aaajj_aMHV(j2,j1,j5,j4,j3,j7,j6,zb,za)
      real_aaajj(2,1,2,1,1,1)=-aaajj_aMHV(j2,j1,j4,j5,j3,j7,j6,zb,za)
      real_aaajj(2,2,1,1,1,1)=-aaajj_aMHV(j2,j1,j3,j4,j5,j7,j6,zb,za)



      real_aaajj(1,2,1,1,1,1)=aaajj_aMHV(j1,j2,j3,j4,j5,j6,j7,zb,za)
      real_aaajj(1,1,2,1,1,1)=aaajj_aMHV(j1,j2,j4,j3,j5,j6,j7,zb,za)
      real_aaajj(1,1,1,2,1,1)=aaajj_aMHV(j1,j2,j5,j4,j3,j6,j7,zb,za)
      real_aaajj(1,2,2,1,1,1)=aaajj_aaNMHV(j1,j2,j3,j4,j5,j6,j7,zb,za)
      real_aaajj(1,2,1,2,1,1)=aaajj_aaNMHV(j1,j2,j3,j5,j4,j6,j7,zb,za)
      real_aaajj(1,1,2,2,1,1)=aaajj_aaNMHV(j1,j2,j5,j4,j3,j6,j7,zb,za)
      real_aaajj(1,2,2,2,1,1)=-aaajj_jjNMHV(j2,j1,j3,j4,j5,j7,j6,za,zb)

      real_aaajj(1,1,1,1,1,2)=aaajj_j7MHV(j1,j2,j3,j4,j5,j6,j7,zb,za)
      real_aaajj(1,2,1,1,1,2)=aaajj_aj7NMHV(j1,j2,j3,j4,j5,j6,j7,zb,za)
      real_aaajj(1,1,2,1,1,2)=aaajj_aj7NMHV(j1,j2,j4,j3,j5,j6,j7,zb,za)
      real_aaajj(1,1,1,2,1,2)=aaajj_aj7NMHV(j1,j2,j5,j4,j3,j6,j7,zb,za)
      real_aaajj(1,2,2,1,1,2)=-aaajj_aj7NMHV(j2,j1,j5,j4,j3,j7,j6,za,zb)
      real_aaajj(1,2,1,2,1,2)=-aaajj_aj7NMHV(j2,j1,j4,j5,j3,j7,j6,za,zb)
      real_aaajj(1,1,2,2,1,2)=-aaajj_aj7NMHV(j2,j1,j3,j4,j5,j7,j6,za,zb)
      real_aaajj(1,2,2,2,1,2)=-aaajj_j7MHV(j2,j1,j3,j4,j5,j7,j6,za,zb)

      real_aaajj(1,1,1,1,2,1)=aaajj_j6MHV(j1,j2,j3,j4,j5,j6,j7,zb,za)
      real_aaajj(1,2,1,1,2,1)=aaajj_aj6NMHV(j1,j2,j3,j4,j5,j6,j7,zb,za)
      real_aaajj(1,1,2,1,2,1)=aaajj_aj6NMHV(j1,j2,j4,j3,j5,j6,j7,zb,za)
      real_aaajj(1,1,1,2,2,1)=aaajj_aj6NMHV(j1,j2,j5,j4,j3,j6,j7,zb,za)
      real_aaajj(1,2,2,1,2,1)=-aaajj_aj6NMHV(j2,j1,j5,j4,j3,j7,j6,za,zb)
      real_aaajj(1,2,1,2,2,1)=-aaajj_aj6NMHV(j2,j1,j4,j5,j3,j7,j6,za,zb)
      real_aaajj(1,1,2,2,2,1)=-aaajj_aj6NMHV(j2,j1,j3,j4,j5,j7,j6,za,zb)
      real_aaajj(1,2,2,2,2,1)=-aaajj_j6MHV(j2,j1,j3,j4,j5,j7,j6,za,zb)

      real_aaajj(1,1,1,1,2,2)=aaajj_jjNMHV(j1,j2,j5,j4,j3,j6,j7,zb,za)
      real_aaajj(1,2,1,1,2,2)=-aaajj_aaNMHV(j2,j1,j5,j4,j3,j7,j6,za,zb)
      real_aaajj(1,1,2,1,2,2)=-aaajj_aaNMHV(j2,j1,j5,j3,j4,j7,j6,za,zb)
      real_aaajj(1,1,1,2,2,2)=-aaajj_aaNMHV(j2,j1,j3,j4,j5,j7,j6,za,zb)
      real_aaajj(1,2,2,1,2,2)=-aaajj_aMHV(j2,j1,j5,j4,j3,j7,j6,za,zb)
      real_aaajj(1,2,1,2,2,2)=-aaajj_aMHV(j2,j1,j4,j5,j3,j7,j6,za,zb)
      real_aaajj(1,1,2,2,2,2)=-aaajj_aMHV(j2,j1,j3,j4,j5,j7,j6,za,zb)



      return
      end
      
      function aaajj_aMHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_aMHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      aaajj_aMHV = (za(i1,i2)**2*za(i2,i3)**2)/
     & (za(i1,i4)*za(i1,i5)*za(i1,i7)*za(i2,i4)*za(i2,i5)*za(i2,i6)*
     &   za(i6,i7))
      return
      end

      function aaajj_j6MHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_j6MHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      aaajj_j6MHV = (za(i1,i2)**2*za(i1,i6)*za(i2,i6)**2)/
     & (za(i1,i3)*za(i1,i4)*za(i1,i5)*za(i1,i7)*za(i2,i3)*za(i2,i4)*
     &   za(i2,i5)*za(i6,i7))
      return
      end

      function aaajj_j7MHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_j7MHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      aaajj_j7MHV = (za(i1,i2)**2*za(i2,i7)**3)/
     & (za(i1,i3)*za(i1,i4)*za(i1,i5)*za(i2,i3)*za(i2,i4)*za(i2,i5)*
     &   za(i2,i6)*za(i6,i7))
      return
      end

      function aaajj_aaNMHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_aaNMHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      real(dp):: t
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)


      aaajj_aaNMHV =
     & -((za(i2,i3)**2*(za(i1,i2)*zb(i6,i2) + za(i1,i3)*zb(i6,i3))*
     &      (za(i2,i4)*zb(i6,i2) + za(i3,i4)*zb(i6,i3))**2)/
     &    (t(i2,i3,i6)*za(i1,i5)*za(i1,i7)*
     &      (za(i2,i5)*zb(i6,i2) + za(i3,i5)*zb(i6,i3))*
     &      (t(i2,i3,i6) - za(i2,i6)*zb(i6,i2) - za(i3,i6)*zb(i6,i3))*
     &      (t(i2,i3,i6)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &        za(i2,i6)*za(i3,i7)*zb(i6,i3)))) - 
     & (za(i2,i4)**2*(za(i1,i2)*zb(i6,i2) + za(i1,i4)*zb(i6,i4))*
     &    (za(i2,i3)*zb(i6,i2) - za(i3,i4)*zb(i6,i4))**2)/
     &  (t(i2,i4,i6)*za(i1,i5)*za(i1,i7)*
     &    (za(i2,i5)*zb(i6,i2) + za(i4,i5)*zb(i6,i4))*
     &    (t(i2,i4,i6) - za(i2,i6)*zb(i6,i2) - za(i4,i6)*zb(i6,i4))*
     &    (t(i2,i4,i6)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &      za(i2,i6)*za(i4,i7)*zb(i6,i4))) + 
     & (za(i2,i3)**2*(za(i2,i3)*zb(i6,i3) + za(i2,i5)*zb(i6,i5))*
     &    (za(i2,i4)*zb(i6,i2) + za(i3,i4)*zb(i6,i3) - 
     &       za(i4,i5)*zb(i6,i5))**2)/
     &  (t(i1,i4,i7)*za(i1,i7)*za(i2,i5)*
     &    (za(i2,i5)*zb(i6,i2) + za(i3,i5)*zb(i6,i3))*
     &    (t(i1,i4,i7) - za(i2,i6)*zb(i6,i2) - za(i3,i6)*zb(i6,i3) - 
     &      za(i5,i6)*zb(i6,i5))*
     &    (t(i1,i4,i7)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &      za(i2,i6)*za(i3,i7)*zb(i6,i3) - 
     &      za(i2,i6)*za(i5,i7)*zb(i6,i5))) + 
     & (za(i2,i4)**2*(za(i2,i4)*zb(i6,i4) + za(i2,i5)*zb(i6,i5))*
     &    (za(i2,i3)*zb(i6,i2) - za(i3,i4)*zb(i6,i4) - 
     &       za(i3,i5)*zb(i6,i5))**2)/
     &  (t(i1,i3,i7)*za(i1,i7)*za(i2,i5)*
     &    (za(i2,i5)*zb(i6,i2) + za(i4,i5)*zb(i6,i4))*
     &    (t(i1,i3,i7) - za(i2,i6)*zb(i6,i2) - za(i4,i6)*zb(i6,i4) - 
     &      za(i5,i6)*zb(i6,i5))*
     &    (t(i1,i3,i7)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &      za(i2,i6)*za(i4,i7)*zb(i6,i4) - 
     &      za(i2,i6)*za(i5,i7)*zb(i6,i5))) - 
     & (za(i1,i2)*(za(i2,i7)*zb(i2,i1) + za(i6,i7)*zb(i6,i1))*zb(i6,i2)*
     &    (s(i3,i4)*za(i2,i7) + za(i2,i3)*za(i2,i7)*zb(i3,i2) + 
     &      za(i2,i4)*za(i2,i7)*zb(i4,i2) - 
     &      za(i2,i3)*za(i6,i7)*zb(i6,i3) - 
     &      za(i2,i4)*za(i6,i7)*zb(i6,i4))*
     &    (za(i2,i5)*zb(i5,i1) + za(i2,i6)*zb(i6,i1) + 
     &      za(i2,i7)*zb(i7,i1)))/
     &  (s(i2,i6)*za(i1,i5)*za(i1,i7)*za(i2,i5)*za(i6,i7)*zb(i3,i1)*
     &    zb(i4,i1)*(za(i2,i7)*zb(i3,i2) - za(i6,i7)*zb(i6,i3))*
     &    (za(i2,i7)*zb(i4,i2) - za(i6,i7)*zb(i6,i4))) + 
     & (za(i2,i4)*za(i2,i7)*zb(i6,i2)*
     &    (za(i2,i6)*zb(i6,i1) + za(i2,i7)*zb(i7,i1))*
     &    (za(i1,i3)*za(i2,i7)*zb(i5,i1) + 
     &      za(i2,i6)*za(i3,i7)*zb(i6,i5) + 
     &      za(i2,i7)*za(i3,i7)*zb(i7,i5)))/
     &  (s(i2,i6)*za(i1,i7)*za(i2,i5)*za(i6,i7)*zb(i3,i1)*
     &    (za(i2,i7)*zb(i4,i2) - za(i6,i7)*zb(i6,i4))*
     &    (s(i4,i5)*za(i2,i7) + za(i2,i4)*za(i2,i7)*zb(i4,i2) + 
     &      za(i2,i5)*za(i2,i7)*zb(i5,i2) - 
     &      za(i2,i4)*za(i6,i7)*zb(i6,i4) - 
     &      za(i2,i5)*za(i6,i7)*zb(i6,i5))) + 
     & (za(i2,i3)*za(i2,i7)*zb(i6,i2)*
     &    (za(i2,i6)*zb(i6,i1) + za(i2,i7)*zb(i7,i1))*
     &    (za(i1,i4)*za(i2,i7)*zb(i5,i1) + 
     &      za(i2,i6)*za(i4,i7)*zb(i6,i5) + 
     &      za(i2,i7)*za(i4,i7)*zb(i7,i5)))/
     &  (s(i2,i6)*za(i1,i7)*za(i2,i5)*za(i6,i7)*zb(i4,i1)*
     &    (za(i2,i7)*zb(i3,i2) - za(i6,i7)*zb(i6,i3))*
     &    (s(i3,i5)*za(i2,i7) + za(i2,i3)*za(i2,i7)*zb(i3,i2) + 
     &      za(i2,i5)*za(i2,i7)*zb(i5,i2) - 
     &      za(i2,i3)*za(i6,i7)*zb(i6,i3) - 
     &      za(i2,i5)*za(i6,i7)*zb(i6,i5))) + 
     & (za(i2,i3)*za(i2,i7)*zb(i5,i1)*zb(i6,i2)*
     &    (za(i1,i4)*za(i2,i6)*zb(i6,i1) - 
     &      za(i2,i6)*za(i4,i5)*zb(i6,i5) + 
     &      za(i1,i4)*za(i2,i7)*zb(i7,i1) - 
     &      za(i2,i7)*za(i4,i5)*zb(i7,i5)))/
     &  (s(i2,i6)*za(i1,i5)*za(i6,i7)*zb(i4,i1)*
     &    (za(i2,i7)*zb(i3,i2) - za(i6,i7)*zb(i6,i3))*
     &    (za(i2,i3)*za(i2,i7)*zb(i3,i2) + 
     &      za(i2,i6)*za(i2,i7)*zb(i6,i2) + 
     &      za(i2,i6)*za(i3,i7)*zb(i6,i3) - 
     &      za(i2,i3)*za(i6,i7)*zb(i6,i3) + za(i2,i7)**2*zb(i7,i2) + 
     &      za(i2,i7)*za(i3,i7)*zb(i7,i3) + 
     &      za(i2,i7)*za(i6,i7)*zb(i7,i6))) + 
     & (za(i2,i4)*za(i2,i7)*zb(i5,i1)*zb(i6,i2)*
     &    (za(i1,i3)*za(i2,i6)*zb(i6,i1) - 
     &      za(i2,i6)*za(i3,i5)*zb(i6,i5) + 
     &      za(i1,i3)*za(i2,i7)*zb(i7,i1) - 
     &      za(i2,i7)*za(i3,i5)*zb(i7,i5)))/
     &  (s(i2,i6)*za(i1,i5)*za(i6,i7)*zb(i3,i1)*
     &    (za(i2,i7)*zb(i4,i2) - za(i6,i7)*zb(i6,i4))*
     &    (za(i2,i4)*za(i2,i7)*zb(i4,i2) + 
     &      za(i2,i6)*za(i2,i7)*zb(i6,i2) + 
     &      za(i2,i6)*za(i4,i7)*zb(i6,i4) - 
     &      za(i2,i4)*za(i6,i7)*zb(i6,i4) + za(i2,i7)**2*zb(i7,i2) + 
     &      za(i2,i7)*za(i4,i7)*zb(i7,i4) + 
     &      za(i2,i7)*za(i6,i7)*zb(i7,i6)))
      return
      end

      function aaajj_aj6NMHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_aj6NMHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      real(dp):: t
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

      aaajj_aj6NMHV = (za(i1,i3)**2*za(i2,i6)**2*zb(i7,i1)**2*
     &    (za(i1,i2)*zb(i7,i1) - za(i2,i3)*zb(i7,i3))*
     &    (za(i1,i6)*zb(i7,i1) + za(i3,i6)*zb(i7,i3)))/
     &  (t(i1,i3,i7)*za(i2,i4)*za(i2,i5)*
     &    (za(i1,i4)*zb(i7,i1) + za(i3,i4)*zb(i7,i3))*
     &    (za(i1,i5)*zb(i7,i1) + za(i3,i5)*zb(i7,i3))*
     &    (-(t(i1,i3,i7)*za(i1,i6)) + za(i1,i6)*za(i1,i7)*zb(i7,i1) + 
     &      za(i1,i7)*za(i3,i6)*zb(i7,i3))*
     &    (t(i1,i3,i7) - za(i1,i7)*zb(i7,i1) - za(i3,i7)*zb(i7,i3))) - 
     & (za(i1,i6)**2*za(i2,i6)*zb(i7,i1)*
     &    (za(i1,i6)*zb(i4,i1) + za(i6,i7)*zb(i7,i4))*
     &    (za(i1,i6)*za(i2,i3)*zb(i5,i2)*zb(i7,i1) + 
     &      za(i1,i6)*za(i3,i6)*zb(i6,i5)*zb(i7,i1) + 
     &      s(i1,i7)*za(i3,i6)*zb(i7,i5)))/
     &  (s(i1,i7)*za(i1,i4)*za(i2,i5)*za(i6,i7)*
     &    (-(za(i1,i6)*zb(i6,i2)*zb(i7,i1)) - s(i1,i7)*zb(i7,i2))*
     &    (za(i1,i6)*zb(i3,i1) + za(i6,i7)*zb(i7,i3))*
     &    (-(s(i3,i4)*za(i1,i6)) - za(i1,i3)*za(i1,i6)*zb(i3,i1) - 
     &      za(i1,i4)*za(i1,i6)*zb(i4,i1) - 
     &      za(i1,i3)*za(i6,i7)*zb(i7,i3) - 
     &      za(i1,i4)*za(i6,i7)*zb(i7,i4))) - 
     & (za(i1,i6)**2*za(i2,i6)*zb(i7,i1)*
     &    (za(i1,i6)*za(i2,i3)*zb(i4,i2)*zb(i7,i1) + 
     &      za(i1,i6)*za(i3,i6)*zb(i6,i4)*zb(i7,i1) + 
     &      s(i1,i7)*za(i3,i6)*zb(i7,i4))*
     &    (za(i1,i6)*zb(i5,i1) + za(i6,i7)*zb(i7,i5)))/
     &  (s(i1,i7)*za(i1,i5)*za(i2,i4)*za(i6,i7)*
     &    (-(za(i1,i6)*zb(i6,i2)*zb(i7,i1)) - s(i1,i7)*zb(i7,i2))*
     &    (za(i1,i6)*zb(i3,i1) + za(i6,i7)*zb(i7,i3))*
     &    (-(s(i3,i5)*za(i1,i6)) - za(i1,i3)*za(i1,i6)*zb(i3,i1) - 
     &      za(i1,i5)*za(i1,i6)*zb(i5,i1) - 
     &      za(i1,i3)*za(i6,i7)*zb(i7,i3) - 
     &      za(i1,i5)*za(i6,i7)*zb(i7,i5))) - 
     & (za(i2,i6)**2*zb(i7,i2)*
     &    (za(i1,i2)*zb(i7,i2) + za(i1,i6)*zb(i7,i6))**2*
     &    (za(i2,i3)*zb(i7,i2) - za(i3,i6)*zb(i7,i6))**2)/
     &  (s(i2,i6)*za(i1,i4)*za(i1,i5)*
     &    (-(s(i2,i6)*za(i1,i6)) - za(i1,i7)*za(i2,i6)*zb(i7,i2))*
     &    zb(i7,i6)*(za(i2,i4)*zb(i7,i2) - za(i4,i6)*zb(i7,i6))*
     &    (za(i2,i5)*zb(i7,i2) - za(i5,i6)*zb(i7,i6))*
     &    (s(i2,i6) + za(i2,i7)*zb(i7,i2) + za(i6,i7)*zb(i7,i6))) - 
     & (za(i2,i6)**2*zb(i4,i2)**2*
     &    (za(i2,i6)*zb(i7,i2) + za(i4,i6)*zb(i7,i4))*
     &    (za(i1,i2)*zb(i7,i2) + za(i1,i4)*zb(i7,i4) + 
     &      za(i1,i6)*zb(i7,i6))*
     &    (za(i2,i3)*zb(i7,i2) - za(i3,i4)*zb(i7,i4) - 
     &       za(i3,i6)*zb(i7,i6))**2)/
     &  (t(i2,i4,i6)*za(i1,i5)*
     &    (t(i2,i4,i6) - za(i2,i6)*zb(i6,i2) - za(i4,i6)*zb(i6,i4))*
     &    (-(t(i2,i4,i6)*za(i1,i6)) - za(i1,i7)*za(i2,i6)*zb(i7,i2) - 
     &      za(i1,i7)*za(i4,i6)*zb(i7,i4))*
     &    (-(t(i2,i4,i6)*zb(i7,i2)) + za(i2,i6)*zb(i6,i2)*zb(i7,i2) + 
     &      za(i4,i6)*zb(i6,i2)*zb(i7,i4))*
     &    (za(i2,i5)*zb(i7,i2) + za(i4,i5)*zb(i7,i4) - 
     &      za(i5,i6)*zb(i7,i6))*
     &    (t(i2,i4,i6) + za(i2,i7)*zb(i7,i2) + za(i4,i7)*zb(i7,i4) + 
     &      za(i6,i7)*zb(i7,i6))) - 
     & (za(i2,i6)**2*zb(i5,i2)**2*
     &    (za(i2,i6)*zb(i7,i2) + za(i5,i6)*zb(i7,i5))*
     &    (za(i1,i2)*zb(i7,i2) + za(i1,i5)*zb(i7,i5) + 
     &      za(i1,i6)*zb(i7,i6))*
     &    (za(i2,i3)*zb(i7,i2) - za(i3,i5)*zb(i7,i5) - 
     &       za(i3,i6)*zb(i7,i6))**2)/
     &  (t(i2,i5,i6)*za(i1,i4)*
     &    (t(i2,i5,i6) - za(i2,i6)*zb(i6,i2) - za(i5,i6)*zb(i6,i5))*
     &    (-(t(i2,i5,i6)*za(i1,i6)) - za(i1,i7)*za(i2,i6)*zb(i7,i2) - 
     &      za(i1,i7)*za(i5,i6)*zb(i7,i5))*
     &    (-(t(i2,i5,i6)*zb(i7,i2)) + za(i2,i6)*zb(i6,i2)*zb(i7,i2) + 
     &      za(i5,i6)*zb(i6,i2)*zb(i7,i5))*
     &    (za(i2,i4)*zb(i7,i2) - za(i4,i5)*zb(i7,i5) - 
     &      za(i4,i6)*zb(i7,i6))*
     &    (t(i2,i5,i6) + za(i2,i7)*zb(i7,i2) + za(i5,i7)*zb(i7,i5) + 
     &      za(i6,i7)*zb(i7,i6))) - 
     & (za(i1,i2)*za(i1,i6)*zb(i7,i1)**2*
     &    (za(i1,i6)*zb(i2,i1) + za(i6,i7)*zb(i7,i2))*
     &    (-(s(i4,i5)*za(i1,i6)) - za(i1,i4)*za(i1,i6)*zb(i4,i1) - 
     &      za(i1,i5)*za(i1,i6)*zb(i5,i1) - 
     &      za(i1,i4)*za(i6,i7)*zb(i7,i4) - 
     &      za(i1,i5)*za(i6,i7)*zb(i7,i5))*
     &    (-(s(i1,i7)*za(i2,i6)) - za(i1,i6)*za(i2,i3)*zb(i3,i1) - 
     &      za(i1,i6)*za(i2,i6)*zb(i6,i1) - 
     &      za(i2,i3)*za(i6,i7)*zb(i7,i3) - 
     &      za(i2,i6)*za(i6,i7)*zb(i7,i6)))/
     &  (s(i1,i7)*za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i6,i7)*
     &    zb(i3,i2)*(-(za(i1,i6)*zb(i6,i2)*zb(i7,i1)) - 
     &      s(i1,i7)*zb(i7,i2))*
     &    (za(i1,i6)*zb(i3,i1) + za(i6,i7)*zb(i7,i3))*
     &    (s(i1,i7) + za(i1,i6)*zb(i6,i1) + za(i6,i7)*zb(i7,i6))) - 
     & (za(i1,i6)**2*za(i2,i3)*
     &    (za(i2,i6)*zb(i5,i2) + za(i3,i6)*zb(i5,i3))*zb(i7,i1)**2*
     &    (za(i1,i6)*zb(i4,i1) + za(i6,i7)*zb(i7,i4)))/
     &  (s(i1,i7)*za(i1,i4)*za(i2,i5)*za(i6,i7)*zb(i3,i2)*
     &    (s(i1,i7) + za(i1,i6)*zb(i6,i1) + za(i6,i7)*zb(i7,i6))*
     &    (s(i1,i7)*za(i1,i6)*zb(i7,i1) + 
     &      za(i1,i4)*za(i1,i6)*zb(i4,i1)*zb(i7,i1) + 
     &      za(i1,i6)**2*zb(i6,i1)*zb(i7,i1) + 
     &      za(i1,i6)*za(i4,i6)*zb(i6,i4)*zb(i7,i1) + 
     &      s(i1,i7)*za(i4,i6)*zb(i7,i4) + 
     &      za(i1,i4)*za(i6,i7)*zb(i7,i1)*zb(i7,i4) + 
     &      za(i1,i6)*za(i6,i7)*zb(i7,i1)*zb(i7,i6))) - 
     & (za(i1,i6)**2*za(i2,i3)*
     &    (za(i2,i6)*zb(i4,i2) + za(i3,i6)*zb(i4,i3))*zb(i7,i1)**2*
     &    (za(i1,i6)*zb(i5,i1) + za(i6,i7)*zb(i7,i5)))/
     &  (s(i1,i7)*za(i1,i5)*za(i2,i4)*za(i6,i7)*zb(i3,i2)*
     &    (s(i1,i7) + za(i1,i6)*zb(i6,i1) + za(i6,i7)*zb(i7,i6))*
     &    (s(i1,i7)*za(i1,i6)*zb(i7,i1) + 
     &      za(i1,i5)*za(i1,i6)*zb(i5,i1)*zb(i7,i1) + 
     &      za(i1,i6)**2*zb(i6,i1)*zb(i7,i1) + 
     &      za(i1,i6)*za(i5,i6)*zb(i6,i5)*zb(i7,i1) + 
     &      s(i1,i7)*za(i5,i6)*zb(i7,i5) + 
     &      za(i1,i5)*za(i6,i7)*zb(i7,i1)*zb(i7,i5) + 
     &      za(i1,i6)*za(i6,i7)*zb(i7,i1)*zb(i7,i6)))
      return
      end

      function aaajj_aj7NMHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_aj7NMHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      real(dp):: t
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

      aaajj_aj7NMHV =
     &  -((za(i2,i3)**2*(za(i1,i2)*zb(i6,i2) + za(i1,i3)*zb(i6,i3))*
     &      (za(i2,i7)*zb(i6,i2) + za(i3,i7)*zb(i6,i3))**3)/
     &    (t(i2,i3,i6)*za(i1,i4)*za(i1,i5)*
     &      (za(i2,i4)*zb(i6,i2) + za(i3,i4)*zb(i6,i3))*
     &      (za(i2,i5)*zb(i6,i2) + za(i3,i5)*zb(i6,i3))*
     &      (t(i2,i3,i6) - za(i2,i6)*zb(i6,i2) - za(i3,i6)*zb(i6,i3))*
     &      (t(i2,i3,i6)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &        za(i2,i6)*za(i3,i7)*zb(i6,i3)))) + 
     & (za(i2,i3)**2*(za(i2,i3)*zb(i6,i3) + za(i2,i4)*zb(i6,i4))*
     &    (za(i2,i7)*zb(i6,i2) + za(i3,i7)*zb(i6,i3) + 
     &       za(i4,i7)*zb(i6,i4))**3)/
     &  (t(i1,i5,i7)*za(i1,i5)*za(i2,i4)*
     &    (za(i2,i4)*zb(i6,i2) + za(i3,i4)*zb(i6,i3))*
     &    (za(i2,i5)*zb(i6,i2) + za(i3,i5)*zb(i6,i3) + 
     &      za(i4,i5)*zb(i6,i4))*
     &    (t(i1,i5,i7) - za(i2,i6)*zb(i6,i2) - za(i3,i6)*zb(i6,i3) - 
     &      za(i4,i6)*zb(i6,i4))*
     &    (t(i1,i5,i7)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &      za(i2,i6)*za(i3,i7)*zb(i6,i3) - 
     &      za(i2,i6)*za(i4,i7)*zb(i6,i4))) + 
     & (za(i2,i3)**2*(za(i2,i3)*zb(i6,i3) + za(i2,i5)*zb(i6,i5))*
     &    (za(i2,i7)*zb(i6,i2) + za(i3,i7)*zb(i6,i3) + 
     &       za(i5,i7)*zb(i6,i5))**3)/
     &  (t(i1,i4,i7)*za(i1,i4)*za(i2,i5)*
     &    (za(i2,i5)*zb(i6,i2) + za(i3,i5)*zb(i6,i3))*
     &    (za(i2,i4)*zb(i6,i2) + za(i3,i4)*zb(i6,i3) - 
     &      za(i4,i5)*zb(i6,i5))*
     &    (t(i1,i4,i7) - za(i2,i6)*zb(i6,i2) - za(i3,i6)*zb(i6,i3) - 
     &      za(i5,i6)*zb(i6,i5))*
     &    (t(i1,i4,i7)*za(i2,i7) - za(i2,i6)*za(i2,i7)*zb(i6,i2) - 
     &      za(i2,i6)*za(i3,i7)*zb(i6,i3) - 
     &      za(i2,i6)*za(i5,i7)*zb(i6,i5))) + 
     & (za(i2,i3)*za(i2,i7)**4*zb(i5,i1)*
     &    (za(i1,i7)*zb(i4,i1) - za(i5,i7)*zb(i5,i4))*zb(i6,i2)**2)/
     &  (s(i2,i6)*za(i1,i5)*za(i2,i4)*za(i6,i7)*
     &    (za(i2,i7)*zb(i3,i2) - za(i6,i7)*zb(i6,i3))*
     &    (s(i3,i4)*za(i2,i7) + za(i2,i3)*za(i2,i7)*zb(i3,i2) + 
     &      za(i2,i4)*za(i2,i7)*zb(i4,i2) - 
     &      za(i2,i3)*za(i6,i7)*zb(i6,i3) - 
     &      za(i2,i4)*za(i6,i7)*zb(i6,i4))*
     &    (s(i2,i6)*zb(i6,i1) + za(i2,i7)*zb(i6,i2)*zb(i7,i1))) + 
     & (za(i2,i3)*za(i2,i7)**4*zb(i4,i1)*
     &    (za(i1,i7)*zb(i5,i1) + za(i4,i7)*zb(i5,i4))*zb(i6,i2)**2)/
     &  (s(i2,i6)*za(i1,i4)*za(i2,i5)*za(i6,i7)*
     &    (za(i2,i7)*zb(i3,i2) - za(i6,i7)*zb(i6,i3))*
     &    (s(i3,i5)*za(i2,i7) + za(i2,i3)*za(i2,i7)*zb(i3,i2) + 
     &      za(i2,i5)*za(i2,i7)*zb(i5,i2) - 
     &      za(i2,i3)*za(i6,i7)*zb(i6,i3) - 
     &      za(i2,i5)*za(i6,i7)*zb(i6,i5))*
     &    (s(i2,i6)*zb(i6,i1) + za(i2,i7)*zb(i6,i2)*zb(i7,i1))) + 
     & (za(i1,i7)**2*za(i2,i3)**2*zb(i6,i1)**3*
     &    (za(i1,i2)*zb(i6,i1) + za(i2,i7)*zb(i7,i6))**2)/
     &  (s(i1,i7)*za(i2,i4)*za(i2,i5)*
     &    (s(i1,i7)*za(i2,i7) + za(i1,i7)*za(i2,i6)*zb(i6,i1))*
     &    zb(i7,i6)*(za(i1,i4)*zb(i6,i1) + za(i4,i7)*zb(i7,i6))*
     &    (za(i1,i5)*zb(i6,i1) + za(i5,i7)*zb(i7,i6))*
     &    (s(i1,i7) + za(i1,i6)*zb(i6,i1) + za(i6,i7)*zb(i7,i6))) - 
     & (za(i1,i2)*za(i2,i7)**2*
     &    (za(i2,i4)*zb(i4,i1) + za(i2,i5)*zb(i5,i1))*
     &    (za(i2,i7)*zb(i2,i1) + za(i6,i7)*zb(i6,i1))*zb(i6,i2)*
     &    (s(i2,i6)*za(i2,i7)*zb(i6,i2) + 
     &      za(i2,i3)*za(i2,i7)*zb(i3,i2)*zb(i6,i2) + 
     &      s(i2,i6)*za(i3,i7)*zb(i6,i3) - 
     &      za(i2,i3)*za(i6,i7)*zb(i6,i2)*zb(i6,i3) + 
     &      za(i2,i7)**2*zb(i6,i2)*zb(i7,i2) + 
     &      za(i2,i7)*za(i3,i7)*zb(i6,i2)*zb(i7,i3) + 
     &      za(i2,i7)*za(i6,i7)*zb(i6,i2)*zb(i7,i6)))/
     &  (s(i2,i6)*za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i6,i7)*
     &    zb(i3,i1)*(za(i2,i7)*zb(i3,i2) - za(i6,i7)*zb(i6,i3))*
     &    (s(i2,i6)*zb(i6,i1) + za(i2,i7)*zb(i6,i2)*zb(i7,i1))*
     &    (s(i2,i6) + za(i2,i7)*zb(i7,i2) + za(i6,i7)*zb(i7,i6))) + 
     & (za(i2,i7)**4*zb(i5,i1)*
     &    (za(i1,i3)*zb(i4,i1) + za(i3,i5)*zb(i5,i4))*zb(i6,i2)**2)/
     &  (s(i2,i6)*za(i1,i5)*za(i2,i4)*za(i6,i7)*zb(i3,i1)*
     &    (s(i2,i6) + za(i2,i7)*zb(i7,i2) + za(i6,i7)*zb(i7,i6))*
     &    (s(i2,i6)*za(i2,i7)*zb(i6,i2) + 
     &      za(i2,i4)*za(i2,i7)*zb(i4,i2)*zb(i6,i2) + 
     &      s(i2,i6)*za(i4,i7)*zb(i6,i4) - 
     &      za(i2,i4)*za(i6,i7)*zb(i6,i2)*zb(i6,i4) + 
     &      za(i2,i7)**2*zb(i6,i2)*zb(i7,i2) + 
     &      za(i2,i7)*za(i4,i7)*zb(i6,i2)*zb(i7,i4) + 
     &      za(i2,i7)*za(i6,i7)*zb(i6,i2)*zb(i7,i6))) + 
     & (za(i2,i7)**4*zb(i4,i1)*
     &    (za(i1,i3)*zb(i5,i1) - za(i3,i4)*zb(i5,i4))*zb(i6,i2)**2)/
     &  (s(i2,i6)*za(i1,i4)*za(i2,i5)*za(i6,i7)*zb(i3,i1)*
     &    (s(i2,i6) + za(i2,i7)*zb(i7,i2) + za(i6,i7)*zb(i7,i6))*
     &    (s(i2,i6)*za(i2,i7)*zb(i6,i2) + 
     &      za(i2,i5)*za(i2,i7)*zb(i5,i2)*zb(i6,i2) + 
     &      s(i2,i6)*za(i5,i7)*zb(i6,i5) - 
     &      za(i2,i5)*za(i6,i7)*zb(i6,i2)*zb(i6,i5) + 
     &      za(i2,i7)**2*zb(i6,i2)*zb(i7,i2) + 
     &      za(i2,i7)*za(i5,i7)*zb(i6,i2)*zb(i7,i5) + 
     &      za(i2,i7)*za(i6,i7)*zb(i6,i2)*zb(i7,i6)))
      return
      end

      function aaajj_jjNMHV(i1,i2,i3,i4,i5,i6,i7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaajj_jjNMHV
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6,i7
      real(dp):: t
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      aaajj_jjNMHV = -(((za(i1,i2)*zb(i2,i6) - za(i3,i1)*zb(i3,i6))*
     &      (t(i2,i3,i6)*za(i7,i2) - za(i6,i2)*za(i7,i2)*zb(i2,i6) - 
     &         za(i6,i2)*za(i7,i3)*zb(i3,i6))**2)/
     &    (t(i2,i3,i6)*za(i3,i2)**2*za(i4,i1)*za(i5,i1)*zb(i2,i6)*
     &      (za(i4,i2)*zb(i2,i6) + za(i4,i3)*zb(i3,i6))*
     &      (za(i5,i2)*zb(i2,i6) + za(i5,i3)*zb(i3,i6)))) - 
     & ((za(i1,i2)*zb(i2,i6) - za(i4,i1)*zb(i4,i6))*
     &    (t(i2,i4,i6)*za(i7,i2) - za(i6,i2)*za(i7,i2)*zb(i2,i6) - 
     &       za(i6,i2)*za(i7,i4)*zb(i4,i6))**2)/
     &  (t(i2,i4,i6)*za(i3,i1)*za(i4,i2)**2*za(i5,i1)*zb(i2,i6)*
     &    (za(i3,i2)*zb(i2,i6) - za(i4,i3)*zb(i4,i6))*
     &    (za(i5,i2)*zb(i2,i6) + za(i5,i4)*zb(i4,i6))) - 
     & ((za(i1,i2)*zb(i2,i6) - za(i5,i1)*zb(i5,i6))*
     &    (t(i2,i5,i6)*za(i7,i2) - za(i6,i2)*za(i7,i2)*zb(i2,i6) - 
     &       za(i6,i2)*za(i7,i5)*zb(i5,i6))**2)/
     &  (t(i2,i5,i6)*za(i3,i1)*za(i4,i1)*za(i5,i2)**2*zb(i2,i6)*
     &    (za(i3,i2)*zb(i2,i6) - za(i5,i3)*zb(i5,i6))*
     &    (za(i4,i2)*zb(i2,i6) - za(i5,i4)*zb(i5,i6))) + 
     & ((s(i1,i7)*za(i7,i2) + za(i6,i2)*za(i7,i1)*zb(i1,i6))**2*
     &    (za(i1,i2)*zb(i1,i6) + za(i7,i2)*zb(i7,i6))**2)/
     &  (s(i1,i7)*za(i3,i2)*za(i4,i2)*za(i5,i2)*za(i7,i1)*zb(i7,i6)*
     &    (za(i3,i1)*zb(i1,i6) - za(i7,i3)*zb(i7,i6))*
     &    (za(i4,i1)*zb(i1,i6) - za(i7,i4)*zb(i7,i6))*
     &    (za(i5,i1)*zb(i1,i6) - za(i7,i5)*zb(i7,i6))) - 
     & ((t(i1,i5,i7)*za(i7,i2) + za(i6,i2)*za(i7,i1)*zb(i1,i6) + 
     &       za(i6,i2)*za(i7,i5)*zb(i5,i6))**2*
     &    (za(i1,i2)*zb(i1,i6) + za(i5,i2)*zb(i5,i6) + 
     &      za(i7,i2)*zb(i7,i6)))/
     &  (t(i1,i5,i7)*za(i3,i2)*za(i4,i2)*za(i5,i1)*
     &    (za(i3,i1)*zb(i1,i6) - za(i5,i3)*zb(i5,i6) - 
     &      za(i7,i3)*zb(i7,i6))*
     &    (za(i4,i1)*zb(i1,i6) - za(i5,i4)*zb(i5,i6) - 
     &      za(i7,i4)*zb(i7,i6))*
     &    (za(i5,i1)*zb(i1,i6) - za(i7,i5)*zb(i7,i6))) - 
     & ((t(i1,i3,i7)*za(i7,i2) + za(i6,i2)*za(i7,i1)*zb(i1,i6) + 
     &       za(i6,i2)*za(i7,i3)*zb(i3,i6))**2*
     &    (za(i1,i2)*zb(i1,i6) + za(i3,i2)*zb(i3,i6) + 
     &      za(i7,i2)*zb(i7,i6)))/
     &  (t(i1,i3,i7)*za(i3,i1)*za(i4,i2)*za(i5,i2)*
     &    (za(i3,i1)*zb(i1,i6) - za(i7,i3)*zb(i7,i6))*
     &    (za(i4,i1)*zb(i1,i6) + za(i4,i3)*zb(i3,i6) - 
     &      za(i7,i4)*zb(i7,i6))*
     &    (za(i5,i1)*zb(i1,i6) + za(i5,i3)*zb(i3,i6) - 
     &      za(i7,i5)*zb(i7,i6))) - 
     & ((t(i1,i4,i7)*za(i7,i2) + za(i6,i2)*za(i7,i1)*zb(i1,i6) + 
     &       za(i6,i2)*za(i7,i4)*zb(i4,i6))**2*
     &    (za(i1,i2)*zb(i1,i6) + za(i4,i2)*zb(i4,i6) + 
     &      za(i7,i2)*zb(i7,i6)))/
     &  (t(i1,i4,i7)*za(i3,i2)*za(i4,i1)*za(i5,i2)*
     &    (za(i3,i1)*zb(i1,i6) - za(i4,i3)*zb(i4,i6) - 
     &      za(i7,i3)*zb(i7,i6))*
     &    (za(i4,i1)*zb(i1,i6) - za(i7,i4)*zb(i7,i6))*
     &    (za(i5,i1)*zb(i1,i6) + za(i5,i4)*zb(i4,i6) - 
     &      za(i7,i5)*zb(i7,i6)))
      return
      end


