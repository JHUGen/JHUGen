C--- Implementation of arXiv:0906.0008v1, Eqs. (A.5), (A.6) and (A.7)
c--- with factor of i removed

      double complex function A0phiAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.5) second line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phiAQggmppm=-za(j1,j4)**2*za(j2,j4)
     .                 /(za(j1,j2)*za(j2,j3)*za(j3,j4))
      
      return
      end
      
      double complex function A0phiAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.5) third line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phiAQggmpmp=za(j1,j3)**3/(za(j1,j2)*za(j3,j4)*za(j4,j1))
      
      return
      end
      
      double complex function A0phidAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.5) fifth line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phidAQggmppm=-zb(j2,j3)**2*zb(j1,j3)
     .                   /(zb(j1,j2)*zb(j3,j4)*zb(j4,j1))
      
      return
      end
      
      double complex function A0phidAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.5) sixth line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phidAQggmpmp=zb(j2,j4)**3
     .                  /(zb(j1,j2)*zb(j2,j3)*zb(j3,j4))
      
      return
      end
      
      double complex function A0phiAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.6)
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double complex zab2
      double precision s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      A0phiAQggmpmm=-zab2(j3,j1,j4,j2)**2*za(j4,j1)
     .               /zb(j2,j4)/s3(j4,j1,j2)*(1d0/s(j1,j2)+1d0/s(j4,j1))
     .              -zab2(j4,j1,j3,j2)**2*za(j1,j3)
     .               /(zb(j2,j3)*s(j1,j2)*s3(j1,j2,j3))
     .              +zab2(j1,j3,j4,j2)**2
     .               /(za(j1,j2)*zb(j2,j4)*zb(j2,j3)*zb(j3,j4))            
     
      
      return
      end
      
      double complex function A0phidAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.7)
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double complex zab2
      double precision s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      A0phidAQggmppp=-zab2(j1,j2,j3,j4)**2*zb(j2,j3)
     .               /za(j1,j3)/s3(j1,j2,j3)*(1d0/s(j1,j2)+1d0/s(j2,j3))
     .              +zab2(j1,j2,j4,j3)**2*zb(j2,j4)
     .               /(za(j1,j4)*s(j1,j2)*s3(j4,j1,j2))
     .              -zab2(j1,j3,j4,j2)**2
     .               /(zb(j1,j2)*za(j1,j3)*za(j1,j4)*za(j3,j4))    
     
      
      return
      end
      
