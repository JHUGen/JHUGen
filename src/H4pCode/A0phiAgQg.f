C--- Implementation of arXiv:0906.0008v1, Eqs. (A.8), (A.9) and (A.10)
c--- with factor of i removed

      double complex function A0phiAgQgmppm(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.8) first line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phiAgQgmppm=-za(j1,j4)**2/(za(j1,j2)*za(j2,j3))
      
      return
      end
      
      double complex function A0phiAgQgmmpp(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.8) second line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phiAgQgmmpp=-za(j1,j2)**2/(za(j3,j4)*za(j4,j1))
      
      return
      end
      
      double complex function A0phidAgQgmppm(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.8) third line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phidAgQgmppm=zb(j2,j3)**2/(zb(j3,j4)*zb(j4,j1))
      
      return
      end
      
      double complex function A0phidAgQgmmpp(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.8) fourth line
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      
      A0phidAgQgmmpp=zb(j3,j4)**2/(zb(j1,j2)*zb(j2,j3))
      
      return
      end
      
      double complex function A0phiAgQgmmpm(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.9)
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double complex zab2
      double precision s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      A0phiAgQgmmpm=-zab2(j4,j1,j2,j3)**2
     .               /(zb(j1,j2)*zb(j2,j3)*s3(j1,j2,j3))
     .              -zab2(j2,j1,j4,j3)**2
     .               /(zb(j3,j4)*zb(j4,j1)*s3(j3,j4,j1))
      
      return
      end
      
      double complex function A0phidAgQgmppp(j1,j2,j3,j4,za,zb)
      implicit none
c--- (A.9)
c----with factor of i removed
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4
      double complex zab2
      double precision s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      A0phidAgQgmppp=+zab2(j1,j2,j3,j4)**2
     .                /(za(j1,j2)*za(j2,j3)*s3(j1,j2,j3))
     .               +zab2(j1,j3,j4,j2)**2
     .                /(za(j3,j4)*za(j4,j1)*s3(j3,j4,j1))
      
      return
      end
      
