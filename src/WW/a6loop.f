      double complex function a6loops(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
c---  DKS Eq. 3.15
      integer j1,j2,j3,j4,j5,j6
      double complex a6loopa
      a6loops=a6loopa(j1,j2,j3,j4,j5,j6,za,zb)
     .       +a6loopa(j1,j2,j6,j5,j4,j3,za,zb)
      return
      end
        
      double complex function a6loopa(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c---  DKS Eq. 2.10 for alpha = a
c---  note that (-i) included in A(alpha),A(tree,alpha)
c---  so no factor of (+i) in front of F(alpha)
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      double complex tree,Vpole,a6treea,fa
      tree=a6treea(j1,j2,j3,j4,j5,j6,za,zb)
      a6loopa=tree*Vpole(s(1,2))+fa(j1,j2,j3,j4,j5,j6,za,zb)
      return 
      end

c      double complex function a6loopb(j1,j2,j3,j4,j5,j6,za,zb)
c      implicit none
c---  DKS Eq. 2.10 for alpha = b
c---  note that (-i) included in A(alpha),A(tree,alpha)
c---  so no factor of (+i) in front of F(alpha)
c      include 'constants.f'
c      include 'sprods_com.f'
c      include 'zprods_decl.f'
c      integer j1,j2,j3,j4,j5,j6
c      double complex tree,Vpole,a6treeb
        
c      tree=a6treeb(j1,j2,j3,j4,j5,j6,za,zb)    
c      a6loopb=tree*Vpole(s(j1,j2))
c
c      return 
c      end

              

        
