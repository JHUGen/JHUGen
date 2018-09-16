      function ggZgam_vec(st,p,j1,j2,j3,j4,j5,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: ggZgam_vec
      
c-----Function for the vector parts of gg->gamma(p3)+Z(p(4)+p(5))
c-----st = 'abc' where a,b label the gluon helicity and c labels the
c-----photon helicity; Z decay is always left-handed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      
      real(dp):: p(mxpart,4)
      integer:: j1,j2,j3,j4,j5
      character*3 st
      complex(dp):: L0,l1,Lsm1,ppp,ppm,pmm
      real(dp):: t
C-----statement functions     
      ppp(j1,j2)=2._dp*(
     & za(j1,j4)**2*zb(j3,j1)/(za(j1,j2)**2*za(j1,j3)*za(j4,j5)) 
     & -0.5_dp*zb(j5,j3)**2/(za(j1,j2)**2*zb(j5,j4)))
     
      ppm(j1,j2)=+2._dp*(
     &  2._dp*L0(-t(j1,j3,j2),-s(j1,j3))
     & *za(j2,j3)*za(j1,j4)*za(j2,j4)*zb(j2,j1)
     & /(zb(j3,j1)*za(j1,j2)**3*za(j4,j5))

     & -L1(-t(j1,j2,j3),-s(j1,j3))
     & *za(j2,j4)**2*zb(j2,j1)**2
     & /(zb(j3,j1)**2*za(j1,j2)**2*za(j4,j5)) 

     & -za(j1,j3)*za(j2,j4)*zb(j2,j1)*zb(j5,j1)
     & /(zb(j3,j1)*za(j1,j2)**2*za(j4,j5)*zb(j5,j4)))

      pmm(j2,j3)=-2._dp*(
     &  2._dp*L0(-t(j1,j2,j3),-s(j1,j3))
     & *zb(j2,j1)*zb(j3,j5)*zb(j2,j5)*za(j2,j3)
     & /(za(j1,j3)*zb(j3,j2)**3*zb(j5,j4))

     & -L1(-t(j1,j2,j3),-s(j1,j3))
     & *zb(j2,j5)**2*za(j2,j3)**2
     & /(za(j1,j3)**2*zb(j3,j2)**2*zb(j5,j4))  
       
     & -zb(j3,j1)*zb(j2,j5)*za(j2,j3)*za(j4,j3)
     & /(za(j1,j3)*zb(j3,j2)**2*zb(j5,j4)*za(j4,j5)))

      ggZgam_vec=(0._dp,0._dp)

      if(st=='+++') then 
      ggZgam_vec=ppp(j1,j2)+ppp(j2,j1) 
      
      elseif(st=='++-') then 
      ggZgam_vec=
     & +2._dp*Lsm1(-s(j2,j3),-t(j1,j2,j3),-s(j1,j3),-t(j1,j2,j3))
     & *(za(j1,j3)**2*za(j2,j4)**2+za(j1,j4)**2*za(j2,j3)**2)
     & /(za(j1,j2)**4*za(j4,j5))
     & +ppm(j1,j2)+ppm(j2,j1) 

      elseif(st=='+--') then 
         
      ggZgam_vec=
     & -2._dp*Lsm1(-s(j1,j3),-t(j1,j2,j3),-s(j1,j2),-t(j1,j2,j3))
     & *(zb(j3,j1)**2*zb(j2,j5)**2+zb(j3,j5)**2*zb(j1,j2)**2)
     & /(zb(j3,j2)**4*zb(j5,j4))
     & +pmm(j2,j3)+pmm(j3,j2) 

      else
         write(6,*) 'Unknown helicity '
         stop
      endif
      return 
      end
    

      
