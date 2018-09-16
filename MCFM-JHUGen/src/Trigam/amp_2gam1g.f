!====== C.Williams March 2013 
!====== Amplitude for q(i1)^-qb(i2)^+gluon(i3)^-gamma(i4)^+gamma(i5)^+ 
      
      function amp_2gam1g(p1,p2,p3,p4,p5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam1g
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      integer:: p1,p2,p3,p4,p5 

      amp_2gam1g=za(p2,p1)*za(p1,p3)**2/
     &     (za(p1,p5)*za(p1,p4)*za(p2,p5)*za(p2,p4))
      
      return 
      end
