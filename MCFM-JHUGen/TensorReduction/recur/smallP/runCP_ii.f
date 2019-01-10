      subroutine runCP_ii(k,i1,i2,f,Gr,Shat3,N0)
      implicit none
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,k,i1,i2,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat3(np,z2max,-2:0)
       
      do ep=-2,0
      Cv(cii(z2(i1,i2))+N0,ep)=
     . (Shat3(k,z2(i1,i2),ep)
     . -2d0*delta(k,i1)*Cv(czzi(i2)+N0,ep)
     . -2d0*delta(k,i2)*Cv(czzi(i1)+N0,ep)
     . -Gr(k,1)*Cv(ciii(z3(1,i1,i2))+N0,ep) 
     . -Gr(k,2)*Cv(ciii(z3(2,i1,i2))+N0,ep))/f(k) 
      enddo
      
      return
      end
