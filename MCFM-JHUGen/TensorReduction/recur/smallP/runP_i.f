      subroutine runP_i(k,i1,f,Gr,Shat2,N0)
      implicit none
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,k,i1,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat2(np,np,-2:0)
       
      do ep=-2,0
      Dv(di(i1)+N0,ep)=
     . (Shat2(k,i1,ep)
     . -2d0*delta(k,i1)*Dv(dd00+N0,ep)
     . -Gr(k,1)*Dv(dii(z2(1,i1))+N0,ep) 
     . -Gr(k,2)*Dv(dii(z2(2,i1))+N0,ep)
     . -Gr(k,3)*Dv(dii(z2(3,i1))+N0,ep))/f(k) 
      enddo
      
      return
      end
