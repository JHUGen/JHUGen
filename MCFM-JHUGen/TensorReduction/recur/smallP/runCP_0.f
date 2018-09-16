      subroutine runCP_0(k,f,Gr,Shat1,N0)
      implicit none
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,k,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat1(np,-2:0)
       
      do ep=-2,0
      Cv(cc0+N0,ep)=
     . (Shat1(k,ep)
     . -Gr(k,1)*Cv(ci(1)+N0,ep) 
     . -Gr(k,2)*Cv(ci(2)+N0,ep))/f(k) 
      enddo
      
      return
      end
