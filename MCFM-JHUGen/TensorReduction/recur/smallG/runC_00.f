      subroutine runC_00(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat1,Shat2,S00,N0)
      implicit none
      include 'TRconstants.f' 
      include 'pvCnames.f' 
      include 'pvCv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      integer ep,N0,k,l,n,m,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex S00(-2:0),Shat1(np,-2:0),Shat2(np,np,-2:0),
     . bit,pole

      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat1(m,ep)-f(n)*f(m)*Cv(cc0+N0,ep))
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Cv(cc00+N0,ep-1)
      enddo
      Cv(cc00+N0,ep)=
     . -(+pole
     . +DetGr*Cv(cii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S00(ep)
     . +Gtwiddle(k,l)*Shat2(1,1,ep)-Gtwiddle(1,l)*Shat2(1,k,ep)
     . +Gtwiddle(k,l)*Shat2(2,2,ep)-Gtwiddle(2,l)*Shat2(2,k,ep)
     . +bit)/(6*Gtwiddle(k,l))

      enddo


      return
      end
  



