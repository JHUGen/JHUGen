      subroutine run_00i(k,l,i1,DetGr,f,Gtwiddle,Gtt,
     . Shat2,Shat3,Shat3zz,S00i,N0)
      implicit none
C---Fixes D00i according to 5.44 Denner-Dittmaier
C---knowing Di,D00 with correction of order Delta*Diii
      include 'TRconstants.f'
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,k,l,n,m,i1
      double precision DetGr,Gtwiddle(3,3),Gtt(3,3,3,3),
     . f(3)
      double complex S00i(3,-2:0),Shat3zz(3,-2:0),
     . Shat2(3,3,-2:0),Shat3(3,z2max,-2:0),pole,bit
 
      do ep=-2,0
      bit=czip
      do n=1,3
      do m=1,3
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat2(m,i1,ep)
     . +2d0*delta(n,i1)*Shat3zz(m,ep)-f(n)*f(m)*Dv(di(i1)+N0,ep)
     . -2d0*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Dv(dd00+N0,ep))
      enddo
      enddo
      pole=0d0
      if (ep.ge.-1) pole=-4d0*Gtwiddle(k,l)*Dv(dzzi(i1)+N0,ep-1)

      Dv(dzzi(i1)+N0,ep)=-(
     . +pole
     . +DetGr*Dv(diii(z3(k,l,i1))+N0,ep)
     . -Gtwiddle(k,l)*S00i(i1,ep)
     . -Gtwiddle(1,l)*Shat3(1,z2(k,i1),ep)
     . -Gtwiddle(2,l)*Shat3(2,z2(k,i1),ep)
     . -Gtwiddle(3,l)*Shat3(3,z2(k,i1),ep)
     . +Gtwiddle(k,l)*(Shat3(1,z2(1,i1),ep)
     .                +Shat3(2,z2(2,i1),ep)
     .                +Shat3(3,z2(3,i1),ep))
     . +bit)/(8d0*Gtwiddle(k,l))
      enddo

      return
      end
  



