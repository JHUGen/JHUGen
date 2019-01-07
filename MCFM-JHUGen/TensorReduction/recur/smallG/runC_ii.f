      subroutine runC_ii(j,i1,i2,DetGr,Xtwiddle0,Gtwiddle,Shat3,N0)
      implicit none
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,j,i1,i2,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),Shat3s(np,z2max,-2:0)
       
      do ep=-2,0
      do n=1,np
      Shat3s(n,z2(i1,i2),ep)=Shat3(n,z2(i1,i2),ep)
     .   -2d0*(delta(n,i1)*Cv(N0+czzi(i2),ep)
     .        +delta(n,i2)*Cv(N0+czzi(i1),ep))
c      if (ep.eq.-1) then
c      write(6,'(a8,3i3,2(f20.15))') 'Shat3s ',n,i1,i2,
c     .    Shat3(n,z2(i1,i2),ep)
c      write(6,'(a8,3i3,2(f20.15))') 'Shat3s ',n,i1,i2,
c     .   -2d0*delta(n,i1)*Cv(N0+czzi(i2),ep)
c      write(6,'(a8,3i3,2(f20.15))') 'Shat3s ',n,i1,i2,
c     .   -2d0*delta(n,i2)*Cv(N0+czzi(i1),ep)
c      endif
      enddo 
      Cv(cii(z2(i1,i2))+N0,ep)=-(   
     . +Gtwiddle(j,1)*Shat3s(1,z2(i1,i2),ep)
     . +Gtwiddle(j,2)*Shat3s(2,z2(i1,i2),ep)
     . -DetGr*Cv(ciii(z3(j,i1,i2))+N0,ep))/Xtwiddle0(j)
     
c      if (ep .eq. -1) then
c      write(6,*) 'Cii ',j,i1,i2,
c     . +Gtwiddle(j,1)*Shat3s(1,z2(i1,i2),ep)/Xtwiddle0(j)
c      write(6,*) 'Cii ',j,i1,i2,
c     . +Gtwiddle(j,2)*Shat3s(2,z2(i1,i2),ep)/Xtwiddle0(j)
c      write(6,*) 'Cii ',j,i1,i2,
c     . -DetGr*Cv(ciii(z3(j,i1,i2))+N0,ep)/Xtwiddle0(j)
c      write(6,*)
c      endif
      
      enddo

c      pause

      return
      end
