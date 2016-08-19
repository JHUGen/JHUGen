      subroutine pvAfill(m1sq,N)
C    N is the offset in the storage
      implicit none
      include 'pvAnames.f'
      include 'TRconstants.f'
      include 'TRonshellcutoff.f'
      include 'pvAv.f'
      include 'TRscale.f'
      integer N,Np,ep,j
      double precision m1sq
      double complex trI1
      logical,save::first=.true.
      logical,save::scaleset=.false.
      double precision,save::id(0:2),idp2(0:2)
!$omp threadprivate(scaleset,first,id,idp2)

      if (first) then
      first=.false.
C--id=1/D
      id(0)=0.25d0
      id(1)=id(0)*0.5d0
      id(2)=id(1)*0.5d0
C--idp2=1/[D+2]
      idp2(0)=1d0/6d0
      idp2(1)=idp2(0)/3d0
      idp2(2)=idp2(1)/3d0
      endif

      if (scaleset .neqv. .true.) then
      scaleset=.true.
      if ((scale .eq. -1d12) .and. (musq .eq. -1d12)) then
      write(6,*) 'Did you forget to call setmudim?'
      write(6,*) 'Setting scale to scale=1d0'
      scale=1d0
      musq=1d0
      endif
      endif

      if (abs(m1sq/musq) .lt. onshellcutoff) then
c      write(6,*) 'setting zero mass, tadpole to zero'
c      write(6,*) 'm1sq=',m1sq
      do Np=N+1,N+Naa
      do ep=-2,0
      Av(Np,ep)=czip
      enddo 
      enddo 
      
      return

      else

      do ep=-2,0
      Av(aa0+N,ep)=trI1(m1sq,musq,ep)
      enddo
      
C Id,A00(m1?)=m1^2*A0(m1)/D;
      do ep=-2,0
        Av(aa00+N,ep)=czip
          do j=0,ep+2
          Av(aa00+N,ep)=Av(aa00+N,ep)+m1sq*Av(aa0+N,ep-j)*id(j)
          enddo
      enddo 
C Id,A0000(m1?)=m1^2*A00(m1)/[D+2];
      do ep=-2,0
      Av(aa0000+N,ep)=czip
          do j=0,ep+2
          Av(aa0000+N,ep)=Av(aa0000+N,ep)+m1sq*Av(aa00+N,ep-j)*idp2(j)
          enddo
      enddo
      endif

      return 
      end
