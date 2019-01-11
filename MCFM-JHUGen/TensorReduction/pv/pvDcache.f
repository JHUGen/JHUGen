      integer function pvDcache(p1s,p2s,p3s,p4s,p1p2,p2p3,
     . m1s,m2s,m3s,m4s)
      implicit none
      include 'pvDnames.f'
      include 'TRclear.f'
      include 'TRonshellcutoff.f'
      include 'pvforcerecalc.f'
      include 'pvDitry.f'
      double precision para(Pdd),p1s,p2s,p3s,p4s,p1p2,p2p3,
     . m1s,m2s,m3s,m4s
      integer j,jtable,Ntrue
      double precision,save:: tableD(Pdd,Ndmax)
      integer,save:: Nstore=0
!$omp threadprivate(tableD,Nstore)

      if (clear(4)) then
      clear(4)=.false.
      Nstore=0
      endif

      if (Nstore .gt. Ndmax) then
      print * 
      print *, 'pvDcache:Nstore .gt. Ndmax'
      print *, 'pvDcache:Nstore,Ndmax',Nstore,Ndmax
      print *, 'Either adjust Ndmax in Dnames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      para(1)=p1s
      para(2)=p2s
      para(3)=p3s
      para(4)=p4s
      para(5)=p1p2
      para(6)=p2p3
      para(7)=m1s
      para(8)=m2s
      para(9)=m3s
      para(10)=m4s
C if parameter set is found set pvDcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Pdd
        if (abs(para(j)-tableD(j,jtable)) .lt. 1d-8) Ntrue=Ntrue+1 
        enddo
      if (Ntrue .eq. Pdd) then
        pvDcache=(jtable-1)*Ndd
      if (pvforcerecalc) then
c--- although integral is cached, need to compute with recursion
          call pvDfill(p1s,p2s,p3s,p4s,p1p2,p2p3,m1s,m2s,m3s,m4s,
     &                 pvDcache)
        endif
        return
      endif
      enddo

C    if parameter set is not found we have to calculate
C    and fill the common block starting at position pvDcache
 20   pvDcache=Nstore*Ndd
      pvDitry(pvDcache)=-1 ! label tensor as unchecked
      Nstore=Nstore+1
      do j=1,Pdd
        if(abs(para(j)) .lt. onshellcutoff) para(j)=0d0
      enddo
      do j=1,Pdd
      tableD(j,Nstore)=para(j)
      enddo
c      call pvDfill(p1s,p2s,p3s,p4s,p1p2,p2p3,m1s,m2s,m3s,m4s,pvDcache)
      call pvDfill(para(1),para(2),para(3),para(4),para(5),
     &             para(6),para(7),para(8),para(9),para(10),pvDcache)
      return
      end
