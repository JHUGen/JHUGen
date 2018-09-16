      subroutine genclust_cone(q,Rmin,qfinal,isub)
      implicit none
      include 'types.f'
c---  clusters momenta using plabel to determine which 
c---  particles should be clustered. Forms 'jets' jets according to
c---  the Run II cone algorithm with cone size Rmin.
c---  Furthermore, the clustered jets are only observed if
c---  pT(jet) > ptjetmin and y(jet) < etajetmax
c--- 
c---  qfinal is the final vector q1,.... q(4+jets)
c---  where non-jet four vectors are set equal to the incoming q 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'kprocess.f'
      include 'plabel.f'
      real(dp):: q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      real(dp):: Rsep,Rmin,aetarap
      integer:: i,j,k,l,nu,iter,maxjet,ajet,jetindex(mxpart),isub
      character*2 finallabel(mxpart)
      real(dp):: protoq(20,4),deltarq,deltarj,et,etmax,net,
     & qshared(4),sharedet,getet
      integer:: maxproto,protoc(20,0:mxpart),eti,shared,
     & sharedc(20),ni
      logical:: jetmerge,failed,first,trackdoubleb,is_hadronic
c--- DEBUG
      parameter (Rsep=1._dp)  ! Default value
c      parameter (Rsep=1.3_dp)  ! Default value
c      parameter (Rsep=2.0_dp)  ! Usual (e.g. Snowmass) definition
      common/jetmerge/jetmerge
      data first/.true./
      save first
!$omp threadprivate(first)
!$omp threadprivate(/jetmerge/)
      
      if (first) then
       write(6,*)
       write(6,*) '*******  Cone algorithm additional parameter *******'
       write(6,*) '*                                                  *'
       write(6,79) '*    parton separation parameter, Rsep : ',Rsep
       write(6,*) '*                                                  *'
       write(6,*) '****************************************************'
       call flush(6)
       first=.false.
      endif
      
      jets=0
      maxjet=0
      jetmerge=.false.

c--- flag to determine whether or not to count "double b" jets
c--- (for Wbb, Wb+X processes only)
      if ( (kcase==kWbbmas) .or. (kcase==kW_bjet)) then
        trackdoubleb=.true.
      else
        trackdoubleb=.false.
      endif
      
      do i=1,mxpart
        do nu=1,4
        qfinal(i,nu)=0._dp
        enddo
      enddo


c--- pick out jets: note that we search to npart+2-isub, to get the
c--- number of particles right. Note that isub=0 for all calls except
c--- the dipole contributions, where isub=1.   
      do i=3,npart+2-isub
      if (is_hadronic(i)) then
        maxjet=maxjet+1
        jetindex(maxjet)=i
        jetlabel(maxjet)=plabel(i)
        do nu=1,4
          qjet(maxjet,nu)=q(i,nu)
        enddo
      endif
      enddo

c--- for no partons, just switch q into qfinal
      if (maxjet == 0) then
        do i=1,mxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        jets=0
        return
      endif

c--- skip clustering if we only have one parton  
      if (maxjet == 1) then
        jets=1
        do nu=1,4
          qfinal(1,nu)=qjet(1,nu)
        enddo
        finallabel(1)=jetlabel(1)
        goto 2
      endif
      
c--- set up the proto-jets
      maxproto=0
      do i=1,maxjet
        maxproto=maxproto+1
        protoc(maxproto,0)=1
        protoc(maxproto,1)=i
        do nu=1,4
          protoq(maxproto,nu)=qjet(i,nu)
        enddo
      enddo
      do i=1,maxjet
        do j=i+1,maxjet
          maxproto=maxproto+1
          protoc(maxproto,0)=2
          protoc(maxproto,1)=i
          protoc(maxproto,2)=j
          do nu=1,4
            protoq(maxproto,nu)=qjet(i,nu)+qjet(j,nu)
          enddo
          if (  (deltarq(maxproto,i,protoq) > Rmin)
     &     .or. (deltarq(maxproto,j,protoq) > Rmin)
     &     .or. (deltarq(i,j,protoq) > Rmin*Rsep) ) then
            maxproto=maxproto-1
          endif
        enddo
      enddo
      if (maxjet > 2) then
      do i=1,maxjet
        do j=i+1,maxjet
          do k=j+1,maxjet
            maxproto=maxproto+1
            protoc(maxproto,0)=3
            protoc(maxproto,1)=i
            protoc(maxproto,2)=j
            protoc(maxproto,3)=k
            do nu=1,4
              protoq(maxproto,nu)=qjet(i,nu)+qjet(j,nu)+qjet(k,nu)
            enddo
            if (  (deltarq(maxproto,i,protoq) > Rmin)
     &       .or. (deltarq(maxproto,j,protoq) > Rmin)
     &       .or. (deltarq(maxproto,k,protoq) > Rmin)
     &       .or. (deltarq(i,j,protoq) > Rmin*Rsep)
     &       .or. (deltarq(i,k,protoq) > Rmin*Rsep)
     &       .or. (deltarq(j,k,protoq) > Rmin*Rsep)) then
              maxproto=maxproto-1
            endif
          enddo
        enddo
      enddo
      endif
      if (maxjet > 3) then
      do i=1,maxjet
        do j=i+1,maxjet
          do k=j+1,maxjet
            do l=k+1,maxjet
              maxproto=maxproto+1
              protoc(maxproto,0)=4
              protoc(maxproto,1)=i
              protoc(maxproto,2)=j
              protoc(maxproto,3)=k
              protoc(maxproto,4)=l
              do nu=1,4
                protoq(maxproto,nu)=qjet(i,nu)+qjet(j,nu)
     &                             +qjet(k,nu)+qjet(l,nu)
              enddo
            if (  (deltarq(maxproto,i,protoq) > Rmin)
     &       .or. (deltarq(maxproto,j,protoq) > Rmin)
     &       .or. (deltarq(maxproto,k,protoq) > Rmin)
     &       .or. (deltarq(maxproto,l,protoq) > Rmin)
     &       .or. (deltarq(i,j,protoq) > Rmin*Rsep)
     &       .or. (deltarq(i,k,protoq) > Rmin*Rsep)
     &       .or. (deltarq(i,l,protoq) > Rmin*Rsep)
     &       .or. (deltarq(j,k,protoq) > Rmin*Rsep)
     &       .or. (deltarq(j,l,protoq) > Rmin*Rsep)
     &       .or. (deltarq(k,l,protoq) > Rmin*Rsep)) then
              maxproto=maxproto-1
            endif
            enddo
          enddo
        enddo
      enddo
      endif
      if (maxjet > 4) then
       write(6,*) 'Too many jets for this version of the cone algorithm'
       stop
      endif
                 
c      write(6,*) 'Found ',maxproto,' proto-jets'
      
      jets=0
      
      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

      if (maxproto == 0) goto 2

c--- find the highest Et proto-jet
      eti=0
      etmax=-1._dp
      do i=1,maxproto
        et=getet(protoq(i,4),protoq(i,1),protoq(i,2),protoq(i,3))
c        et=sqrt(protoq(i,1)**2+protoq(i,2)**2)
        if (et > etmax) then
          eti=i
          etmax=et
        endif
      enddo
      
c      write(6,*) 'Max Et proto-jet is ',eti
      
c--- check to see if any partons are shared by this proto-jet
      shared=0
      do i=1,maxproto
        sharedc(i)=0
        if (i .ne. eti) then
          do j=1,protoc(i,0)
            do k=1,protoc(eti,0)
              if (protoc(i,j) == protoc(eti,k)) then
                shared=shared+1
                sharedc(i)=1
              endif
            enddo
          enddo
        endif
      enddo
      
      if (shared == 0) then
c-- proto-jet does not share any partons - move it to qfinal and repeat
        jets=jets+1
        do nu=1,4
          qfinal(jets,nu)=protoq(eti,nu)
        enddo
        finallabel(jets)='pp'
        do i=1,protoc(eti,0)
          if (jetlabel(protoc(eti,i)) == 'bq') finallabel(jets)='bq'
          if (jetlabel(protoc(eti,i)) == 'ba') finallabel(jets)='ba'
        enddo
c--- special combination for W+heavy quarks
        if ((trackdoubleb) .and. (protoc(eti,0) == 2)) then
          if ( ((jetlabel(protoc(eti,1)) == 'bq') .and. 
     &          (jetlabel(protoc(eti,2)) == 'ba')) 
     &    .or. ((jetlabel(protoc(eti,1)) == 'ba') .and. 
     &          (jetlabel(protoc(eti,2)) == 'bq')) ) then
             finallabel(jets)='bb'
           endif 
      endif
        if ((trackdoubleb) .and. (protoc(eti,0) == 3)) then
             finallabel(jets)='bb'
        endif
c      if (protoc(eti,0) == 3) then
c        write(6,*) jetlabel(protoc(eti,1)),jetlabel(protoc(eti,2)),
c     &               jetlabel(protoc(eti,3)),' ->',finallabel(jets)
c      endif
c--- shuffle down the proto-jets
        do i=eti+1,maxproto
          do nu=1,4
            protoq(i-1,nu)=protoq(i,nu)
          enddo
          do j=0,mxpart
            protoc(i-1,j)=protoc(i,j)
          enddo
        enddo
        maxproto=maxproto-1
c        write(6,*) 'Found jet number ',jets
        goto 1
      endif

c--- a parton is shared: perform split/merge procedure
c      write(6,*) 'Need to do split/merge'

c--- calculate which proto-jet that shares has the highest Et      
      ni=0
      net=-1._dp
      do i=1,maxproto
        et=getet(protoq(i,4),protoq(i,1),protoq(i,2),protoq(i,3))
        if ((sharedc(i) == 1) .and. (et > net)) then
          ni=i
          net=et
        endif
      enddo
     
c--- calculate the shared Et
      do nu=1,4
        qshared(nu)=0._dp
      enddo
      do j=1,protoc(eti,0)
        do k=1,protoc(ni,0)
          if (protoc(eti,j) == protoc(ni,k)) then
            do nu=1,4
              qshared(nu)=qshared(nu)+qjet(protoc(eti,j),nu)
            enddo
          endif
        enddo
      enddo
      sharedet=getet(qshared(4),qshared(1),qshared(2),qshared(3))
      
c      write(6,*) 'Proto-jet is',eti
c      write(6,*) 'Highest et neighbour is',ni
c      write(6,*) 'Shared Et is',sharedet
c      write(6,*) 'Neighbour Et is',net
      
      if (sharedet/net > 0.5_dp) then
c---  we should merge the proto-jets
        do i=1,protoc(ni,0)
          shared=0
          do j=1,protoc(eti,0)
            if (protoc(ni,i) == protoc(eti,j)) shared=1
          enddo
c--- add cells that are not shared
          if (shared == 0) then
            protoc(eti,0)=protoc(eti,0)+1
            protoc(eti,protoc(eti,0))=protoc(ni,i)
            do nu=1,4
              protoq(eti,nu)=protoq(eti,nu)+qjet(protoc(ni,i),nu)
            enddo
          endif
        enddo
c--- shuffle down the proto-jets
        do i=ni+1,maxproto
          do nu=1,4
            protoq(i-1,nu)=protoq(i,nu)
          enddo
          do j=0,mxpart
            protoc(i-1,j)=protoc(i,j)
          enddo
        enddo
        maxproto=maxproto-1
c        write(6,*) 'Merged proto-jets',eti,' and ',ni
      else
c---  we should split the proto-jets
        do i=1,protoc(ni,0)
          shared=0
          do j=1,protoc(eti,0)
            if (protoc(ni,i) == protoc(eti,j)) shared=j
          enddo
c--- if a cell is shared, decide where to put it based on distance in Delta_R
c--- update the contents list and momentum of the protojet it's removed from
          if (shared > 0) then
            if (deltarj(protoc(ni,i),ni,qjet,protoq)
     &     < deltarj(protoc(ni,i),eti,qjet,protoq)) then
c--- shared cell is closer to neighbour, ni
              do j=shared+1,protoc(eti,0)
              protoc(eti,j-1)=protoc(eti,j)
              enddo
              protoc(eti,0)=protoc(eti,0)-1     
              do nu=1,4
                protoq(eti,nu)=protoq(eti,nu)-qjet(protoc(ni,i),nu)
              enddo
            else            
c--- shared cell is closer to original proto-jet, eti     
              do j=i+1,protoc(ni,0)
              protoc(ni,j-1)=protoc(ni,j)
              enddo
              protoc(ni,0)=protoc(ni,0)-1     
              do nu=1,4
                protoq(ni,nu)=protoq(ni,nu)-qjet(protoc(ni,i),nu)
              enddo
            endif
          endif
        enddo
c        write(6,*) 'Split proto-jets',eti,' and ',ni
      endif
      
c      pause
      goto 1                   
      
 2    continue 
 
c---- transfer qfinal --> qjet
      do i=1,jets
        jetlabel(i)=finallabel(i)
        do nu=1,4
          qjet(i,nu)=qfinal(i,nu)
        enddo
      enddo
                        
c      write(6,*) 'Finished finding jets: got ',jets
c      pause
      
c--- restore incoming partons
      do i=1,2
        do nu=1,4
          qfinal(i,nu)=q(i,nu)
        enddo
      enddo
c--- set all other momenta to zero and restore leptons
      do i=3,npart+2
        do nu=1,4
          qfinal(i,nu)=0._dp
          if (.not.(is_hadronic(i))) then
            qfinal(i,nu)=q(i,nu)
          endif
        enddo
      enddo
      
      
c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'

c--- flag whether or not any jets have been merged
      if (jets == maxjet) then
        jetmerge=.false.
      else
        jetmerge=.true.
      endif
      
c--- restore jets
      ajet=0
      do i=1,jets
c        write(*,*) 'Jet ',i,'(',jetlabel(i),')',jetindex(i)
c        write(*,*) 'pt: ',getet(qjet(i,4),qjet(i,1),
c     &               qjet(i,2),qjet(i,3)),' vs min. ',ptjetmin
c        write(*,*) 'ay: ',aetarap(i,qjet),' vs max. ',etajetmax
        if ((getet(qjet(i,4),qjet(i,1),qjet(i,2),qjet(i,3))
     &          > ptjetmin) .and.
     &      (aetarap(i,qjet) > etajetmin) .and.
     &      (aetarap(i,qjet) < etajetmax)) then  
        ajet=ajet+1
        do nu=1,4
          qfinal(jetindex(ajet),nu)=qjet(i,nu)
        enddo
        jetlabel(ajet)=jetlabel(i)
        endif
      enddo
      
c--- if no jets are removed by eta and pt cuts, then jets=ajet
      if (ajet < jets) then
        do i=ajet+1,jets
          do nu=1,4
            qfinal(jetindex(i),nu)=0._dp
          enddo
        enddo
        jets=ajet
      endif
      
c      write(*,*) '... and ',jets,' jets after pt and eta cuts'
c      do i=1,jets
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause

c-- check jets for heavy quark content and invariant mass
      call checkjets(jets,qfinal,isub,failed)
      if (failed) jets=-1

      return

   79 format(a42,f6.3,'    *')

      end
