      subroutine genclust_kt_ab(q,Rmin,qfinal,isub,ipow)
      implicit none
      include 'types.f'
c--- Clusters momenta using plabel to determine which 
c--- particles should be clustered. Forms 'jets' jets according to
c--- the standard kT algorithm with cone size Rmin.
c--- Furthermore, the clustered jets are only observed if
c--- pT(jet) >= ptjetmin and etajetmin <= |eta(jet)| <= etajetmax
c--- 
c--- qfinal is the final vector q1,.... q(4+jets)
c--- where non-jet four vectors are set equal to the incoming q 
c---
c--- Modified 23/11/09: generalized to include whole class of kt-like
c---                    algorithms with measures raised to the power
c---                    "ipow" passed into routine. In particular,
c---                    ipow = +1 (normal kt), ipow = -1 ("anti-kt")
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'is_functions_com.f'

      real(dp):: q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      real(dp):: pt,Rmin,dijmin,dkmin,ayrap
      integer:: i,nu,iter,nmin1,nmin2,maxjet,nk,
     & ajet,jetindex(mxpart),isub,ipow
      logical:: jetmerge,failed,is_spechadronic
      common/jetmerge/jetmerge
!$omp threadprivate(/jetmerge/)

      jets=0
      maxjet=0
      jetmerge=.false.
      
      do i=1,mxpart
        do nu=1,4
        qfinal(i,nu)=0._dp
        enddo
      enddo

c--- pick out jets: note that we search to npart+2-isub, to get the
c--- number of particles right. Note that isub=0 for all calls except
c--- the dipole contributions, where isub=1.   
      do i=3,npart+2-isub
      if (is_spechadronic(i)) then
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
      if (maxjet == 1) goto 2

c--- for W+bbj, skip if b and b-bar are too close together
c      if ( ((nproc==292) .or. (nproc==297))
c     &     .and. (isub ==0) .and. (R(q,5,6) < Rmin) ) then
c        jets=-1
c        return
c      endif

      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

c      write(*,*) 'iter ',iter
c      write(*,*) 'jets ',jets
c      write(*,*) 'maxjet ',maxjet

c--- step1: find (i,j) pair with lowest measure of all non-jets so far
      call findmind(q,qjet,iter,maxjet,dijmin,nmin1,nmin2,ipow)
 
c--- step2: find jet K with lowest Et
      call findminet(q,qjet,iter,maxjet,dkmin,nk,ipow)
      dkmin=dkmin*Rmin

c      write(*,*) 'Comparing pair (',nmin1,',',nmin2,') value of'
c      write(*,*) 'dijmin = ',dijmin,' with ',nk,' value of dk = ',dkmin
      
c--- step3: compare the two ...      
      if (dijmin < dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2
        jetmerge=.true.
        call combine_ab(qjet,nmin1,nmin2)
c--- combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet(qjet,jetindex,nmin2,maxjet)        
        maxjet=maxjet-1
        iter=iter-1
c        do i=1,maxjet
c          do j=1,4
c            write(*,*) 'qjet(',i,',',nu,') = ',qjet(i,nu)
c          enddo
c        enddo
      else
c---  ... we've finished a jet
        jets=jets+1
c        write(*,*) 'Now swapping ',jets,' and ',nk
        call swapjet(qjet,jetindex,jets,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter < maxjet-1) goto 1
      
 2    continue      
      jets=jets+1

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
          if (.not.(is_spechadronic(i))) then
            qfinal(i,nu)=q(i,nu)
          endif
        enddo
      enddo
      
c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'
     
c--- restore jets

      ajet=0

      do i=1,jets

c        write(*,*) 'Jet ',i,'(',jetlabel(i),')',jetindex(i)
c        write(*,*) 'pt: ',pt(i,qjet),' vs min. ',ptjetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs min. ',etajetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs max. ',etajetmax

        if ((pt(i,qjet) >= ptjetmin) .and.
     &      (ayrap(i,qjet) >= etajetmin) .and.
     &      (ayrap(i,qjet) <= etajetmax)) then 
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
      end
