c--- This subroutine is a basic jet clustering algorithim for 
c--- f(p1)+f(p2) --> gamma(p3)+f(p4)+f(p5) 


c---- Takes in pin, which has passed frixione cuts will cluster partons and determine whether 
c---- any partons lie in cone Rij < delta_0
      subroutine photgenclust(pin,Rmin,pfinal,isub,ipow) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'npart.f'
      include 'frag.f'
      real(dp):: pin(mxpart,4),Rmin,pjet(mxpart,4),
     & pfinal(mxpart,4)
      real(dp):: dijmin,dkmin,aetarap,pt,Rgen
      integer:: isub,i,nu,iter,ipow,nmin1,nmin2,maxjet,jetindex(mxpart) 
      integer:: nk,ajet
      logical:: jetmerge,insideacone,inclusive_inside_cone,
     & is_hadronic,is_photon
      integer:: photindex(npart),nphotons,j
      integer:: softjet
      common/jetmerge/jetmerge
!$omp threadprivate(/jetmerge/)

c--- this flag tells the algorithm whether or not to be inclusive of jets
c--- inside the photon isolation cone; this should be set to TRUE for an
c--- implementation of the original Frixione algorithm
c--- (c.f. Step 4 of hep-ph/9801442);
c--- other implementations (e.g. VBFNLO) correspond to setting this flag FALSE
      inclusive_inside_cone=.true.
      
      jets=0
      nphotons=0
      maxjet=0
      jetmerge =.false.


c---- Pick out jets 
      do i=3,npart+2-isub
         if (is_hadronic(i))then
            maxjet=maxjet+1
            jetindex(maxjet)=i
            jetlabel(maxjet)=plabel(i) 
            do nu=1,4
               pjet(maxjet,nu)=pin(i,nu)
            enddo
         elseif (is_photon(i)) then 
            nphotons=nphotons+1
            photindex(nphotons)=i
         endif
      enddo

      if (maxjet == 0 ) then 
         do i =1,mxpart 
            do nu=1,4
               pfinal(i,nu)=pin(i,nu)
            enddo 
         enddo
         jets=0
         return 
      endif

      if (maxjet == 1) goto 2

      iter = 0
           
 1    iter =iter +1 
      
      call findmind(pin,pjet,iter,maxjet,dijmin,nmin1,nmin2,ipow)
      
      call findminet(pin,pjet,iter,maxjet,dkmin,nk,ipow) 
      dkmin=dkmin*Rmin

      if (dijmin < dkmin) then 
         jetmerge = .true. 
         call combine(pjet,nmin1,nmin2) 
         
         call swapjet(pjet,jetindex,nmin2,maxjet)
         maxjet=maxjet-1
         iter =iter-1

      else
         jets=jets+1
         call swapjet(pjet,jetindex,jets,nk)
      endif

      if (iter < maxjet-1) goto 1


 2    continue 
      jets=jets+1

      pfinal(:,:)=0._dp
      do i=1,2
         do nu=1,4
            pfinal(i,nu)=pin(i,nu) 
         enddo 
      enddo

      do i=3,npart+2
         do nu=1,4
            pfinal(i,nu)=0._dp
             if (is_hadronic(i) .eqv. .false.)then
                pfinal(i,nu)=pin(i,nu)
             endif
          enddo
       enddo
       
      

c       if(isub == 0) then 
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'
c      endif
       
      
       ajet=0
       softjet=0

c--- loop over all partons
       do i=1,jets 

c--- check to see whether parton is inside one of the photon cones
         insideacone=.false.
         do j=1,nphotons           

c           write(6,*) i,j,Rgen(pjet,i,pin,photindex(j)),cone_ang
           if(Rgen(pjet,i,pin,photindex(j)) < cone_ang) then 
           insideacone=.true.
           endif
         enddo

         if (insideacone .and. inclusive_inside_cone) then
c--- if passed frix and is in isolation cone then do not apply cuts
c--- will add to jet tally
           softjet=softjet+1
         else
c--- if jet doesnt lie within photon cone apply cuts
           if ((pt(i,pjet) >= ptjetmin) .and. 
     &         (aetarap(i,pjet) >= etajetmin) .and.
     &         (aetarap(i,pjet) <= etajetmax)) then
             ajet=ajet+1
             do nu=1,4
               pfinal(jetindex(ajet),nu)=pjet(i,nu)
             enddo
           endif
         endif

       enddo

c       write(6,*) 'isub,jets,ajet,softjet',isub,jets,ajet,softjet
c       if (softjet == 1) pause
         
c--- if no jets are removed by eta and pt cuts, then jets=ajet
       if (ajet < jets) then
          do i=ajet+1,jets
             do nu=1,4
                pfinal(jetindex(i),nu)=0._dp
             enddo
          enddo        
       endif

       jets=ajet
       
c      write(6,*) 'input momenta' 
c      call writeout(pin)
c      write(6,*) '*****************************************'
c      write(6,*) ' output momenta' 
c      call writeout(pfinal) 
c      pause
      
      return 
      end
