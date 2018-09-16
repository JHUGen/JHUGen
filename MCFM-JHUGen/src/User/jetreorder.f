      subroutine jetreorder(p,q,isub)
      implicit none
      include 'types.f'
c--- Routine to inspect the entries of 'jetlabel' and compare with 'plabel'
c--- in order to put clustered jets back in the correct positions.
c--- This is important if one wants to identify particles in the plotting
c--- routine with those in the process.DAT file.
c--- At present this routine is aimed at top processes with top decay
c--- (i.e. putting b and anti-b jets in correct positions)
c---
c--- Input array of momenta is p, output re-ordered array is q
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'npart.f'
      integer:: i,j,k,nu,isub,maxjet,ijet(mxpart),jetindex(mxpart)
      real(dp):: p(mxpart,4),q(mxpart,4)
      logical:: alreadyfound,is_hadronic
      
c      call writeout(p)
      
      if (jets == 0) then
        do i=1,mxpart
      do nu=1,4
        q(i,nu)=p(i,nu)
      enddo
      enddo
      return
      endif
      
      maxjet=0
      do i=3,npart+2-isub
        ijet(i)=0
        if (is_hadronic(i)) then
          maxjet=maxjet+1
          jetindex(maxjet)=i
        ijet(i)=-1 ! If no match is found, this signifies a parton is lost
        if (jets > 0) then
          do j=1,jets
            if (plabel(i) == jetlabel(j)) then
              alreadyfound=.false.
              if (maxjet > 0) then
              do k=1,maxjet      
                if (ijet(jetindex(k)) == j) then
                  alreadyfound=.true.
                endif
              enddo
            endif
              if (.not.(alreadyfound)) ijet(i)=j
            endif
          enddo
        endif
      endif
      enddo
      
c      call writeout(p)
      
c      do i=3,npart+2-isub
c      write(6,*) i,ijet(i)
c      enddo
c      write(6,*) 'maxjet=',maxjet
      
c      do j=1,jets
c      write(6,*) j,jetlabel(j)
c      enddo
      
      do i=3,npart+2-isub
        if     (ijet(i) == 0) then
          j=i ! Not a jet, no swap required
      elseif (ijet(i) == -1) then
        j=-1 ! Jet lost, need to set components to zero
        else
          j=jetindex(ijet(i))
        endif
      if (j > 0) then
        do nu=1,4
          q(i,nu)=p(j,nu)
          enddo
      else
        do nu=1,4
          q(i,nu)=zip
          enddo
      endif
      enddo
      
      do i=1,2
        do nu=1,4
        q(1,nu)=p(1,nu)
        q(2,nu)=p(2,nu)
      enddo
      enddo
      
      do i=npart+3-isub,mxpart
      do nu=1,4
         q(i,nu)=zip
      enddo
      enddo
      
c      call writeout(q)
c      pause
 
c--- reset jetlabel to correspond to plabel      
      do i=1,jets
      jetlabel(i)=plabel(jetindex(i))
c      write(6,*) 'i,jetindex(i),jetlabel(i)',i,jetindex(i),jetlabel(i)
      enddo
c      write(6,*)
      
c      call writeout(q)
c      write(6,*) 'jets =',jets
c      write(6,*) 'npart,isub',npart,isub
c      pause
      
      return
      end
      
