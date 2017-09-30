      subroutine genclust_hqrk(q,Rmin,qfinal,isub)
c--- performs the appropriate clustering with cones for the
c--- W+b+jet processes
c--- This is a much simpler version than the full routines, but
c--- the logic is much easier to follow
c--- To reject an event, return jets=-1
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'first.f'
      include 'nproc.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision pt,Rmin,R,Ry,Rbbmin,aetarap,ayrap
      double precision rtest,rab,raptest
      integer isub,i,j,k,nu,icandj,maxjet,jetindex(mxpart),ia,ib
      logical jetmerge,verbalg,pseudo,bjetmerge,is_hadronic
      common/Rbbmin/Rbbmin
      common/jetmerge/jetmerge
!$omp threadprivate(/jetmerge/)

      verbalg=.false.

      pseudo=.true.

      if (first) then
        if (pseudo) then
      write(6,*) 'Heavy quark clustering routine using pseudorapidity'
      else
      write(6,*) 'Heavy quark clustering routine using rapidity'
      endif
      write(6,*) '     (for both cuts and definition of Delta-R)'
      first=.false.
      endif

c--- transfer incoming and lepton momenta to the final array
      do j=1,4
        do nu=1,4
          qfinal(j,nu)=q(j,nu)
        enddo
      enddo

c--- simple clustering algorithm designed for WQj and ZQj
      jetmerge=.false.
      bjetmerge=.false.
      maxjet=0
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

c--- for the WQj(+Q) and ZQj(+Q) processes, the dipole subtraction
c--- terms contain just a b-quark and an untagged jet
      if ( (nproc .eq. 312) .or. (nproc .eq. 317)
     . .or.(nproc .eq. 322) .or. (nproc .eq. 327)
     . .or.(nproc .eq. 342) .or. (nproc .eq. 352) ) then
        if (isub .eq. 1)  jetlabel(2)='pp'
      endif

************************************************************************
c--- To check the ALPGEN total inclusive number
c      jets=-1
c      if ((pt(3,qjet) .lt. ptjetmin).or.(ayrap(3,qjet) .gt. etajetmax))
c     .  return
c      if ((Ry(qjet,1,3) .gt. Rmin) .and.
c     . (pt(1,qjet) .gt. ptjetmin).and.(ayrap(1,qjet) .lt. etajetmax))
c     . jets=2
c      if ((Ry(qjet,2,3) .gt. Rmin) .and.
c     . (pt(2,qjet) .gt. ptjetmin).and.(ayrap(2,qjet) .lt. etajetmax))
c     . jets=2
c      if (Ry(qjet,1,2) .gt. Rbbmin) return
c      do nu=1,4
c        qjet(4,nu)=qjet(1,nu)+qjet(2,nu)
c      enddo
c      if ((Ry(qjet,4,3) .gt. Rmin) .and.
c     . (pt(4,qjet) .gt. ptjetmin).and.(ayrap(4,qjet) .lt. etajetmax))
c     . jets=2
c      return
************************************************************************

c--- this is the number of candidate jets in the array qjet
      icandj=maxjet
c--- skip clustering if there are less than 2 jets
      if (icandj .lt. 2) then
        goto 77
c        write(6,*) 'Error in genclust_wbjt.f - icandj too small'
c        pause
      endif

      if (verbalg) write(6,*) 'Rmin=',Rmin
      if (verbalg) write(6,*) 'Rbbmin=',Rbbmin

c--- start another iteration
   66 continue

      if (verbalg) write(6,*) 'Starting iteration, # of cands =',icandj

c--- find the jet pair that meets the merging requirements and has
c--- the minimum separation in R
      rab=999d0
      ia=0
      ib=0
      jetmerge=.false.
      do j=1,icandj
      do k=j+1,icandj
        if (pseudo) then
          rtest=R(qjet,j,k)   ! Uses pseudorapidity
      else
          rtest=Ry(qjet,j,k)  ! Uses rapidity
      endif
        if (verbalg) write(6,*) jetlabel(j),',',jetlabel(k),' R =',rtest
c--- check to see if this pair should be merged - pair of b quarks
        if ( ((jetlabel(j).eq.'bq') .and. (jetlabel(k).eq.'ba'))
     .   .or.((jetlabel(j).eq.'ba') .and. (jetlabel(k).eq.'bq')) ) then
          if (rtest .lt. Rbbmin) then
            jetmerge=.true.
c--- ...... reset minimum separation if necessary
            if (rtest .lt. rab) then
              rab=rtest
              ia=j
              ib=k
            endif
          endif
c--- check to see if this pair should be merged - non b quarks
        else
          if (rtest .lt. Rmin) then
            jetmerge=.true.
          endif
c--- ...... reset minimum separation if necessary
          if (rtest .lt. rab) then
            rab=rtest
            ia=j
            ib=k
          endif
        endif
      enddo
      enddo

c--- if no jets need to be merged, we're done here
      if (jetmerge .eqv. .false.) goto 77

c--- check for an error
      if ((ia .eq. 0) .or. (ib .eq. 0)) then
        write(6,*) 'Error in genclust_wbjt: one of ia and ib is 0'
        stop
      endif

c--- if two b's will be merged, set the flag
      if ( ((jetlabel(ia).eq.'bq') .and. (jetlabel(ib).eq.'ba'))
     . .or.((jetlabel(ia).eq.'ba'). and. (jetlabel(ib).eq.'bq')) ) then
            bjetmerge=.true.
      endif
c--- do merging
c--- ...... set label for merged jet
      if     ((jetlabel(ia) .eq. 'bq').or.(jetlabel(ib) .eq. 'bq')) then
        jetlabel(ia)='bq'
      elseif ((jetlabel(ia) .eq. 'ba').or.(jetlabel(ib) .eq. 'ba')) then
        jetlabel(ia)='ba'
      else
        jetlabel(ia)='pp'
      endif
c--- ...... merge momenta for merged jet
      do nu=1,4
        qjet(ia,nu)=qjet(ia,nu)+qjet(ib,nu)
      enddo
c--- ...... shuffle down the other jets
      do j=ib,icandj-1
        jetlabel(ib)=jetlabel(ib+1)
        do nu=1,4
          qjet(j,nu)=qjet(j+1,nu)
        enddo
      enddo
c--- ...... remove a jet from the end
      do nu=1,4
        qjet(icandj,nu)=0d0
      enddo
      icandj=icandj-1

c--- return to the next iteration if there is more than one jet left
      if (icandj .gt. 1) goto 66

c--- done with the clustering - just need to check pt and y now
   77 continue

      if (verbalg) write(6,*) 'After merging # of jets =',icandj

c--- set the "jetmerge" flag properly now
      if (icandj .lt. maxjet) then
        jetmerge=.true.
      else
        jetmerge=.false.
      endif

c--- now check jet pt and rapidity
      do i=1,icandj
        if (pseudo) then
        raptest=aetarap(i,qjet)
      else
        raptest=ayrap(i,qjet)
      endif
c--- .... for b-quarks
        if ((jetlabel(i) .eq. 'bq') .or. (jetlabel(i) .eq. 'ba')) then
          if ( (pt(i,qjet) .lt. ptbjetmin)
     .    .or. (raptest .gt. etabjetmax)) jetlabel(i)='ff'
c--- .... for generic jets
        else
          if ( (pt(i,qjet) .lt. ptjetmin)
     .    .or. (raptest .gt. etajetmax)) jetlabel(i)='ff'
        endif
        if (verbalg) write(6,*) jetlabel(i),pt(i,qjet),raptest
      enddo

c--- the final reckoning: transfer to qfinal
      jets=0
      do i=1,icandj
        if (jetlabel(i) .ne. 'ff') then
          jets=jets+1
          do nu=1,4
            qfinal(jetindex(i),nu)=qjet(i,nu)
          enddo
        endif
      enddo

c--- shuffle down jetlabel
      i=1
   88 continue
      if (jetlabel(i) .eq. 'ff') then
        do j=i,icandj-1
          jetlabel(j)=jetlabel(j+1)
        enddo
      endif
      if (i .lt. jets) then
        i=i+1
        goto 88
      endif

      if (verbalg) write(6,*) 'After pt and y # of jets =',jets
c      if (verbalg) pause

c--- now we have to do process-specific cuts

************************************************************************
c--- W(QQ) for inclusive Wb
************************************************************************
      if ((nproc .eq. 301) .or. (nproc .eq. 306)) then
c--- In this category, there should be just one jet, formed by merging
        if ((jets .ne. 1) .or. (jetmerge .eqv. .false.)) then
          jets=-1
        endif
        return
      endif

************************************************************************
c--- WQ(+Q) for inclusive Wb
************************************************************************
      if ((nproc .eq. 302) .or. (nproc .eq. 307)) then
c--- In this category, there should be just one jet, formed by merging
        if ((jets .ne. 1) .or. (jetmerge)) then
          jets=-1
        endif
        return
      endif

************************************************************************
c--- WQQ/ZQQ + extra jet in real
************************************************************************
      if ((nproc .eq. 21) .or. (nproc .eq. 26)
     ..or.(nproc .eq. 51) .or. (nproc .eq. 56)) then
c--- In this category, there should be at least two jets, both b-quarks
        if (jets .lt. 2) then
          jets=-1
        endif
c--- only need to check jet identities if there are just 2 jets
        if (jets .eq. 2) then
        if (((jetlabel(1) .eq. 'bq').and.(jetlabel(2) .eq. 'ba'))
     . .or. ((jetlabel(1) .eq. 'ba').and.(jetlabel(2) .eq. 'bq'))) then
          continue
        else
          jets=-1
        endif
        endif
        return
      endif

************************************************************************
c--- WQQj/ZQjj/ZQQj
************************************************************************
      if ( (nproc .eq.  24) .or. (nproc .eq.  29)
     . .or.(nproc .eq. 346) .or. (nproc .eq. 347)
     . .or.(nproc .eq. 356) .or. (nproc .eq. 357)) then
c--- In this category, there should be three jets
        if (jets .lt. 3) then
          jets=-1
        endif
        return
      endif

************************************************************************
c--- WQj/ZQj - basic LO process + extra jet in real
************************************************************************
      if ( (nproc .eq. 311) .or. (nproc .eq. 316)
     . .or.(nproc .eq. 321) .or. (nproc .eq. 326)
     . .or.(nproc .eq. 341) .or. (nproc .eq. 351) ) then
c--- In this category, there should be two jets, one of which is a b
        if (jets .lt. 2) then
          jets=-1
        endif
c--- only need to check jet identities if there are just 2 jets
        if (jets .eq. 2) then
        if (((jetlabel(1) .eq. 'bq').and.(jetlabel(2) .eq. 'pp'))
     . .or. ((jetlabel(1) .eq. 'ba').and.(jetlabel(2) .eq. 'pp'))
     . .or. ((jetlabel(1) .eq. 'pp').and.(jetlabel(2) .eq. 'bq'))
     . .or. ((jetlabel(1) .eq. 'pp').and.(jetlabel(2) .eq. 'ba'))) then
          continue
        else
          jets=-1
        endif
        endif
        return
      endif

************************************************************************
c--- WQj/ZQj - basic LO process + extra Q in real
************************************************************************
      if ( (nproc .eq. 312) .or. (nproc .eq. 317)
     . .or.(nproc .eq. 322) .or. (nproc .eq. 327)
     . .or.(nproc .eq. 342) .or. (nproc .eq. 352) ) then
c--- In this category, there should be two jets, one of which is a b
c--- The cut on the Delta_R of the original b-bbar, to split the PS,
c--- is performed in "realint.f", so that neither the real contribution
c--- nor the counterterms are included in that case
c        if ((Ry(q,5,6) .lt. Rbbmin) .and. (isub .eq. 0)) then
c          jets=-1
c        endif
        if (jets .lt. 2) then
          jets=-1
        endif
c--- only need to check jet identities if there are just 2 jets
        if (jets .eq. 2) then
        if (((jetlabel(1) .eq. 'bq').and.(jetlabel(2) .eq. 'pp'))
     . .or. ((jetlabel(1) .eq. 'ba').and.(jetlabel(2) .eq. 'pp'))
     . .or. ((jetlabel(1) .eq. 'pp').and.(jetlabel(2) .eq. 'bq'))
     . .or. ((jetlabel(1) .eq. 'pp').and.(jetlabel(2) .eq. 'ba'))) then
          continue
        else
          jets=-1
        endif
        endif
        return
      endif

************************************************************************
c--- To check the ALPGEN approach to producing WQj
************************************************************************
c      if ((nproc .eq. 313) .or. (nproc .eq. 318)) then
c--- In this category, there should be two jets, one of which is a b
c--- jets should not have been merged
c        if ((jets .ne. 2) .or. (jetmerge .eqv. .true.)) then
c          jets=-1
c        endif
c        if (((jetlabel(1) .eq. 'bq').and.(jetlabel(2) .eq. 'ba'))
c     . .or. ((jetlabel(1) .eq. 'ba').and.(jetlabel(2) .eq. 'bq'))) then
c          jets=-1
c        endif
c        return
c      endif

************************************************************************
c--- W(QQ)j
************************************************************************
      if ((nproc .eq. 313) .or. (nproc .eq. 318)) then
c--- In this category, there should be two jets, one of which is a b
c--- also, jets should have been merged, not missed because of pt or y
        if ((jets .ne. 2) .or. (jetmerge .eqv. .false.)) then
          jets=-1
        endif
        if (((jetlabel(1) .eq. 'bq').and.(jetlabel(2) .eq. 'ba'))
     . .or. ((jetlabel(1) .eq. 'ba').and.(jetlabel(2) .eq. 'bq'))) then
          jets=-1
        endif
        return
      endif

************************************************************************
c--- WQjj
************************************************************************
      if ( (nproc .eq. 314) .or. (nproc .eq. 319)
     . .or.(nproc .eq. 324) .or. (nproc .eq. 329) ) then
c--- In this category, there should be three jets
        if (jets .lt. 3) then
          jets=-1
        endif
        return
      endif

      return
      end


      double precision function ry(p,i,j)
c---- calculate the jets separation between p(i) and p(j)
c--- Modification of "r.f" to use rapidity rather than pseudorapidity
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),r1,r2,dely,delphi,ei,ej
      integer i,j

      ei=p(i,4)
      ej=p(j,4)

      r1= (ei+p(i,3))*(ej-p(j,3))/
     .   ((ej+p(j,3))*(ei-p(i,3)))
      dely=0.5d0*dlog(r1)

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))
     .   /dsqrt((p(i,1)**2+p(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 .gt. +0.999999999D0) r2=+1D0
      if (r2 .lt. -0.999999999D0) r2=-1D0
      delphi=dacos(r2)

      ry=dsqrt(dely**2+delphi**2)

      return
      end

