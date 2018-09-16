      function filterW_bjet()
       implicit none
      include 'types.f'
      logical:: filterW_bjet
c--- this routine is specific to the "W_bjet" (and related)
c--- processes 411 and 416; it inspects the jets to check whether an
c--- event should be included
c--- routine returns FALSE if event does not pass the process-specific cuts

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'clustering.f'
      include 'jetlabel.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'notag.f'
      logical:: jetmerge,veto3jets
      common/jetmerge/jetmerge
!$omp threadprivate(/jetmerge/)

c--- setting this to true vetoes 3-jet events (c.f. CDF analysis)
      veto3jets=.true.

c--- IMPORTANT: generically, we don't want this cut. Here it can be applied
c---            for comparison with analysis of C. Neu et al. (CDF)
      if (veto3jets .and. (jets == 3)) then
        filterW_bjet=.false.
        return
      endif

c--- default: all cuts are passed and event is accepted
      filterW_bjet=.true.

c--- first check the number of jets (just as in includedipole)
      if ( ((jets .ne. nqcdjets-notag).and.(inclusive.eqv..false.))
     & .or.((jets < nqcdjets-notag).and.(inclusive.eqv..true.))) then
        filterW_bjet=.false.
      return
      endif

************************************************************************
c--- W(QQ) for inclusive Wb
************************************************************************
      if ((nproc == 301) .or. (nproc == 306)) then
c--- In this category, there should be just one jet, formed by merging
        if ((jets .ne. 1) .or. (jetmerge .eqv. .false.)) then
          filterW_bjet=.false.
        endif
        return
      endif

************************************************************************
c--- WQ(+Q) for inclusive Wb
************************************************************************
      if ((nproc == 302) .or. (nproc == 307)) then
c--- In this category, there should be just one jet, the other outside
c--- the acceptance (i.e. failing pt or rapidity requirements)
        if ((jets < 1) .or. (jetmerge)) then
          filterW_bjet=.false.
        endif
        return
      endif

************************************************************************
c--- WQQ/ZQQ + extra jet in real
************************************************************************
      if ((nproc == 21) .or. (nproc == 26)
     ..or.(nproc == 51) .or. (nproc == 56)) then
c--- In this category, there should be at least two jets, both b-quarks
        if (jets < 2) then
          filterW_bjet=.false.
        endif
c--- only need to check jet identities if there are just 2 jets
        if (jets == 2) then
        if (((jetlabel(1) == 'bq').and.(jetlabel(2) == 'ba'))
     & .or. ((jetlabel(1) == 'ba').and.(jetlabel(2) == 'bq'))) then
          continue
        else
          filterW_bjet=.false.
        endif
        endif
        return
      endif

************************************************************************
c--- WQQj/ZQjj/ZQQj
************************************************************************
      if ( (nproc ==  24) .or. (nproc ==  29)
     & .or.(nproc == 346) .or. (nproc == 347)
     & .or.(nproc == 356) .or. (nproc == 357)) then
c--- In this category, there should be three jets
        if (jets < 3) then
          filterW_bjet=.false.
        endif
        return
      endif

************************************************************************
c--- WQj/ZQj - basic LO process + extra jet in real
************************************************************************
      if ( (nproc == 411) .or. (nproc == 416)
     & .or.(nproc == 321) .or. (nproc == 326)
     & .or.(nproc == 341) .or. (nproc == 351)) then
c--- In this category, there should be at least one b-jet
        if (jets < 1) then
          filterW_bjet=.false.
        endif
c--- check jet identity if there is just 1 jet
c--- (this is not normally the case, unless looking at WQ+X final state)
        if (jets == 1) then
          if ((jetlabel(1) == 'bq').or.(jetlabel(1) == 'ba')) then
            continue
          else
            filterW_bjet=.false.
          endif
      endif
c--- only need to check jet identities if there are just 2 jets
        if (jets == 2) then
        if (((jetlabel(1) == 'bq').and.(jetlabel(2) == 'pp'))
     & .or. ((jetlabel(1) == 'ba').and.(jetlabel(2) == 'pp'))
     & .or. ((jetlabel(1) == 'pp').and.(jetlabel(2) == 'bq'))
     & .or. ((jetlabel(1) == 'pp').and.(jetlabel(2) == 'ba'))) then
          continue
        else
          filterW_bjet=.false.
        endif
        endif
        return
      endif

************************************************************************
c--- WQj/ZQj - basic LO process + extra Q in real
************************************************************************
      if ( (nproc == 312) .or. (nproc == 317)
     & .or.(nproc == 322) .or. (nproc == 327)
     & .or.(nproc == 342) .or. (nproc == 352) ) then
c--- In this category, there should be two jets, one of which is a b
c--- The cut on the Delta_R of the original b-bbar, to split the PS,
c--- is performed in "realint.f", so that neither the real contribution
c--- nor the counterterms are included in that case
c        if ((Ry(q,5,6) < Rbbmin) .and. (isub == 0)) then
c          filterW_bjet=.false.
c        endif
        if (jets < 2) then
          filterW_bjet=.false.
        endif
c--- only need to check jet identities if there are just 2 jets
        if (jets == 2) then
        if (((jetlabel(1) == 'bq').and.(jetlabel(2) == 'pp'))
     & .or. ((jetlabel(1) == 'ba').and.(jetlabel(2) == 'pp'))
     & .or. ((jetlabel(1) == 'pp').and.(jetlabel(2) == 'bq'))
     & .or. ((jetlabel(1) == 'pp').and.(jetlabel(2) == 'ba'))) then
          continue
        else
          filterW_bjet=.false.
        endif
        endif
        return
      endif

************************************************************************
c--- WQj/ZQj - basic LO process + extra Q in real
************************************************************************
      if ((nproc == 441) .or. (nproc == 446)) then
c--- In this category, there should be at least one b-jet
        if (jets < 1) then
          filterW_bjet=.false.
        endif
c--- check jet identity if there is just 1 jet
        if (jets == 1) then
          if ((jetlabel(1) == 'bq').or.(jetlabel(1) == 'ba')) then
            continue
          else
            filterW_bjet=.false.
          endif
      endif
c--- check jet identities if there are just 2 jets
        if (jets == 2) then
          if (((jetlabel(1) == 'bq').or.(jetlabel(2) == 'bq'))
     & .  or. ((jetlabel(1) == 'ba').or.(jetlabel(2) == 'ba'))) then
            continue
          else
            filterW_bjet=.false.
          endif
        endif
        return
      endif


************************************************************************
c--- WQj/ZQj - basic LO process + extra Q in real
************************************************************************
      if ((nproc == 443) .or. (nproc == 448)) then
c--- In this category, there should be a b-jet formed by merging 2 b's
        if (jets < 1) then
          filterW_bjet=.false.
        endif
c--- check jet identity if there is just 1 jet
        if (jets == 1) then
          if (jetlabel(1) == 'bb') then
            continue
          else
            filterW_bjet=.false.
          endif
        endif
c--- check jet identities if there are just 2 jets
        if (jets == 2) then
          if (((jetlabel(1) == 'bb').and.(jetlabel(2) == 'pp'))
     & .  or. ((jetlabel(1) == 'pp').and.(jetlabel(2) == 'bb')))then
            continue
          else
            filterW_bjet=.false.
          endif
        endif
        return
      endif

************************************************************************
c--- To check the ALPGEN approach to producing WQj
************************************************************************
c      if ((nproc == 433) .or. (nproc == 438)) then
c--- In this category, there should be two jets, one of which is a b
c--- jets should not have been merged
c        if ((jets .ne. 2) .or. (jetmerge .eqv. .true.)) then
c          filterW_bjet=.false.
c        endif
c        if (((jetlabel(1) == 'bq').and.(jetlabel(2) == 'ba'))
c     & .or. ((jetlabel(1) == 'ba').and.(jetlabel(2) == 'bq'))) then
c          filterW_bjet=.false.
c        endif
c        return
c      endif

************************************************************************
c--- W(QQ)j
************************************************************************
      if ((nproc == 433) .or. (nproc == 438)) then
c--- In this category, there should be two jets, one of which is a b
c--- also, jets should have been merged, not missed because of pt or y
        if ((jets .ne. 2) .or. (jetmerge .eqv. .false.)) then
          filterW_bjet=.false.
        endif
        if (((jetlabel(1) == 'bq').and.(jetlabel(2) == 'ba'))
     & .or. ((jetlabel(1) == 'ba').and.(jetlabel(2) == 'bq'))) then
          filterW_bjet=.false.
        endif
        return
      endif

************************************************************************
c--- WQjj
************************************************************************
      if ( (nproc == 434) .or. (nproc == 439)
     & .or.(nproc == 324) .or. (nproc == 329) ) then
c--- In this category, there should be three jets
        if (jets < 3) then
          filterW_bjet=.false.
        endif
        return
      endif


      return
      end

