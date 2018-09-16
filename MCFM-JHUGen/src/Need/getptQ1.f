      subroutine getptQ1(pt5,pt6,eta5,eta6,ptQ1,etaQ1,ptoreta)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: ptoreta
      real(dp):: pt5,pt6,eta5,eta6,ptQ1,ptQ2,etaQ1,etaQ2

c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this
c--- If two b-quarks are found, the one returned is either:
c---     the one with the highest pt (ptoreta=1)
c---     the most central one (ptoreta=2)

      if (jets == 1) then
        if ((jetlabel(1) == 'bq') .or. (jetlabel(1) == 'ba')) then
          ptQ1=pt5
          etaQ1=eta5
          return
        else
          write(6,*) 'Error in getptQ1: only 1 jet and it'
          write(6,*) ' is not a heavy quark!'
          stop
        endif
      endif

      if (jets .ne. 2) then
        write(6,*) 'Error in getptQ1: strange number of jets, ',jets,'!'
        stop
      endif

c--- now we know that we have 2 jets
      ptQ1=-1._dp
      ptQ2=-1._dp
      etaQ1=99._dp
      etaQ2=99._dp

      if ((jetlabel(1) == 'bq') .or. (jetlabel(1) == 'ba')) then
        ptQ1=pt5
        etaQ1=eta5
      endif
      if ((jetlabel(2) == 'bq') .or. (jetlabel(2) == 'ba')) then
        ptQ2=pt6
        etaQ2=eta6
      endif

      if     (ptoreta == 1) then
        ptQ1=max(ptQ1,ptQ2)
      elseif (ptoreta == 2) then
        if (abs(etaQ2) < abs(etaQ1)) then
          ptQ1=ptQ2
          etaQ1=etaQ2
        endif
      else
        write(6,*) 'The value of ptoreta in getptQ1.f is incorrect'
        stop
      endif

      if (ptQ1 < 0._dp) then
        write(6,*) 'Error in getptQ1: 2 jets, but no heavy quarks!'
        stop
      endif

      return
      end
