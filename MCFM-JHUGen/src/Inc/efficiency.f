c--- common block for the calculation of cut efficiencies
c
c--- ntotshot : total number of shots in the run
c--- njetzero : number of shots that failed the jet cuts
c--- ncutzero : number of shots that passed the jet cuts
c---            but then failed the process-specific cuts
c--- ntotzero : total number of events that automatically
c---            returned zero weight. Should be approximately
c---            njetzero+ncutzero, with a small extra number
c---            (dependent on 'cutoff') due to 'masscuts' and 'smalls'

      integer(kind=8) njetzero,ncutzero,ntotzero,ntotshot
      common/efficiency/njetzero,ncutzero,ntotzero,ntotshot
