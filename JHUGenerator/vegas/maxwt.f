c --- Common block for keeping track of weights, used
c --- if unweighting is selected :
      double precision wtmax,newwt
      logical evtgen
      logical unweight
      logical skipnt
      common/maxwt/wtmax,newwt,evtgen,skipnt,unweight

c --- Useful local variables where weights are being checked :
      double precision wtabs
