c --- Common block for keeping track of weights, used
c --- if unweighting is selected :
      double precision wtmin,wtmax,newwt
      logical evtgen
      logical unweight
      logical skipnt
      common/maxwt/wtmin,wtmax,newwt,evtgen,skipnt,unweight

c --- Useful local variables where weights are being checked :
      double precision wtabs
