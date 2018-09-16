c --- Common block for keeping track of weights, used
c --- if unweighting is selected :
      real(dp):: wtmax,newwt
      logical:: evtgen
      logical:: unweight
      logical:: skipnt
      integer:: nevtrequested
      common/maxwt/wtmax,newwt,nevtrequested,evtgen,skipnt,unweight

c --- Useful local variables where weights are being checked :
      real(dp):: wtabs
