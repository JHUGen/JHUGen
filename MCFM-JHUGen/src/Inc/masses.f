c      real(dp):: md,mu,ms,mc,mb,mt
c      common/qmass1/md,mu,ms,mc,mb,mt

c      real(dp):: mel,mmu,mtau
c      common/lmass/mel,mmu,mtau

c      real(dp):: hmass,hwidth
c      common/hmass/hmass,hwidth

c      real(dp):: wmass,wwidth
c      common/wmass/wmass,wwidth

c      real(dp):: zmass,zwidth
c      common/zmass/zmass,zwidth

c      real(dp):: twidth
c      common/twidth/twidth

c      real(dp):: tauwidth
c      common/tauwidth/tauwidth

c      real(dp):: mtausq,mcsq,mbsq
c      common/qmassq/mtausq,mcsq,mbsq

      real(dp):: 
     & md,mu,ms,mc,mb,mt,
     & mel,mmu,mtau,
     & hmass,hwidth,
     & wmass,wwidth,
     & zmass,zwidth,
     & twidth,
     & tauwidth,
     & mtausq,mcsq,mbsq
      common/masses/
     & md,mu,ms,mc,mb,mt,
     & mel,mmu,mtau,
     & hmass,hwidth,
     & wmass,wwidth,
     & zmass,zwidth,
     & twidth,
     & tauwidth,
     & mtausq,mcsq,mbsq
!$omp threadprivate(/masses/)
