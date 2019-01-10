c      double precision md,mu,ms,mc,mb,mt
c      common/qmass1/md,mu,ms,mc,mb,mt

c      double precision mel,mmu,mtau
c      common/lmass/mel,mmu,mtau

c      double precision hmass,hwidth
c      common/hmass/hmass,hwidth

c      double precision wmass,wwidth
c      common/wmass/wmass,wwidth

c      double precision zmass,zwidth
c      common/zmass/zmass,zwidth

c      double precision twidth
c      common/twidth/twidth

c      double precision tauwidth
c      common/tauwidth/tauwidth

c      double precision mtausq,mcsq,mbsq
c      common/qmassq/mtausq,mcsq,mbsq

      double precision 
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
