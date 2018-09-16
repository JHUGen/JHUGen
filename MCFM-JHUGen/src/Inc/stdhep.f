************************************************************************
*  Implementation of STDHEP v4.09 common block                         *
************************************************************************
      integer:: nmxhep
      parameter (nmxhep=4000)
      integer:: nevhep,nhep,isthep,idhep,jmohep,jdahep
      real(dp):: phep,vhep
      common /hepevt/ nevhep, nhep, isthep(nmxhep), idhep(nmxhep),
     &jmohep(2,nmxhep), jdahep(2,nmxhep), phep(5,nmxhep), vhep(4,nmxhep)
