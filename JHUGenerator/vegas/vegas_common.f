      integer mxdim
      parameter(mxdim=30)
      integer ndim,ncall,itmx,nprn
      double precision xl(mxdim),xu(mxdim),acc
      double precision xi(250,mxdim),si,si2,swgt,schi
      integer ndo,it,idum
      logical readin,stopadapt,writeout,stopvegas
      character*(500) ingridfile,outgridfile
      common/gridinfo_logic/readin,stopadapt,writeout,stopvegas
      common/gridinfo_char/ingridfile,outgridfile
      common/bveg1/xl,xu,acc,ndim,ncall,itmx,nprn
      common/bveg2/xi,si,si2,swgt,schi,ndo,it
      common/ranno/idum
