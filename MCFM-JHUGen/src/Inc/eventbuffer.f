c --- Author: D. Waters, September 2002
c --- Event buffer, used during unweighting to store
c --- the unweighted events generated during a single
c --- VEGAS pass.
      integer buffersize
      parameter(buffersize=10000)
      integer numstored,numused
      double precision eventbuffer(buffersize,mxpart,4)
      double precision wtbuffer(buffersize)
      integer indexlist(buffersize)
      integer pflavbuffer(buffersize),pbarflavbuffer(buffersize)
      common/eventlist/numstored,numused,eventbuffer,wtbuffer,indexlist,
     +     pflavbuffer,pbarflavbuffer
