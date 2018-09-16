      subroutine genclust2(q,R,qfinal,isub)
      implicit none
      include 'types.f'
c--- this is a wrapper routine for the jet clustering algorithm
c--- which re-routes according to the value of 'algorithm'  to:
c---  ('ktal') genclust_kt.f     for kt clustering
c---  ('ankt') genclust_kt.f     for "anti-kt" clustering
c---  ('cone') genclust_cone.f   for cone algorithm
c---  ('hqrk') genclust_hqrk.f   for simplified heavy-quark algorithm
c---  ('none') to perform no clustering at all

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'clustering.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'bbproc.f'
      include 'kprocess.f'
      include 'kpart.f'
      include 'nqcdjets.f'
      include 'notag.f'
      include 'first.f'
      include 'mpicommon.f'
      integer mykpart
      real(dp):: q(mxpart,4),qfinal(mxpart,4),
     & qreorder(mxpart,4),R,Rbbmin
      integer:: isub,i,nu,njetsmin,njetsmax
      common/Rbbmin/Rbbmin
      common/mykpart/mykpart

      if ((first) .and.
     &   ((nqcdjets > 0).or.(kpart==kreal).or.(notag>0))) then
        first=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
!$omp master
      if (rank == 0) then
      write(6,*)
      write(6,*) '*********** Basic jet-defining parameters **********'
      if     (jetalgorithm == kt) then
      write(6,*) '*          (Run II kT clustering algorithm)        *'
      elseif (jetalgorithm == antikt) then
      write(6,*) '*     (Anti-kt algorithm - see arXiv:0802.1189)    *'
      elseif (jetalgorithm == Rsepcone) then
      write(6,*) '*              (Run II cone algorithm)             *'
      elseif (jetalgorithm == hqrk) then
      write(6,*) '*        (Simple cone algorithm for W/Z+Q+j)       *'
      elseif (jetalgorithm == noclustering) then
      write(6,*) '*             (no clustering algorithm)            *'
      else
      write(6,*)
      write(6,*) 'Invalid selection of algorithm in input file.'
      write(6,*) 'Please select either ktal, ankt, cone, hqrk or none'
      stop
      endif
      write(6,*) '*                                                  *'
      write(6,79) ' *     pt(jet)         > ',ptjetmin
      write(6,79) ' *   |pseudo-rap(jet)| > ',etajetmin
      write(6,79) ' *   |pseudo-rap(jet)| < ',etajetmax
      endif
!$omp end master
      if (bbproc) then
        ptbjetmin=max(ptjetmin,ptbjetmin)
        etabjetmax=min(etajetmax,etabjetmax)
      write(6,79) ' *   pt(b-jet)         > ',ptbjetmin
      write(6,79) ' * |pseudo-rap(b-jet)| < ',etabjetmax
      endif
      if (jetalgorithm == hqrk) then
      write(6,79) ' *   b-bbar separation : ',Rbbmin
      write(6,79) ' *        cone size, R : ',R
      else
!$omp master
      if (rank == 0) then
      write(6,79) ' * pseudo-cone size, R : ',R
      endif
!$omp end master
      endif
!$omp master
      if (rank == 0) then
      write(6,*) '*                                                  *'
      endif
!$omp end master
      if ((kcase==kW_twdk) .or. (kcase==kWtdkay)) then
      write(6,79) ' *   pt(b-jet @ NLO)   < ',ptbjetmin
      write(6,*) '*                                                  *'
      endif
      njetsmin=nqcdjets-notag
      if (inclusive) then
        if( (kpart==kreal) .or. (mykpart==ktota)
     &  .or.(mykpart==ktodk) )then
          njetsmax=nqcdjets+1
        elseif (mykpart==knnlo) then
          njetsmax=nqcdjets+2
        else
          njetsmax=nqcdjets
        endif
      else
        njetsmax=nqcdjets-notag
      endif
!$omp master
      if (rank == 0) then
      write(6,78) njetsmin,njetsmax
      write(6,*) '****************************************************'
      call flush(6)
      endif
!$omp end master
      endif
   78 format(' *    Cross-section defined by:  ',i2,' <= jets <=',
     &        i2,'    *')
   79 format(a25,f8.4,'                   *')

      if     (jetalgorithm == kt) then
        call genclust_kt(q,R,qfinal,isub,+1)
      elseif (jetalgorithm == antikt) then
        call genclust_kt(q,R,qfinal,isub,-1)
      elseif (jetalgorithm == Rsepcone) then
        call genclust_cone(q,R,qfinal,isub)
      elseif (jetalgorithm == hqrk) then
        call genclust_hqrk(q,R,qfinal,isub)
      elseif (jetalgorithm == noclustering) then
        do i=1,mxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        jets=nqcdjets
      if ((kpart==kreal) .and. (isub == 0)) jets=jets+1
        return
      else
        write(6,*) 'Invalid choice of jet algorithm, must be'
        write(6,*) '   ktal, ankt, cone, hqrk, none'
        stop
      endif

c--- reorder jets for some special cases, to preserve identities of
c--- particles for use in the plotting routines
      if (    (kcase==kqq_ttg)
     &   .or. (kcase==ktt_bbl)
     &   .or. (kcase==ktt_bbh)
     &   .or. (kcase==ktt_ldk)
     &   .or. (kcase==ktt_hdk)
     &   .or. (kcase==ktt_udk)
     &   .or. (kcase==ktthWdk)
     &   .or. (kcase==ktt_bbu)
     &   .or. (kcase==k4ftwdk)
     &   .or. (kcase==kdk_4ft)
     &   .or. (kcase==kqq_ttw)
     &   .or. (kcase==kttwldk)) then
      call jetreorder(qfinal,qreorder,isub)
      do i=1,mxpart
        do nu=1,4
          qfinal(i,nu)=qreorder(i,nu)
        enddo
      enddo
      endif

      return
      end

