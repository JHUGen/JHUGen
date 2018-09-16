      subroutine genclust2(q,R,qfinal,isub)
c--- this is a wrapper routine for the jet clustering algorithm
c--- which re-routes according to the value of 'algorithm'  to:
c---  ('ktal') genclust_kt.f     for kt clustering
c---  ('ankt') genclust_kt.f     for "anti-kt" clustering
c---  ('cone') genclust_cone.f   for cone algorithm
c---  ('hqrk') genclust_hqrk.f   for simplified heavy-quark algorithm
c---  ('none') to perform no clustering at all
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'bbproc.f'
      include 'process.f'
      include 'part.f'
      include 'nqcdjets.f'
      include 'notag.f'
      character*4 mypart
      double precision q(mxpart,4),qfinal(mxpart,4),
     & qreorder(mxpart,4),R,Rbbmin
      integer isub,i,nu,njetsmin,njetsmax
      logical first
      common/Rbbmin/Rbbmin
      common/mypart/mypart
      data first/.true./
      save first
!$omp threadprivate(first)      
 
      if ((first) .and.
     &   ((nqcdjets .gt. 0).or.(part .eq. 'real').or.(notag.gt.0))) then
        first=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
      write(6,*)
      write(6,*) '*********** Basic jet-defining parameters **********'
      if     (algorithm .eq. 'ktal') then
      write(6,*) '*          (Run II kT clustering algorithm)        *'
      elseif (algorithm .eq. 'ankt') then
      write(6,*) '*     (Anti-kt algorithm - see arXiv:0802.1189)    *'
      elseif (algorithm .eq. 'cone') then
      write(6,*) '*              (Run II cone algorithm)             *'
      elseif (algorithm .eq. 'hqrk') then
      write(6,*) '*        (Simple cone algorithm for W/Z+Q+j)       *'
      elseif (algorithm .eq. 'none') then
      write(6,*) '*             (no clustering algorithm)            *'
      else
      write(6,*)
      write(6,*) 'Invalid selection of algorithm in input file.'
      write(6,*) 'Please select either ktal, cone, hqrk or none'
      stop
      endif
      write(6,*) '*                                                  *'
      write(6,79) ' *     pt(jet)         > ',ptjetmin
      write(6,79) ' *   |pseudo-rap(jet)| > ',etajetmin   
      write(6,79) ' *   |pseudo-rap(jet)| < ',etajetmax   
      if (bbproc) then
        ptbjetmin=max(ptjetmin,ptbjetmin)
        etabjetmax=min(etajetmax,etabjetmax)
      write(6,79) ' *   pt(b-jet)         > ',ptbjetmin
      write(6,79) ' * |pseudo-rap(b-jet)| < ',etabjetmax   
      endif
      if (algorithm .eq. 'hqrk') then
      write(6,79) ' *   b-bbar separation : ',Rbbmin
      write(6,79) ' *        cone size, R : ',R      
      else
      write(6,79) ' * pseudo-cone size, R : ',R
      endif
      write(6,*) '*                                                  *'
      if ((case .eq. 'W_twdk') .or. (case .eq. 'Wtdkay')) then
      write(6,79) ' *   pt(b-jet @ NLO)   < ',ptbjetmin
      write(6,*) '*                                                  *'
      endif
      njetsmin=nqcdjets-notag
      if (inclusive) then
        if( (part.eq.'real') .or. (mypart.eq.'tota')
     &  .or.(mypart.eq.'todk') )then
          njetsmax=nqcdjets+1
        else
          njetsmax=nqcdjets
        endif
      else
        njetsmax=nqcdjets-notag
      endif
      write(6,78) njetsmin,njetsmax
      write(6,*) '****************************************************'
      call flush(6)
      endif
   78 format(' *    Cross-section defined by:  ',i2,' <= jets <=',
     &        i2,'    *')
   79 format(a25,f8.4,'                   *')

      if     (algorithm .eq. 'ktal') then
        call genclust_kt(q,R,qfinal,isub,+1)
      elseif (algorithm .eq. 'ankt') then
        call genclust_kt(q,R,qfinal,isub,-1)
      elseif (algorithm .eq. 'cone') then
        call genclust_cone(q,R,qfinal,isub)
      elseif (algorithm .eq. 'hqrk') then
        call genclust_hqrk(q,R,qfinal,isub)
      elseif (algorithm .eq. 'none') then
        do i=1,mxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        jets=nqcdjets
      if ((part .eq. 'real') .and. (isub .eq. 0)) jets=jets+1
        return
      else
        write(6,*) 'Invalid choice of jet algorithm, must be'
        write(6,*) '   ktal, ankt, cone, hqrk, none'
        stop
      endif

c--- reorder jets for some special cases, to preserve identities of
c--- particles for use in the plotting routines
      if (    (case .eq. 'qq_ttg')
     &   .or. (case .eq. 'tt_bbl') 
     &   .or. (case .eq. 'tt_bbh')
     &   .or. (case .eq. 'tt_ldk') 
     &   .or. (case .eq. 'tt_hdk')
     &   .or. (case .eq. 'tt_udk') 
     &   .or. (case .eq. 'tthWdk')
     &   .or. (case .eq. 'tt_bbu')
     &   .or. (case .eq. '4ftwdk')
     &   .or. (case .eq. 'dk_4ft')
     &   .or. (case .eq. 'qq_ttw')
     &   .or. (case .eq. 'ttwldk')) then
      call jetreorder(qfinal,qreorder,isub)
      do i=1,mxpart
        do nu=1,4
          qfinal(i,nu)=qreorder(i,nu)
        enddo
      enddo
      endif

      return
      end
      
