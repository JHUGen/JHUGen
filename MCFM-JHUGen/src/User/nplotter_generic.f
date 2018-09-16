      subroutine nplotter_generic(p,wt,wt2,switch)
      implicit none
      include 'types.f'
      include 'vegas_common.f'
      include 'bbproc.f'
      include 'clustering.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'kprocess.f'
      include 'removebr.f'
      include 'nodecay.f'
      include 'useet.f'
      include 'plabel.f'
      include 'outputflags.f'
      include 'nqcdjets.f'
      include 'nproc.f'
cz Add single-top b fraction, Z. Sullivan 1/25/05
      real(dp):: bwgt
      common/btagging/ bwgt

cz //

      integer:: n,switch,i5,i6,i7,nu,nplotmax
      integer tag
      real(dp):: DETABB,DPHIZJ,DPHIBB
     & ,ETAB1,ETAB2,ETANOB,ETARAP,ETARAPTWO
     & ,ETARAPTHREE
     & ,ETBIN
     & ,ETdoubleBIN
     & ,GETET
     & ,M56CLUST
     & ,MBB
     & ,PT
     & ,PTB1
     & ,PTB2
     & ,PTNOB
     & ,PTQ1
     & ,ETAQ1
     & ,PTOTHER
     & ,ETAOTHER
     & ,PTTWO
     & ,PTTHREE
     & ,R
     & ,R57
     & ,R67
     & ,RBB
     & ,SWAP
     & ,WT
     & ,WT2
     & ,YRAP
     & ,YRAPTWO
     & ,YRAPTHREE
     & ,DPHI_LL
     & ,M_LL
     & ,MTRANS
     & ,SCUT1
     & ,SCUT2
     & ,PT34_VETO
     & ,PTQ1_VETO,ETAQ1_VETO,PHI56,PTLEADINGB,PTLEADINGNONB
c     & ,HT_STOP
c     & ,MERECON
c     & ,MLBNU
c     & ,QETA
c     & ,RECONCORR

cz Added by Z. Sullivan 1/25/05
cz jet(3) will be filled with pt-ordered jets in all processes
cz ibbar  will hold an index to the extra b in t-channel single-top for
cz          use with bwgt
      integer:: jet(mxpart),jetstart,ibbar,iz,izj,iztmp,
     & inotb,ilight1,ilight2
cz //

      real(dp):: wtbbar,wtnotb,wtlight1,wtlight2
      real(dp):: m34,etmiss,misset,m345678,
     & p(mxpart,4),fphi,HT,etcharm,deltaeta,cosdeltaphi
      real(dp):: eta3,eta4,eta5,eta6,eta7,eta8,eta9,eta10,
     & eta34,eta56,y34
      real(dp):: r45,r56,m345,m348,m567,m678
      real(dp):: pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,
     & pt34,pt56,oldpt(5:7)
      real(dp):: pt345,eta345,y345,y3,y4
      real(dp):: bclustmass,etvec(4),tmp5(4),tmp6(4),tmp7(4)
c      real(dp):: langle,plep(4),plep_wrest(4),pw(4),pw_wrest(4)
      integer:: ssi(4),tmpi,j,k
      real(dp):: ssd(4),tmpd
      integer:: eventpart,ib1,ib2
      logical:: first,jetmerge
      logical:: jetevent
      character*2 ptet
      common/nplotmax/nplotmax
      common/jetmerge/jetmerge
      common/hwwvars/dphi_ll,m_ll,mtrans,scut1,scut2
c      real(dp):: ycut001,ycut005,ycut010,ycut015,ycut020,
c     & y32,y43,z3,z4,z5,z6
      data first/.true./
      save first
ccccc!$omp threadprivate(first,/nplotmax/)
!$omp threadprivate(/jetmerge/)
      
c--- Check for special plotting instructions
c      if (runstring(1:4) == 'cfmt') then
c        call plots_stop_cfmt(p,wt,wt2)
c      return
c      endif

c--- Set up string for pt or Et
      if (useEt) then
        ptet='Et'
      else
        ptet='pt'
      endif
      
      if (first) then
        tag=tagbook
c--- ensure we initialize all possible histograms
        eventpart=npart+3
        eta3=1d3
        pt3=zip
        eta4=1d3
        pt4=zip
        eta5=1d3
        pt5=zip
        eta6=1d3
        pt6=zip
        eta7=1d3
        pt7=zip
        eta8=1d3
        pt8=zip
        eta9=1d3
        pt9=zip
        eta10=1d3
        pt10=zip
        eta34=1d3
        y34=1d3
        pt34=zip
        r45=zip
        eta56=1d3
        pt56=zip        
        m56clust=zip
        r56=zip
        r57=zip
        r67=zip
        misset=zip
        etbin=zip
        mbb=zip
        etab1=1d3
        etab2=1d3
        etanob=1d3
        ptb1=zip
        ptb2=zip
        ptQ1=zip
        etaQ1=zip
        ptother=zip
        etaother=1d3
        ptnob=zip
        rbb=zip
        detabb=1d3
        dphibb=1d3
        m345=zip
        m348=zip
        m567=zip
        m678=zip
        m345678=zip
        HT=zip
      deltaeta=99._dp
      cosdeltaphi=99._dp
      ssi(1)=5
      ssi(2)=6
        jetmerge=.true.
        jets=nqcdjets
        goto 99
      else
        tag=tagplot
      endif

c--- 'eventpart' will contain the number of actual particles that have
c--- a defined momentum. For most processes, this is calculated as follows:
c----  for lowest order and virtual terms switch=0 and eventpart=npart+2
c---   for real events switch=0 and eventpart=npart+2
c---   for real counter-events switch=1 and eventpart=npart+1
c--- There are some processes for which this is not correct and these
c---  are handled with reference to nproc  
      eventpart=npart-switch+2

      if (jets > 0) then
        eventpart=4+jets
      endif
      if (    (kcase==kW_only) .or. (kcase==kZ_only)
     &    .or.(kcase==kggfus0) .or. (kcase==kqg_tbq)) then
        eventpart=4+jets
      elseif (kcase==kbq_tpq) then
        eventpart=5+jets
      elseif ((kcase==kWWqqbr) .or. (kcase==kWWnpol)
     &    .or.(kcase==kWZbbar) .or. (kcase==kZZlept)
     &    .or.(kcase==kHWW_4l) .or. (kcase==kHZZ_4l)
     &    .or.(kcase==kHWW2lq) .or. (kcase==kggWWbx) 
     &    .or.(kcase==kHWW_tb) .or. (kcase==kHWWint)
     &    .or.(kcase==kHWWHpi) .or. (kcase==kggWW4l)
     &    .or.(kcase==kHZZ_tb) .or. (kcase==kHZZint)
     &    .or.(kcase==kHWWjet) .or. (kcase==kqq_HWW)
     &    .or.(kcase==kWW_jet) .or. (kcase==kZZ_jet)
     &    .or.(kcase==kHWWjet) .or. (kcase==kHZZjet)
     &    .or.(kcase==kHWW2jt) .or. (kcase==kHZZ2jt)
     &    .or.(kcase==kHWW3jt) .or. (kcase==kHZZ3jt)) then
        eventpart=6+jets
      elseif (((kcase==kWHbbar) .or. (kcase==kZHbbar))
     &         .and. (removebr)) then
        eventpart=6+jets
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)
     &   .or. (kcase==ktt_bbu) .or. (kcase==ktt_ldk)
     &   .or. (kcase==ktt_udk)) then
        eventpart=min(8+jets,9)
      elseif (kcase==kqq_ttg) then
        eventpart=9
      elseif ((kcase==kW_twdk) .or. (kcase==kWtdkay)
     &    .or.(kcase==kWtbwdk)) then
        eventpart=6+jets
        if (removebr) eventpart=eventpart+1
      elseif ((kcase==kWH__WW) .or. (kcase==kZH__WW)) then
        eventpart=8+jets
      endif
      if (nproc == 73) eventpart=4+jets
      if (kcase==kHi_Zga) eventpart=5+jets
      if (kcase==kqq_ttw) eventpart=11
      
c--- this variable should be set to .true. when the jets are reordered
c---  according to their pt (or Et)
      jetevent=.false.
c--- re-order jets according to pt, for a W/Z/H+jet event        
      if ((kcase==kW_1jet) .or. (kcase==kZ_1jet) .or.
     &    (kcase==kW_2jet) .or. (kcase==kZ_2jet) .or.
     &    (kcase==kW_3jet) .or. (kcase==kZ_3jet) .or.
     &    (kcase==kggfus1) .or. (kcase==kggfus2) .or.
     &    (kcase==kggfus3) .or. (kcase==kqq_Hqq) .or.
     &    (kcase==kqqHqqg) .or. (kcase==kqg_tbq) .or.
     &    (kcase==kgagajj)) then
        jetevent=.true.
        if (jetalgorithm == Rsepcone) then
          if (jets > 0) pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))
          if (jets > 1) pt6=getet(p(6,4),p(6,1),p(6,2),p(6,3))
          if (jets > 2) pt7=getet(p(7,4),p(7,1),p(7,2),p(7,3))
        else
          if (jets > 0) pt5=pt(5,p)
          if (jets > 1) pt6=pt(6,p)
          if (jets > 2) pt7=pt(7,p)        
        endif
        i5=5
        i6=6
        i7=7
      if (jets > 0) oldpt(5)=pt5
      if (jets > 1) oldpt(6)=pt6
      if (jets > 2) oldpt(7)=pt7
c--- sort for 2 jets 
        if (jets == 2) then          
          if (pt6 > pt5) then
            i5=6
            i6=5
          endif
        endif
c--- sort for 3 jets 
        if (jets == 3) then
          if ((pt5 > pt6) .and. (pt5 > pt7)) then
             i5=5
            if (pt6 > pt7) then
              i6=6
              i7=7
            else
              i6=7
              i7=6
            endif
          endif
          if ((pt6 > pt5) .and. (pt6 > pt7)) then
             i5=6
            if (pt5 > pt7) then
              i6=5
              i7=7
            else
              i6=7
              i7=5
            endif
          endif
          if ((pt7 > pt5) .and. (pt7 > pt6)) then
             i5=7
            if (pt5 > pt6) then
              i6=5
              i7=6
            else
              i6=6
              i7=5
            endif
          endif
        endif
c--- perform exchange
        do nu=1,4
          tmp5(nu)=p(i5,nu)
          tmp6(nu)=p(i6,nu)
          tmp7(nu)=p(i7,nu)
        enddo
        do nu=1,4
          p(5,nu)=tmp5(nu)
          p(6,nu)=tmp6(nu)
          p(7,nu)=tmp7(nu)
        enddo
      if (jets > 0) pt5=oldpt(i5)
      if (jets > 1) pt6=oldpt(i6)
      if (jets > 2) pt7=oldpt(i7)
      endif      

      m34=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)

      if (bbproc .and. clustering) then
c--- returns zero cluster mass if two b's are in one jet
        m56clust=sqrt(bclustmass(p))
      else
        m56clust=sqrt(max(zip,
     &   (p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &  -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2))
      endif  

      eta3=etarap(3,p)
      y3=yrap(3,p)
      pt3=pt(3,p)
      eta4=etarap(4,p)
      y4=yrap(4,p)
      pt4=pt(4,p)        
      eta34=etaraptwo(3,4,p)
      y34=yraptwo(3,4,p)
      pt34=pttwo(3,4,p)
      HT=pt3+pt4
      
      if (eventpart > 4) then        
      eta5=etarap(5,p)
      if (jetevent .eqv. .false.) pt5=pt(5,p)
      r45=R(p,4,5)
      m345=sqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     &          -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
      pt345=ptthree(3,4,5,p)
      eta345=etarapthree(3,4,5,p)
      y345=yrapthree(3,4,5,p)
      HT=HT+pt5
      endif
      
      if (eventpart > 5) then
      eta6=etarap(6,p)
      if (jetevent .eqv. .false.) pt6=pt(6,p)
      eta56=etaraptwo(5,6,p)
      eta56=etaraptwo(5,6,p)
      pt56=pttwo(5,6,p)
      r56=R(p,5,6)
      phi56=fphi(5,6,p)
      HT=HT+pt6
      endif

      if (eventpart > 6) then        
      eta7=etarap(7,p)
      if (jetevent .eqv. .false.) pt7=pt(7,p)
      r57=R(p,5,7)
      r67=R(p,6,7)
      m567=sqrt((p(5,4)+p(6,4)+p(7,4))**2-(p(5,1)+p(6,1)+p(7,1))**2
     &          -(p(5,2)+p(6,2)+p(7,2))**2-(p(5,3)+p(6,3)+p(7,3))**2)
      HT=HT+pt7
      endif

      if (eventpart > 7) then        
      eta8=etarap(8,p)
      pt8=pt(8,p)
      m348=sqrt((p(3,4)+p(4,4)+p(8,4))**2-(p(3,1)+p(4,1)+p(8,1))**2
     &          -(p(3,2)+p(4,2)+p(8,2))**2-(p(3,3)+p(4,3)+p(8,3))**2)
      m678=sqrt((p(6,4)+p(7,4)+p(8,4))**2-(p(6,1)+p(7,1)+p(8,1))**2
     &          -(p(6,2)+p(7,2)+p(8,2))**2-(p(6,3)+p(7,3)+p(8,3))**2)
      m345678=sqrt((p(3,4)+p(4,4)+p(5,4)+p(6,4)+p(7,4)+p(8,4))**2
     &             -(p(3,1)+p(4,1)+p(5,1)+p(6,1)+p(7,1)+p(8,1))**2
     &             -(p(3,2)+p(4,2)+p(5,2)+p(6,2)+p(7,2)+p(8,2))**2
     &             -(p(3,3)+p(4,3)+p(5,3)+p(6,3)+p(7,3)+p(8,3))**2)
      HT=HT+pt8
      endif      
           
      if (eventpart > 8) then        
      eta9=etarap(9,p)
      pt9=pt(9,p)
      HT=HT+pt9
      endif
             
      if (eventpart > 9) then        
      eta10=etarap(10,p)
      pt10=pt(10,p)
      HT=HT+pt10
      endif
          
      misset=etmiss(p,etvec)

c--- set-up variables to catch b's
        if (bbproc) then
          if     (jets == 1) then
            write(6,*) 'Error: bbproc set, but only 1 jet in nplotter.f'
            stop
          elseif (jets == 2) then
            mbb=m56clust
            ptb1=pt5
            ptb2=pt6
            etab1=eta5
            etab2=eta6
            rbb=r56
            if (ptb2 > ptb1) then
              swap=ptb1
              ptb1=ptb2
              ptb2=swap
              swap=etab1
              etab1=etab2
              etab2=swap
            endif
          elseif (jets == 3) then
            call getbs(p,ib1,ib2)
            if     (ib1 == 5) then
              ptb1=pt5
              etab1=eta5
            elseif (ib1 == 6) then
              ptb1=pt6
              etab1=eta6
            elseif (ib1 == 7) then
              ptb1=pt7
              etab1=eta7
            endif
            if     (ib2 == 5) then
              ptb2=pt5
              etab2=eta5
            elseif (ib2 == 6) then
              ptb2=pt6
              etab2=eta6
            elseif (ib2 == 7) then
              ptb2=pt7
              etab2=eta7
            endif
            if (ptb2 > ptb1) then
              swap=ptb1
              ptb1=ptb2
              ptb2=swap
              swap=etab1
              etab1=etab2
              etab2=swap
            endif
            if     (ib1+ib2 == 11) then
              mbb=sqrt((p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &                 -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2)
              ptnob=pt7
              etanob=eta7
            elseif (ib1+ib2 == 12) then
              mbb=sqrt((p(5,4)+p(7,4))**2-(p(5,1)+p(7,1))**2
     &                 -(p(5,2)+p(7,2))**2-(p(5,3)+p(7,3))**2)
              ptnob=pt6
              etanob=eta6
            elseif (ib1+ib2 == 13) then
              mbb=sqrt((p(6,4)+p(7,4))**2-(p(6,1)+p(7,1))**2
     &                 -(p(6,2)+p(7,2))**2-(p(6,3)+p(7,3))**2)
              ptnob=pt5
              etanob=eta5
            endif
            rbb=r(p,ib1,ib2)
          else
            if     (nproc == 146) then
              call getbs(p,ib1,ib2)
              ptb1=pt(ib1,p)
              ptb2=pt(ib2,p)
              etab1=etarap(ib1,p)
              etab2=etarap(ib2,p)
              rbb=r(p,ib1,ib2)
            elseif (nproc == 103) then
c--- set b-quarks to 5 and 6 arbitrarily: if you're interested in
c--- this process, these should be set up properly 
              ib1=5
              ib2=6
              ptb1=pt(ib1,p)
              ptb2=pt(ib2,p)
              etab1=etarap(ib1,p)
              etab2=etarap(ib2,p)
              rbb=r(p,ib1,ib2)
            elseif ((kcase==ktt_bbh) .or. (kcase==kqq_ttw)) then
c--- set b-quarks to 5 and 6: this should have been ensured by jetreorder.f
              ib1=5
              ib2=6
              ptb1=pt(ib1,p)
              ptb2=pt(ib2,p)
              etab1=etarap(ib1,p)
              etab2=etarap(ib2,p)
              rbb=r(p,ib1,ib2)
            else
            write(6,*) 'Unforeseen # of jets and b-quarks in nplotter.f'
            stop
            endif
          endif
          detabb=abs(etab1-etab2)
          dphibb=sqrt(abs(rbb**2-detabb**2))
        endif

c--- special variables for ordering of b-quarks in ttH production
      if (kcase==kqq_tth) then
        eventpart=10 ! ensure all histograms are filled
c--- index into pt-ordered list of b-quarks using ssi(n), where n=1,..,4
c--- ordered list is obtained by shuffle-sort
        ssi(1)=5
      ssi(2)=6
      ssi(3)=9
      ssi(4)=10
      do j=1,4
        ssd(j)=pt(ssi(j),p)
      enddo
      do j=4,1,-1
      do k=j-1,1,-1
        if(ssd(j) > ssd(k)) then
          tmpd=ssd(k)
          tmpi=ssi(k)
          ssd(k)=ssd(j)
          ssi(k)=ssi(j)
          ssd(j)=tmpd
          ssi(j)=tmpi
        endif
      enddo
      enddo
c--- example: invariant mass of two highest-pt jets
        mbb=sqrt(max(zip,
     &   (p(ssi(1),4)+p(ssi(2),4))**2-(p(ssi(1),1)+p(ssi(2),1))**2
     &  -(p(ssi(1),2)+p(ssi(2),2))**2-(p(ssi(1),3)+p(ssi(2),3))**2))
c      do j=1,4
c      write(6,*) j,ssd(j),ssi(j)
c      enddo
c      pause
      endif

      if     (kcase==kgQ__ZQ) then
c--- for the Z+b process, catch the highest pt heavy quark
        call getptQ1(pt5,pt6,eta5,eta6,ptQ1,etaQ1,1)
      elseif (kcase==kH_1jet) then
c--- for the H+b process, catch the highest pt heavy quark
        call getptQ1(pt5,pt6,eta5,eta6,ptQ1,etaQ1,1)
        if (jets > 1) then
          ptother=pt5+pt6-ptQ1
          etaother=eta5+eta6-etaQ1
        else
          ptother=-1._dp
          etaother=99._dp
        endif
        if (jets == 1) then
          pt34_veto=pt34
          ptQ1_veto=ptQ1
          etaQ1_veto=etaQ1
        else
          pt34_veto=-1._dp
          ptQ1_veto=-1._dp
          etaQ1_veto=99._dp
        endif
      else
        ptQ1=-1._dp
        etaQ1=99._dp
        ptother=-1._dp
        etaother=99._dp
      endif

c--- for The Z+Q+jet process, make histograms of the leading jet pt
c---  when it is a b-quark and when it is not
      ptleadingb=-1._dp
      ptleadingnonb=-1._dp
      if (kcase==kZ_bjet) then
        if     ((jetlabel(1) == 'bq') .and. (pt5 > pt6)) then
        ptleadingb=pt5
      elseif ((jetlabel(2) == 'bq') .and. (pt6 > pt5)) then
        ptleadingb=pt6
      elseif ((jetlabel(1) == 'pp') .and. (pt5 > pt6)) then
        ptleadingnonb=pt5
      elseif ((jetlabel(2) == 'pp') .and. (pt6 > pt5)) then
        ptleadingnonb=pt6
      endif
      endif
      
      if (nproc == 61) then
        etbin=etdoublebin(pt4,pt5)
      endif

c--- find largest rapidity difference between the jets
      deltaeta=99._dp
      cosdeltaphi=99._dp
      if ((jets == 2) .or. (jets == 3)) then
        if     (jets == 2) then
        i5=5
        i6=6
        elseif (jets == 3) then
          if (abs(eta5-eta6)>max(abs(eta5-eta7),abs(eta6-eta7))) then
          i5=5
          i6=6
        else
          if (abs(eta5-eta7)>abs(eta6-eta7)) then
            i5=5
            i6=7
          else
            i5=6
            i6=7
          endif
        endif
        endif
      deltaeta=abs(etarap(i5,p)-etarap(i6,p))
      cosdeltaphi=(p(i5,1)*p(i6,1)+p(i5,2)*p(i6,2))
     &           /sqrt((p(i5,1)**2+p(i5,2)**2)*(p(i6,1)**2+p(i6,2)**2))
        if (cosdeltaphi < -0.999999999_dp) cosdeltaphi=-1._dp
c      cosdeltaphi=abs(cosdeltaphi)
      endif
        
   99 continue

c--- make plots for H(->WW)+jet paper
c      if (runstring(1:6) == 'hjetww') then
c        call hwwjetplots(eventpart,tag,p,wt,wt2)
c      first=.false.
c      return
c      endif
      
cz Added by Z. Sullivan, 1/25/05
c--- further amended by JC, 19/8/08
cz jet(i) will hold a pt-ordered index to jet momenta.  jet(i)=0 by default
cz This is designed to work for processes where jet momenta come last,
cz  specifically single-top processes 160-179.  Check process.DAT for momenta
cz  ordering before using in other processes.  (E.g. ttbar, WW, ZZ are more
cz  complicated, and need an index sorting method.)
      do iz = 1,mxpart
         jet(iz) = 0
      enddo
      if (jets>0) then
c Figure out where jets start. Uses same algorithm as genclust_kt.f
c--- pick out jets: note that we search to npart+2-isub, to get the
c--- number of particles right. Note that isub=0 for all calls except
c--- the dipole contributions, where isub=1.   
         do iz=3,eventpart
            if ( (plabel(iz) == 'pp') .or. (plabel(iz) == 'pj')
     &           .or.(plabel(iz) == 'bq') .or. (plabel(iz) == 'ba')
     &           .or.(plabel(iz) == 'qj') ) then
               jetstart=iz
               goto 221
            endif
         enddo

 221     do iz=1,jets                        ! loop over jets
            jet(iz)=jetstart + iz-1          ! append jet
c           Compare iz vs. ordered list, and move up if needed
            if (iz>1) then
               do izj = iz,2,-1
                  if (pt(jet(izj),p)>pt(jet(izj-1),p)) then
c                    swap izj,izj-1
                     iztmp=jet(izj-1)
                     jet(izj-1)=jet(izj)
                     jet(izj)=iztmp
                  endif
               enddo
            endif
         enddo
      endif
cz // end jet(3) filling


cz Extract index to b~ (b for t~) in t-channel single-top
cz Events should be plotted with weight = wt*bwgt instead of just wt.
      ibbar=0
      inotb=0
      if (((nproc==161).or.(nproc==166)).and.(jets>0)) then
c         if (bwgt>1.e-10_dp) then    ! some reasonable minimal threshold
c   There are 2 cases: 
c   If merging did not occur, then only if jetlabel(i)=='pp' is it the b~.
            do iz=1,jets
               if     (jetlabel(jet(iz)-jetstart+1)=='pp') then
               ibbar=jet(iz)
               elseif (jetlabel(jet(iz)-jetstart+1)=='qj') then
               inotb=jet(iz)
               endif        
            enddo
c   If merging occurred, then it matters whether the merged jet was already
c    a bq or not.  There is a dirty test for the simple case if removebr=true
            if (jetmerge) then
               if (jetstart==6) ibbar=jet(1)   ! b~ was merged into jet(1)
               if (jetstart==5) then
c          To do this correctly requires changing the clustering routine.
c           Therefore, we are ignoring this case for now.  You should always
c           use removebr = .true. for comparisons to event generators anyway.
               endif
            endif
c         endif  ! end bwgt > 1.e-10_dp
      endif
cz // end finding index to b~ in t-channel single-top
cz // Can now plot pt-ordered jets, and b~ fraction (with wt*bwgt)

      if ((kcase==kbq_tpq) .and. (removebr)) then
      ilight1=0
      ilight2=0
      wtbbar=zip
      wtnotb=zip
      wtlight1=zip
      wtlight2=zip
c--- there are three cases
      if     (jets == 1) then
        if     (jetlabel(1) == 'qj') then
c---   1) only one jet, that is definitely not a b
c---           jetlabel = (qj)
        ibbar=0
c         inotb already set above
        wtnotb=wt
      elseif (jetlabel(1) == 'pp') then
c---   2) only one jet, that could be a b
c---           jetlabel = (pp)
c         ibbar already set above
        inotb=ibbar
        wtbbar=wt*bwgt
        wtnotb=wt*(1._dp-bwgt)
        endif        
c---   3) two jets, one of which could be a b
c---           jetlabel = (qj,pp) OR (pp,qj)
      elseif (jets == 2) then
c       ibbar already set above
c       inotb already set above
        wtbbar=wt*bwgt
      wtnotb=wt*bwgt
        ilight1=jet(1) ! they are ordered according to pt
        ilight2=jet(2) ! they are ordered according to pt
      wtlight1=wt*(1._dp-bwgt)
      wtlight2=wtlight1
      elseif (jets == 0) then
c--- nothing to set: all variables=0 and should only occur when notag=1
        continue
      else
        write(6,*) 'Error: there should be 1 or 2 jets, instead ',jets
      stop
      endif
      if (first) then
c--- on first call, use dummies to ensure histograms are initialized
      ibbar=6
      inotb=7
      ilight1=6
      ilight2=7
      endif
      endif

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histos
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/real(itmx,dp))  
c--- REMOVED - to produce normal histograms as well
c        return    
      endif

c--- Otherwise, fill the histograms 
      n=nextnplot 
      
c--- SPECIAL PLOTS FOR T-CHANNEL SINGLE TOP
      if ((kcase==kbq_tpq) .and. (removebr)) then
c--- RAPIDITY FIRST
c--- leading jet that isn't a b
c---    first account for events with only one non b jet (always occurs)
c---    (always occurs, except if notag=1)
      if (inotb .ne. 0) then
      call bookplot(n,tag,'eta notb 1',
     & etarap(inotb,p),wtnotb,wt2,-8._dp,8._dp,0.5_dp,'lin')
      endif
c---    then account for events with two non b jets (only if ilight1>0)
      if (ilight1 .ne. 0) then
      call bookplot(n,tag,'eta notb 1',
     & etarap(ilight1,p),wtlight1,wt2,-8._dp,8._dp,0.5_dp,'lin')
      endif
      n=n+1
c--- subleading jet that isn't a b (only if ilight1>0)
      if (ilight1 .ne. 0) then
      call bookplot(n,tag,'eta notb 2',
     & etarap(ilight2,p),wtlight2,wt2,-8._dp,8._dp,0.5_dp,'lin')
      endif
      n=n+1
c--- b jet (only if ibbar>0)
      if (ibbar .ne. 0) then
      call bookplot(n,tag,'eta b',
     & etarap(ibbar,p),wtbbar,wt2,-8._dp,8._dp,0.5_dp,'lin')
      endif
      n=n+1

c--- PT SECOND
c--- leading jet that isn't a b
c---    first account for events with only one non b jet (always occurs)
c---    (always occurs, except if notag=1)
      if (inotb .ne. 0) then
      call bookplot(n,tag,'pt notb 1',
     & pt(inotb,p),wtnotb,wt2,zip,200._dp,5._dp,'log')
      endif
c---    then account for events with two non b jets (only if ilight1>0)
      if (ilight1 .ne. 0) then
      call bookplot(n,tag,'pt notb 1',
     & pt(ilight1,p),wtlight1,wt2,zip,200._dp,5._dp,'log')
      endif
      n=n+1
c--- subleading jet that isn't a b (only if ilight1>0)
      if (ilight1 .ne. 0) then
      call bookplot(n,tag,'pt notb 2',
     & pt(ilight2,p),wtlight2,wt2,zip,200._dp,5._dp,'log')
      endif
      n=n+1
c--- b jet (only if ibbar>0)
      if (ibbar .ne. 0) then
      call bookplot(n,tag,'pt b',
     & pt(ibbar,p),wtbbar,wt2,0._dp,200._dp,5._dp,'log')
      endif
      n=n+1
      endif
      
c      if (runstring(1:5) == 'carlo') then
cc--- Special histograms for comparison with C. Oleari's
cc---  e+e- --> QQbg calculation (Durham jet rates)
c        ycut001=zip
c      ycut005=zip
c      ycut010=zip
c      ycut015=zip
c      ycut020=zip
c      if (tag == tagplot) then
c        call durhamalg(p,npart-switch,y32,y43,z3,z4,z5,z6)
c      endif
cc      write(6,'(i3,3e16.6)') switch,y43,y32,wt
c      if ((y32 > 0.01_dp) .and. (y43 < 0.01_dp)) ycut001=1._dp
c      if ((y32 > 0.05_dp) .and. (y43 < 0.05_dp)) ycut005=1._dp
c      if ((y32 > 0.10_dp) .and. (y43 < 0.10_dp)) ycut010=1._dp
c      if ((y32 > 0.15_dp) .and. (y43 < 0.15_dp)) ycut015=1._dp
c      if ((y32 > 0.20_dp) .and. (y43 < 0.20_dp)) ycut020=1._dp
c        call bookplot(n,tag,'ycut=0.01',ycut001,wt,wt2,
c     &                0.5_dp,1.5_dp,1._dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'ycut=0.05',ycut005,wt,wt2,
c     &                0.5_dp,1.5_dp,1._dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'ycut=0.10',ycut010,wt,wt2,
c     &                0.5_dp,1.5_dp,1._dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'ycut=0.15',ycut015,wt,wt2,
c     &                0.5_dp,1.5_dp,1._dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'ycut=0.20',ycut020,wt,wt2,
c     &                0.5_dp,1.5_dp,1._dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'y43      ',y43,wt,wt2,
c     &                0._dp,1._dp,0.05_dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'y43>0.01 ',y43,wt,wt2,
c     &                0.01_dp,1.01_dp,0.05_dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'y43>0.05 ',y43,wt,wt2,
c     &                0.05_dp,1._dp,0.05_dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'y43>0.10 ',y43,wt,wt2,
c     &                0.10_dp,1._dp,0.05_dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'y43>0.15 ',y43,wt,wt2,
c     &                0.15_dp,1._dp,0.05_dp,'lin')
c        n=n+1
c        call bookplot(n,tag,'y43>0.20 ',y43,wt,wt2,
c     &                0.20_dp,1._dp,0.05_dp,'lin')
c        n=n+1
c      endif
      
c      if (runstring(1:4) == 'stop') then
cc--- Special histograms for single top search
c        if (tag == tagplot) then
cc        call stopcuts(p,eventpart,ht,qeta,mlbnu,merecon,reconcorr)
cc these variables should already be in the common-block
c        else
c          ht_stop=-1._dp
c          qeta=99._dp
c          mlbnu=-1._dp
c          merecon=-1._dp
c          reconcorr=-99._dp
c        endif
c      call bookplot(n,tag,'HT',HT_stop,wt,wt2,100._dp,500._dp,10._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'qeta',qeta,wt,wt2,-2.8_dp,3.2_dp,0.4_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'mlbnu',mlbnu,wt,wt2,0._dp,300._dp,10._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'reconcorr?',reconcorr,
c     & wt,wt2,-1._dp,1._dp,1._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'MEn recon',merecon,
c     & wt,wt2,0._dp,300._dp,10._dp,'lin')
c      n=n+1
c      endif
      
c      if (runstring(1:3) == 'hww') then
cc--- Special histograms for H->WW search
c        if (tag == tagplot) then
cc these variables should already be in the common-block
c        else
c          dphi_ll=-1._dp
c          m_ll=-1._dp
c          mtrans=-1._dp
c          scut1=-1._dp
c          scut2=-1._dp
c        endif
c      call bookplot(n,tag,'dphi_ll',dphi_ll,wt,wt2,
c     &              0._dp,3.142_dp,0.1571_dp,'lin')
cc      write(6,*) switch,dphi_ll,wt,wt2
c      n=n+1
c      call bookplot(n,tag,'m_ll',m_ll,wt,wt2,0._dp,200._dp,5._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'mtrans',mtrans,wt,wt2,0._dp,250._dp,5._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'scut1',scut1,wt,wt2,0.5_dp,1.5_dp,1._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'scut2',scut2,wt,wt2,0.5_dp,1.5_dp,1._dp,'lin')
c      n=n+1
c      endif

c--- calculate lepton 3 angle in W rest frame
c      do nu=1,4
c      pw(nu)=p(3,nu)+p(4,nu)
c      plep(nu)=p(3,nu)
c      pw_wrest(nu)=0._dp
c      enddo
c      pw_wrest(4)=sqrt(pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2)
c      call boostx(plep,pw,pw_wrest,plep_wrest)
c      langle=plep_wrest(3)/
c     & sqrt(plep_wrest(1)**2+plep_wrest(2)**2+plep_wrest(3)**2)
cc      langle=acos(langle)
c      call bookplot(n,tag,'lepton 3',langle,wt,wt2,
c     & -1._dp,1.01_dp,0.1_dp,'lin')
c      n=n+1
            
c--- calculate lepton 4 angle in W rest frame
c      do nu=1,4
c      pw(nu)=p(3,nu)+p(4,nu)
c      plep(nu)=p(4,nu)
c      pw_wrest(nu)=0._dp
c      enddo
c      pw_wrest(4)=sqrt(pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2)
c      call boostx(plep,pw,pw_wrest,plep_wrest)
c      langle=plep_wrest(3)/
c     & sqrt(plep_wrest(1)**2+plep_wrest(2)**2+plep_wrest(3)**2)
cc      langle=acos(langle)
c      call bookplot(n,tag,'lepton 4',langle,wt,wt2,
c     & -1._dp,1.01_dp,0.1_dp,'lin')
c      n=n+1

c--- special plots for W+c
      if (kcase==kW_cjet) then
        etcharm=-1._dp
        if     (jets == 1) then
        if ((jetlabel(1) == 'bq') .or. (jetlabel(1) == 'ba')) then
          etcharm=sqrt(p(5,1)**2+p(5,2)**2)
     &                  *p(5,4)/sqrt(p(5,1)**2+p(5,2)**2+p(5,3)**2)
          endif
      elseif (jets == 2) then
        if ((jetlabel(1) == 'bq') .or. (jetlabel(1) == 'ba')) then
          etcharm=sqrt(p(5,1)**2+p(5,2)**2)
     &                  *p(5,4)/sqrt(p(5,1)**2+p(5,2)**2+p(5,3)**2)
          endif
        if ((jetlabel(2) == 'bq') .or. (jetlabel(2) == 'ba')) then
          etcharm=sqrt(p(6,1)**2+p(6,2)**2)
     &                  *p(6,4)/sqrt(p(6,1)**2+p(6,2)**2+p(6,3)**2)
          endif
      endif
      if ((etcharm < 0._dp) .and. (tag == tagplot)) then
        write(6,*) 'Error in nplotter: no charm jet for Et plotting'
        stop
      endif
        call bookplot(n,tag,'Et_c',etcharm,wt,wt2,0._dp,100._dp,10._dp,'log')
        n=n+1
      endif
            
c--- special plots for W+cc~ (sign of quark works for W+ only)
      if (kcase==kWbbmas) then
        etcharm=-1._dp
        if     (jets == 1) then
        if (jetlabel(1) == 'ba') then
          etcharm=sqrt(p(5,1)**2+p(5,2)**2)
     &                  *p(5,4)/sqrt(p(5,1)**2+p(5,2)**2+p(5,3)**2)
          endif
      elseif (jets == 2) then
        if (jetlabel(1) == 'ba') then
          etcharm=sqrt(p(5,1)**2+p(5,2)**2)
     &                  *p(5,4)/sqrt(p(5,1)**2+p(5,2)**2+p(5,3)**2)
          endif
        if (jetlabel(2) == 'ba') then
          etcharm=sqrt(p(6,1)**2+p(6,2)**2)
     &                  *p(6,4)/sqrt(p(6,1)**2+p(6,2)**2+p(6,3)**2)
          endif
      endif
c--- Not an error here: only plot the quark correlated with the W
c      if ((etcharm < 0._dp) .and. (tag == tagplot)) then
c        write(6,*) 'Error in nplotter: no charm jet for Et plotting'
c        stop
c      endif
        call bookplot(n,tag,'Et_c',etcharm,wt,wt2,0._dp,100._dp,10._dp,'log')
        n=n+1
      endif

c--- special plots for computing ttbar asymmetry 

c--- inclusive case
      call bookplot(n,tag,'y(t)-y(tbar)',y3-y4,wt,wt2,
     & -4._dp,4._dp,4._dp,'lin')
      n=n+1

c--- for bins of |y(t)-y(tbar)| of (0,0.5), (0.5,1), (1,1.5)
      call bookplot(n,tag,'y(t)-y(tbar)',y3-y4,wt,wt2,
     & -1.5_dp,1.5_dp,0.5_dp,'lin')
      n=n+1

c--- for bins of |y(t)-y(tbar)| of (1.5,4.5)
      call bookplot(n,tag,'y(t)-y(tbar)',y3-y4,wt,wt2,
     & -4.5_dp,4.5_dp,3._dp,'lin')
      n=n+1

    
c--- added extra plot here, for the angle analysis of G. Hesketh et al.
      if ((kcase==kZ_1jet) .or. (kcase==kZ_2jet)) then
      dphizj=atan2(p(3,1)+p(4,1),p(3,2)+p(4,2))-atan2(p(5,1),p(5,2))
      if (dphizj > pi) dphizj=twopi-dphizj
      if (dphizj < -pi) dphizj=twopi+dphizj
      call bookplot(n,tag,'dphi_zj',dphizj,wt,wt2,-3.1416_dp,3.1416_dp,
     &              1.5708e-1_dp,'log')
      n=n+1
      endif
      
      call bookplot(n,tag,'HT',HT,wt,wt2,0._dp,500._dp,20._dp,'lin')
      n=n+1
c --- Histograms to monitor exclusive/inclusive cross-sections:
      if ((jets == nqcdjets)  .or. (tag == tagbook)) then
      call bookplot(n,tag,'= #LO j ',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
      endif
      n=n+1
      if ((jets >= nqcdjets)  .or. (tag == tagbook)) then
      call bookplot(n,tag,'>= #LO j',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
      endif
      n=n+1
      if (nodecay .eqv. .false.) then
      call bookplot(n,tag,'eta3',eta3,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y3',y3,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'3',pt3,wt,wt2,0._dp,150._dp,5._dp,'log')
      n=n+1      
      call bookplot(n,tag,'eta4',eta4,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y4',y4,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'4',pt4,wt,wt2,25._dp,50._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0._dp,150._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0._dp,600._dp,20._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta34',eta34,wt,wt2,-4._dp,4._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-4._dp,4._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'34',pt34,wt,wt2,0._dp,100._dp,2._dp,'log')
      call ebookplot(n,tag,pt34,wt)
      n=n+1
      call bookplot(n,tag,ptet//'34',pt34,wt,wt2,0._dp,40._dp,2._dp,'log')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,0._dp,200._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      endif
      call bookplot(n,tag,'misset',misset,wt,wt2,0._dp,100._dp,2._dp,'lin')
      n=n+1
      
      if (eventpart > 4) then
      call bookplot(n,tag,'eta5',eta5,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta5',eta5,wt,wt2,-10._dp,10._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'5',pt5,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
c      call bookplot(n,tag,ptet//'5',pt5,wt,wt2,0._dp,100._dp,2._dp,'log')
c      n=n+1
      call bookplot(n,tag,'pt345',pt345,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta345',eta345,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y345',y345,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'m345',m345,wt,wt2,0._dp,200._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'R45',R45,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1

c--- The selected heavy quark jet for Z+b and H+b
c      call bookplot(n,tag,'pt34_veto',
c     & pt34_veto,wt,wt2,0._dp,200._dp,10._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'ptQ_central',ptQ1,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'ptQ_veto',
c     & ptQ1_veto,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'etaQ_central',etaQ1,
c     & wt,wt2,-5._dp,5._dp,0.4_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'etaQ_veto',
c     & etaQ1_veto,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
c      n=n+1
c--- The other jet for H+b (which is a heavy quark for process 143)
c      call bookplot(n,tag,'pt other',ptother,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1

      endif
      
      if (eventpart > 5) then
c      call bookplot(n,tag,'eta other',etaother,
c     & wt,wt2,-5._dp,5._dp,0.4_dp,'lin')
c      n=n+1
      call bookplot(n,tag,'eta6',eta6,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'6',pt6,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
c      call bookplot(n,tag,ptet//'6',pt6,wt,wt2,0._dp,200._dp,2._dp,'log')
c      n=n+1
      call bookplot(n,tag,'eta56',eta56,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'56',pt56,wt,wt2,10._dp,250._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'r56',r56,wt,wt2,0._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'phi56',phi56,wt,wt2,0._dp,3.14_dp,0.314_dp,'lin')
      n=n+1
      call bookplot(n,tag,'m56',m56clust,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'m56',m56clust,wt,wt2,0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta5-eta6',eta5-eta6,
     & wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      endif

c      if (bbproc) then
cc--- Leading b Et, 5 GeV bins from 15 to 200
c      call bookplot(n,tag,'bjet1 Et',ptb1,wt,wt2,15._dp,200._dp,5._dp,'log')
c      n=n+1
cc--- Second b Et, 5 GeV bins from 15 to 200
c      call bookplot(n,tag,'bjet2 Et',ptb2,wt,wt2,15._dp,200._dp,5._dp,'log')
c      n=n+1
cc--- Dijet invariant mass, 2 b-jets, 10 GeV bins from 0 to 250
c      call bookplot(n,tag,'bb mass',mbb,wt,wt2,0._dp,350._dp,10._dp,'log')
c      n=n+1
cc--- Delta_R(b,b), 2 b-jets, bins of 0.2 from 0.35 to 4.95
c      call bookplot(n,tag,'deltaRbb',rbb,
c     & wt,wt2,0.35_dp,4.95_dp,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'etab1',etab1,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,ptet//'b1',ptb1,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'etab2',etab2,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,ptet//'b2',ptb2,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'mbb',mbb,wt,wt2,0._dp,200._dp,10._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'mbb',mbb,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'rbb',rbb,wt,wt2,0._dp,4._dp,0.1_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'detabb',detabb,wt,wt2,0._dp,4._dp,0.1_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'dphibb',dphibb,wt,wt2,0._dp,4._dp,0.1_dp,'lin')
c      n=n+1
c      endif

      if (nqcdjets >= 2) then
      call bookplot(n,tag,'delta(eta) >= 2 jets',deltaeta,
     & wt,wt2,0._dp,10._dp,0.5_dp,'lin')
      n=n+1      
      endif

      if (eventpart > 6) then

c      if (bbproc) then
cc--- Non b-jet Et, 5 GeV bins from 15 to 200
c      call bookplot(n,tag,'non-b Et',ptnob,wt,wt2,15._dp,200._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'etanob',etanob,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,ptet//'nob',ptnob,wt,wt2,0._dp,200._dp,5._dp,'log')
c      n=n+1
c      endif

      call bookplot(n,tag,'eta7',eta7,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'7',pt7,wt,wt2,0._dp,100._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'7',pt7,wt,wt2,0._dp,600._dp,20._dp,'log')
      n=n+1
      call bookplot(n,tag,'r57',r57,wt,wt2,0._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'r67',r67,wt,wt2,0._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'m567',m567,wt,wt2,108._dp,308._dp,2._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta diff',eta7-(eta5+eta6)/2._dp,
     & wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      endif      

      if (eventpart > 7) then
      call bookplot(n,tag,'eta8',eta8,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'8',pt8,wt,wt2,0._dp,100._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'m348',m348,wt,wt2,108._dp,308._dp,2._dp,'log')
      n=n+1
      call bookplot(n,tag,'m678',m678,wt,wt2,108._dp,308._dp,2._dp,'log')
      n=n+1
      call bookplot(n,tag,'m345678',m345678,wt,wt2,0._dp,600._dp,10._dp,'lin')
      n=n+1
c--- extra plots for azimuthal angle difference between two charged
c--- leptons, in the top decay case  
      dphi_ll=(p(4,1)*p(7,1)+p(4,2)*p(7,2))
     &     /sqrt((p(4,1)**2+p(4,2)**2)*(p(7,1)**2+p(7,2)**2))
      if (dphi_ll > +0.9999999_dp) dphi_ll=+1._dp
      if (dphi_ll < -0.9999999_dp) dphi_ll=-1._dp
      dphi_ll=acos(dphi_ll)  
      call bookplot(n,tag,'dphi_ll',dphi_ll,wt,wt2,
     &              0._dp,3.142_dp,0.1571_dp,'lin')
      n=n+1
      if ((m345678 < 400._dp) .or. (tag == tagbook)) then
      call bookplot(n,tag,'dphi_ll + cut mtt',dphi_ll,wt,wt2,
     &              0._dp,3.142_dp,0.1571_dp,'lin')
      endif
      n=n+1
      endif      

      if (eventpart > 8) then
      call bookplot(n,tag,'eta9',eta9,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'9',pt9,wt,wt2,0._dp,100._dp,5._dp,'lin')
      n=n+1
      endif      

      if (eventpart > 9) then
      call bookplot(n,tag,'eta10',eta10,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,ptet//'10',pt10,wt,wt2,0._dp,100._dp,5._dp,'lin')
      n=n+1
      endif      

c--- special plots for the ttH process
      if (kcase==kqq_tth) then
      call bookplot(n,tag,'mb1b2',mbb,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'ptb1',pt(ssi(1),p),
     &              wt,wt2,0._dp,100._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ptb2',pt(ssi(2),p),
     &              wt,wt2,0._dp,100._dp,5._dp,'lin')
      n=n+1
      endif

c--- additional plot for WW process      
      if (nproc == 61) then
      call bookplot(n,tag,'etbin',etbin,wt,wt2,0.5_dp,25.5_dp,1._dp,'lin')
      n=n+1
      endif
      
c--- special plots for Z+b jet processes      
      if (kcase==kZ_bjet) then
      call bookplot(n,tag,ptet//' b leading',ptleadingb,
     &               wt,wt2,0._dp,500._dp,20._dp,'log')
      n=n+1
      call bookplot(n,tag,ptet//' non-b leading',ptleadingnonb,
     &               wt,wt2,0._dp,500._dp,20._dp,'log')
      n=n+1
      endif
      
c--- plots useful for processes with 2 or more jets only
      if (nqcdjets >= 2) then
      call bookplot(n,tag,'delta(eta) >= 2 jets',deltaeta,
     & wt,wt2,0._dp,10._dp,0.5_dp,'lin')
      n=n+1      
      if ((jets == 2) .or. (tag == tagbook)) then
      call bookplot(n,tag,'delta(eta) - 2 jets',deltaeta,
     & wt,wt2,0._dp,10._dp,0.5_dp,'lin')
      endif
      n=n+1      
      if ((jets == 2) .or. (tag == tagbook)) then
      call bookplot(n,tag,'for phi: delta(eta) - 2 jets',deltaeta,
     & wt*cosdeltaphi,wt2,0._dp,10._dp,0.5_dp,'lin')
      endif
      n=n+1      
      if ((jets == 2) .or. (tag == tagbook)) then
      call bookplot(n,tag,'cos delta(phi) - 2 jets',cosdeltaphi,
     & wt,wt2,0._dp,1._dp,0.05_dp,'lin')
      endif
      n=n+1      
      if ((jets == 3) .or. (tag == tagbook)) then
      call bookplot(n,tag,'delta(eta) - 3 jets',deltaeta,
     & wt,wt2,0._dp,10._dp,0.5_dp,'lin')
      endif
      n=n+1      
      if ((jets == 3) .or. (tag == tagbook)) then
      call bookplot(n,tag,'for phi: delta(eta) - 3 jets',deltaeta,
     & wt*cosdeltaphi,wt2,0._dp,10._dp,0.5_dp,'lin')
      endif
      n=n+1      
      if ((jets == 3) .or. (tag == tagbook)) then
      call bookplot(n,tag,'cos delta(phi) - 3 jets',cosdeltaphi,
     & wt,wt2,0._dp,1._dp,0.05_dp,'lin')
      endif
      n=n+1
      endif
            
      n=n-1
      
c--- ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
cz      
      bwgt=0._dp  ! for safety
cz //

      return 
      end
      
      subroutine cross(p,i,j,r)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: i,j
      real(dp):: p(mxpart,4),r(3)
      
      r(1)=p(i,2)*p(j,3)-p(j,2)*p(i,3)
      r(2)=p(i,3)*p(j,1)-p(j,3)*p(i,1)
      r(3)=p(i,1)*p(j,2)-p(j,1)*p(i,2)
      
      return
      end
            
      function fphi(n1,n2,p)
      implicit none
      include 'types.f'
      real(dp):: fphi
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: n1,n2
      real(dp):: p(mxpart,4)
    
      fphi=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)
      fphi=fphi/sqrt(p(n1,1)**2+p(n1,2)**2)
      fphi=fphi/sqrt(p(n2,1)**2+p(n2,2)**2)
      if     (fphi > +0.9999999_dp) then
        fphi=0._dp
      elseif (fphi < -0.9999999_dp) then
        fphi=pi
      else
        fphi=acos(fphi)
      endif

      return
      end
                     
      function coslpairet(n1,n2,nm1,nm2,p)
      implicit none
      include 'types.f'
      real(dp):: coslpairet
      include 'nf.f'
      include 'mxpart.f'
      integer:: n1,n2,nm1,nm2,i
      real(dp):: p(mxpart,4),misset(4),pp(4)
      
      if (nm2 == 0) then
        do i=1,4
          misset(i)=p(nm1,i)
        enddo
      else
        do i=1,4
          misset(i)=p(nm1,i)+p(nm2,i)
        enddo
      endif
      
      do i=1,4
        pp(i)=p(n1,i)+p(n2,i)
      enddo
      
      coslpairet=pp(1)*misset(1)+pp(2)*misset(2)
      coslpairet=coslpairet/sqrt(pp(1)**2+pp(2)**2)
      coslpairet=coslpairet/sqrt(misset(1)**2+misset(2)**2)
      
      return
      end
            
      function deltar(i,j,p)
      implicit none
      include 'types.f'
      real(dp):: deltar
      include 'constants.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),phi1,phi2,etarap,dphi
      integer:: i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(p(j,1),p(j,2))
      dphi=phi1-phi2
      if (dphi > pi) dphi=twopi-dphi
      if (dphi < -pi) dphi=twopi+dphi
      deltar=(etarap(i,p)-etarap(j,p))**2+dphi**2
      deltar=sqrt(deltar)
      
      return
      end
      
      subroutine checkmaxhisto(n)
      implicit none
c--- ensure the built-in maximum number of histograms is not exceeded    
      include 'types.f'
      include 'histo.f'

      integer:: n
      
      if (n > maxhisto) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > ',maxhisto,', which is the built-in maximum.'
      write(6,*) 'To use more histograms, change the value of the'
      write(6,*) 'constant MAXHISTO in src/Inc/nplot.f then do:'
      write(6,*)
      write(6,*) ' make clean; make        to recompile from scratch.'
        write(6,*)
      stop
      endif
      
      return
      end
      
