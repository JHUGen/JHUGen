      subroutine writeinput(unitno,lstring,rstring,tag)
      implicit none
      include 'types.f'
c--- This routine echoes the line in the input file specified by the
c--- input parameter 'tag'; if 'tag' is set to 'WRITEALL' then all
c--- of the input file is echoed.
c--- Output is written to the unit 'unitno', with lines bracketed
c--- by the strings 'lstring' and 'rstring'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'maxwt.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'flags.f'
      include 'clustering.f'
      include 'anomcoup.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'lhapdf.f'
      include 'pdlabel.f'
      include 'removebr.f'
      include 'dynamicscale.f'
      include 'stopscales.f'
      include 'vanillafiles.f'
      include 'alfacut.f'
      include 'betacut.f'
      include 'verbose.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'frag.f'
      include 'initialscales.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'kpart.f'
      include 'anomHiggs.f'
      include 'anom_higgs.f'
      include 'vdecayid.f'
      include 'runstring.f'
      include 'energy.f'
      include 'nproc.f'
      include 'taucut.f'
      include 'iterat.f'
      include 'mpicommon.f'
c--- APPLgrid - flag using grid
c      include 'ptilde.f'
c      include 'APPLinclude.f'
c--- APPLgrid - end
      character*(*) tag,lstring,rstring
      character*72 f93,f94,f95,f96,f97,f98,f99
      character*15 kpartstring
      logical:: dryrun,makecuts,writeall,spira,writerefs
      integer:: unitno, nmin,nmax
      integer:: ih1,ih2,origij
      integer:: NPTYPE,NGROUP,NSET
      real(dp):: rtsmin,Rcut

      common/writerefs/writerefs
      common/spira/spira
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin

      common/density/ih1,ih2
      common/dryrun/dryrun

      common/pdflib/NPTYPE,NGROUP,NSET

      common/Rcut/Rcut
      common/makecuts/makecuts

      common/origij/origij

c--- f93 scientific format
      f93='('''//lstring//''',es20.2E2,12x,''['',a,'']'','''
     & //rstring//''')'
c--- f94 integer::.XXYY format
      f94='('''//lstring//''',11x,i4,''.'',2a2,12x,''['',a,'']'','''
     & //rstring//''')'
c--- f95 2xfloating point format
      f95='('''//lstring//''',f8.3,'','',f8.3,16x,''['',a,'']'','''
     & //rstring//''')'
c--- f96 character format
      f96='('''//lstring//''',a20,12x,''['',a,'']'','''//rstring//''')'
c--- f97 integer:: format
      f97='('''//lstring//''',i20,12x,''['',a,'']'','''//rstring//''')'
c--- f98 logical:: format
      f98='('''//lstring//''',L20,12x,''['',a,'']'','''//rstring//''')'
c--- f99 floating point format
      f99='('''//lstring//''',f20.4,12x,''['',a,'']'','''
     & //rstring//''')'

      writeall=.false.
      if (tag == 'WRITEALL') writeall=.true.

      if (writeall) then
      write(unitno,*) lstring//' Run corresponds to this input file)'
      write(unitno,*)
      write(unitno,*)
     & lstring//' [Flags to specify the mode in which MCFM is run] )'
      endif

      if ((tag == 'nevtrequested') .or. (writeall)) then
      write(unitno,fmt=f97) nevtrequested,'nevtrequested'
      endif
      if ((tag == 'creatent') .or. (writeall)) then
      write(unitno,fmt=f98) creatent,'creatent'
      endif
!      if ((tag == 'skipnt') .or. (writeall)) then
!      write(unitno,fmt=f98) skipnt,'skipnt'
!      endif
      if ((tag == 'dswhisto') .or. (writeall)) then
      write(unitno,fmt=f98) dswhisto,'dswhisto'
      endif
c      if ((tag == 'creategrid') .or. (writeall)) then
c      write(unitno,fmt=f98) creategrid,'creategrid'
c      endif
      if ((tag == 'writerefs') .or. (writeall)) then
      write(unitno,fmt=f98) writerefs,'writerefs'
      endif
      if ((tag == 'writetop') .or. (writeall)) then
      write(unitno,fmt=f98) writetop,'writetop'
      endif
      if ((tag == 'writedat') .or. (writeall)) then
      write(unitno,fmt=f98) writedat,'writedat'
      endif
      if ((tag == 'writegnu') .or. (writeall)) then
      write(unitno,fmt=f98) writegnu,'writegnu'
      endif
      if ((tag == 'writeroot') .or. (writeall)) then
      write(unitno,fmt=f98) writeroot,'writeroot'
      endif
      if ((tag == 'writepwg') .or. (writeall)) then
      write(unitno,fmt=f98) writepwg,'writepwg'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*) lstring//
     & ' [General options to specify the process and execution] )'
      endif
      if (vdecayid) then
        if ((tag == 'nproc') .or. (writeall)) then
        write(unitno,fmt=f94) nproc,v34id,v56id,'nproc'
        endif
      else
        if ((tag == 'nproc') .or. (writeall)) then
        write(unitno,fmt=f97) nproc,'nproc'
        endif
      endif
      if ((tag == 'part') .or. (writeall)) then
      write(unitno,fmt=f96) trim(kpartstring(kpart)),'part'
      endif
      if ((tag == 'runstring') .or. (writeall)) then
      write(unitno,fmt=f96) runstring,'runstring'
      endif
      if ((tag == 'sqrts') .or. (writeall)) then
      write(unitno,fmt=f99) sqrts,'sqrts'
      endif
      if ((tag == 'ih1') .or. (writeall)) then
      write(unitno,fmt=f97) ih1,'ih1'
      endif
      if ((tag == 'ih2') .or. (writeall)) then
      write(unitno,fmt=f97) ih2,'ih2'
      endif
      if ((tag == 'hmass') .or. (writeall)) then
      write(unitno,fmt=f99) hmass,'hmass'
      endif

c--- catch special scale choices for stop+b process
      if ((nproc >= 231) .and. (nproc <= 240)) then
         if ((tag == 'renscale_L') .or. (writeall)) then
         write(unitno,fmt=f99) renscale_L,'renscale_L'
         endif
         if ((tag == 'facscale_L') .or. (writeall)) then
         write(unitno,fmt=f99) facscale_L,'facscale_L'
         endif
         if ((tag == 'renscale_H') .or. (writeall)) then
         write(unitno,fmt=f99) renscale_H,'renscale_H'
         endif
         if ((tag == 'facscale_H') .or. (writeall)) then
         write(unitno,fmt=f99) facscale_H,'facscale_H'
         endif
      else
         if ((tag == 'scale') .or. (writeall)) then
         write(unitno,fmt=f99) initscale,'scale'
         endif
         if ((tag == 'facscale') .or. (writeall)) then
         write(unitno,fmt=f99) initfacscale,'facscale'
         endif
      endif


      if ((tag == 'dynamicscale') .or. writeall) then
      write(unitno,fmt=f96) dynstring,'dynamicscale'
      endif
      if ((tag == 'zerowidth') .or. (writeall)) then
      write(unitno,fmt=f98) zerowidth,'zerowidth'
      endif
      if ((tag == 'removebr') .or. (writeall)) then
      write(unitno,fmt=f98) removebr,'removebr'
      endif
      if ((tag == 'itmx1') .or. (writeall)) then
      write(unitno,fmt=f97) itmx1,'itmx1'
      endif
      if ((tag == 'ncall1') .or. (writeall)) then
      write(unitno,fmt=f97) ncall1,'ncall1'
      endif
      if ((tag == 'itmx2') .or. (writeall)) then
      write(unitno,fmt=f97) itmx2,'itmx2'
      endif
      if ((tag == 'ncall2') .or. (writeall)) then
      write(unitno,fmt=f97) ncall2,'ncall2'
      endif
      if ((tag == 'taucut') .or. (writeall)) then
      write(unitno,fmt=f93) taucut,'taucut'
      endif
      if ((tag == 'ij') .or. (writeall)) then
      write(unitno,fmt=f97) origij,'ij'
      endif
      if ((tag == 'dryrun') .or. (writeall)) then
      write(unitno,fmt=f98) dryrun,'dryrun'
      endif
      if ((tag == 'Qflag') .or. (writeall)) then
      write(unitno,fmt=f98) Qflag,'Qflag'
      endif
      if ((tag == 'Gflag') .or. (writeall)) then
      write(unitno,fmt=f98) Gflag,'Gflag'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     & lstring//' [Heavy quark masses] )'
      endif
      if ((tag == 'top mass') .or. (writeall)) then
      write(unitno,fmt=f99) mt,'top mass'
      endif
      if ((tag == 'bottom mass') .or. (writeall)) then
      write(unitno,fmt=f99) mb,'bottom mass'
      endif
      if ((tag == 'charm mass') .or. (writeall)) then
      write(unitno,fmt=f99) mc,'charm mass'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     & lstring//' [Pdf selection] )'
      endif
      if ((tag == 'pdlabel') .or. (writeall)) then
      write(unitno,fmt=f96) pdlabel,'pdlabel'
      endif
      if ((tag == 'NGROUP') .or. (writeall)) then
      write(unitno,fmt=f97) NGROUP,'NGROUP'
      endif
      if ((tag == 'NSET') .or. (writeall)) then
      write(unitno,fmt=f97) NSET,'NSET'
      endif
      if ((tag == 'LHAPDF group') .or. (writeall)) then
      write(unitno,fmt=f96) PDFname,'LHAPDF group'
      endif
      if ((tag == 'LHAPDF set') .or. (writeall)) then
      write(unitno,fmt=f97) PDFmember,'LHAPDF set'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     & lstring//' [Jet definition and event cuts] )'
      endif
      if ((tag == 'm34min') .or. (writeall)) then
      write(unitno,fmt=f99) sqrt(wsqmin),'m34min'
      endif
      if ((tag == 'm34max') .or. (writeall)) then
      write(unitno,fmt=f99) sqrt(wsqmax),'m34max'
      endif
      if ((tag == 'm56min') .or. (writeall)) then
      write(unitno,fmt=f99) sqrt(bbsqmin),'m56min'
      endif
      if ((tag == 'm56max') .or. (writeall)) then
      write(unitno,fmt=f99) sqrt(bbsqmax),'m56max'
      endif
      if ((tag == 'm3456min') .or. (writeall)) then
      write(unitno,fmt=f99) m3456min,'m3456min'
      endif
      if ((tag == 'm3456max') .or. (writeall)) then
      write(unitno,fmt=f99) m3456max,'m3456max'
      endif
      if ((tag == 'inclusive') .or. (writeall)) then
      write(unitno,fmt=f98) inclusive,'inclusive'
      endif
      if ((tag == 'algorithm') .or. (writeall)) then
      write(unitno,fmt=f96) algorithm,'algorithm'
      endif
      if ((tag == 'ptjetmin') .or. (writeall)) then
      write(unitno,fmt=f99) ptjetmin,'ptjetmin'
      endif
      if ((tag == 'etajetmin') .or. (writeall)) then
      write(unitno,fmt=f99) etajetmin,'etajetmin'
      endif
      if ((tag == 'etajetmax') .or. (writeall)) then
      write(unitno,fmt=f99) etajetmax,'etajetmax'
      endif
      if ((tag == 'Rcut') .or. (writeall)) then
      write(unitno,fmt=f99) Rcut,'Rcut'
      endif
      if ((tag == 'makecuts') .or. (writeall)) then
      write(unitno,fmt=f98) makecuts,'makecuts'
      endif
      if ((tag == 'leptpt') .or. (writeall)) then
      write(unitno,fmt=f99) leptpt,'leptpt'
      endif
      if ((tag == 'leptrap') .or. (writeall)) then
      write(unitno,fmt=f99) leptrap,'leptrap'
      endif
      if ((tag == 'leptveto') .or. (writeall)) then
      write(unitno,fmt=f95) leptveto1min,leptveto1max,'leptveto'
      endif
      if ((tag == 'misspt') .or. (writeall)) then
      write(unitno,fmt=f99) misspt,'misspt'
      endif
      if ((tag == 'leptpt2') .or. (writeall)) then
      write(unitno,fmt=f99) leptpt2,'leptpt2'
      endif
      if ((tag == 'leptrap2') .or. (writeall)) then
      write(unitno,fmt=f99) leptrap2,'leptrap2'
      endif
      if ((tag == 'leptveto2') .or. (writeall)) then
      write(unitno,fmt=f95) leptveto2min,leptveto2max,'leptveto2'
      endif
      if ((tag == 'mtrans34cut') .or. (writeall)) then
      write(unitno,fmt=f99) mtrans34cut,'mtrans34cut'
      endif
      if ((tag == 'Rjlmin') .or. (writeall)) then
      write(unitno,fmt=f99) Rjlmin,'Rjlmin'
      endif
      if ((tag == 'Rllmin') .or. (writeall)) then
      write(unitno,fmt=f99) Rllmin,'Rllmin'
      endif
      if ((tag == 'delyjjmin') .or. (writeall)) then
      write(unitno,fmt=f99) delyjjmin,'delyjjmin'
      endif
      if ((tag == 'jetsopphem') .or. (writeall)) then
      write(unitno,fmt=f98) jetsopphem,'jetsopphem'
      endif
      if ((tag == 'lbjscheme') .or. (writeall)) then
      write(unitno,fmt=f97) lbjscheme,'lbjscheme'
      endif
      if ((tag == 'ptbjetmin') .or. (writeall)) then
      write(unitno,fmt=f99) ptbjetmin,'ptbjetmin'
      endif
      if ((tag == 'etabjetmax') .or. (writeall)) then
      write(unitno,fmt=f99) etabjetmax,'etabjetmax'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     & lstring//' [Settings for photon processes] )'
      endif
      if ((tag == 'frag') .or. (writeall)) then
      write(unitno,fmt=f98) frag,'frag'
      endif
      if ((tag == 'fragset') .or. (writeall)) then
      write(unitno,fmt=f96) fragset,'fragset'
      endif
      if ((tag == 'frag_scale') .or. (writeall)) then
      write(unitno,fmt=f99) frag_scale,'frag_scale'
      endif
      if ((tag == 'gammpt') .or. (writeall)) then
      write(unitno,fmt=f99) gammpt,'gammpt'
      endif
      if ((tag == 'gammrap') .or. (writeall)) then
      write(unitno,fmt=f99) gammrap,'gammrap'
      endif
      if ((tag == 'gammpt2') .or. (writeall)) then
      write(unitno,fmt=f99) gammpt2,'gammpt2'
      endif
      if ((tag == 'gammpt3') .or. (writeall)) then
      write(unitno,fmt=f99) gammpt2,'gammpt3'
      endif
      if ((tag == 'Rgalmin') .or. (writeall)) then
      write(unitno,fmt=f99) Rgalmin,'Rgalmin'
      endif
      if ((tag == 'Rgagamin') .or. (writeall)) then
      write(unitno,fmt=f99) Rgagamin,'Rgagamin'
      endif
      if ((tag == 'Rgajetmin') .or. (writeall)) then
      write(unitno,fmt=f99) Rgajetmin,'Rgajetmin'
      endif
      if ((tag == 'cone_ang') .or. (writeall)) then
      write(unitno,fmt=f99) cone_ang,'cone_ang'
      endif
      if ((tag == 'epsilon_h') .or. (writeall)) then
      write(unitno,fmt=f99) epsilon_h,'epsilon_h'
      endif
      if ((tag == 'n_pow') .or. (writeall)) then
      write(unitno,fmt=f99) n_pow,'n_pow'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     & lstring//' [Anomalous couplings of the W and Z] )'
      endif
      if ((tag == 'delg1_z') .or. (writeall)) then
      write(unitno,fmt=f99) delg1_z,'delg1_z'
      endif
      if ((tag == 'delk_z') .or. (writeall)) then
      write(unitno,fmt=f99) delk_z,'delk_z'
      endif
      if ((tag == 'delk_g') .or. (writeall)) then
      write(unitno,fmt=f99) delk_g,'delk_g'
      endif
      if ((tag == 'lambda_z') .or. (writeall)) then
      write(unitno,fmt=f99) lambda_z,'lambda_z'
      endif
      if ((tag == 'lambda_g') .or. (writeall)) then
      write(unitno,fmt=f99) lambda_g,'lambda_g'
      endif
      if ((tag == 'h1Z') .or. (writeall)) then
      write(unitno,fmt=f99) h1Z,'h1Z'
      endif
      if ((tag == 'h1gam') .or. (writeall)) then
      write(unitno,fmt=f99) h1gam,'h1gam'
      endif
      if ((tag == 'h2Z') .or. (writeall)) then
      write(unitno,fmt=f99) h2Z,'h2Z'
      endif
      if ((tag == 'h2gam') .or. (writeall)) then
      write(unitno,fmt=f99) h2gam,'h2gam'
      endif
      if ((tag == 'h3Z') .or. (writeall)) then
      write(unitno,fmt=f99) h3Z,'h3Z'
      endif
      if ((tag == 'h3gam') .or. (writeall)) then
      write(unitno,fmt=f99) h3gam,'h3gam'
      endif
      if ((tag == 'h4Z') .or. (writeall)) then
      write(unitno,fmt=f99) h4Z,'h4Z'
      endif
      if ((tag == 'h4gam') .or. (writeall)) then
      write(unitno,fmt=f99) h4gam,'h4gam'
      endif
      if ((tag == 'tevscale') .or. (writeall)) then
      write(unitno,fmt=f99) tevscale,'tevscale'
      endif
      if ((tag == 'cttH') .or. (writeall)) then
      write(unitno,fmt=f99) cttH,'cttH'
      endif
      if ((tag == 'cWWH') .or. (writeall)) then
      write(unitno,fmt=f99) cWWH,'cWWH'
      endif
      if ((tag == 'hwidth_ratio') .or. (writeall)) then
      write(unitno,fmt=f99) hwidth_ratio,'Gamma_H/Gamma_H(SM)'
      endif


      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     & lstring//' [How to resume/save a run] )'
      endif
      if ((tag == 'readin') .or. (writeall)) then
      write(unitno,fmt=f98) readin,'readin'
      endif
      if ((tag == 'writeout') .or. (writeall)) then
      write(unitno,fmt=f98) writeout,'writeout'
      endif
      if ((tag == 'ingridfile') .or. (writeall)) then
      write(unitno,fmt=f96) ingridfile,'ingridfile'
      endif
      if ((tag == 'outgridfile') .or. (writeall)) then
      write(unitno,fmt=f96) outgridfile,'outgridfile'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*) lstring//
     & ' [Technical parameters that should not normally be changed]'
      endif
      if ((tag == 'debug') .or. (writeall)) then
      write(unitno,fmt=f98) debug,'debug'
      endif
      if ((tag == 'verbose') .or. (writeall)) then
      write(unitno,fmt=f98) verbose,'verbose'
      endif
      if ((tag == 'new_pspace') .or. (writeall)) then
      write(unitno,fmt=f98) new_pspace,'new_pspace'
      endif
      if ((tag == 'virtonly') .or. (writeall)) then
      write(unitno,fmt=f98) virtonly,'virtonly'
      endif
      if ((tag == 'realonly') .or. (writeall)) then
      write(unitno,fmt=f98) realonly,'realonly'
      endif
      if ((tag == 'spira') .or. (writeall)) then
      write(unitno,fmt=f98) spira,'spira'
      endif
      if ((tag == 'noglue') .or. (writeall)) then
      write(unitno,fmt=f98) noglue,'noglue'
      endif
      if ((tag == 'ggonly') .or. (writeall)) then
      write(unitno,fmt=f98) ggonly,'ggonly'
      endif
      if ((tag == 'gqonly') .or. (writeall)) then
      write(unitno,fmt=f98) gqonly,'gqonly'
      endif
      if ((tag == 'omitgg') .or. (writeall)) then
      write(unitno,fmt=f98) omitgg,'omitgg'
      endif
      if ((tag == 'vanillafiles') .or. (writeall)) then
      write(unitno,fmt=f98) vanillafiles,'vanillafiles'
      endif
!      if ((tag == 'nmin') .or. (writeall)) then
!      write(unitno,fmt=f97) nmin,'nmin'
!      endif
!      if ((tag == 'nmax') .or. (writeall)) then
!      write(unitno,fmt=f97) nmax,'nmax'
!      endif
      if ((tag == 'clustering') .or. (writeall)) then
      write(unitno,fmt=f98) clustering,'clustering'
      endif
      if ((tag == 'realwt') .or. (writeall)) then
      write(unitno,fmt=f98) realwt,'realwt'
      endif
      if ((tag == 'colourchoice') .or. (writeall)) then
      write(unitno,fmt=f97) colourchoice,'colourchoice'
      endif
      if ((tag == 'rtsmin') .or. (writeall)) then
      write(unitno,fmt=f93) rtsmin,'rtsmin'
      endif
      if ((tag == 'cutoff') .or. (writeall)) then
      write(unitno,fmt=f93) cutoff,'cutoff'
      endif
      if ((tag == 'aii') .or. (writeall)) then
      write(unitno,fmt=f99) aii,'aii'
      endif
      if ((tag == 'aif') .or. (writeall)) then
      write(unitno,fmt=f99) aif,'aif'
      endif
      if ((tag == 'afi') .or. (writeall)) then
      write(unitno,fmt=f99) afi,'afi'
      endif
      if ((tag == 'aff') .or. (writeall)) then
      write(unitno,fmt=f99) aff,'aff'
      endif
      if ((tag == 'bfi') .or. (writeall)) then
      write(unitno,fmt=f99) bfi,'bfi'
      endif
      if ((tag == 'bff') .or. (writeall)) then
      write(unitno,fmt=f99) bff,'bff'
      endif

      if (writeall) then
      write(unitno,*)
      endif

      return

c--- 96 character format
c   96 format(' (',a20,12x,'[',a,']',' )')
c--- 97 integer:: format
c   97 format(' (',i20,12x,'[',a,']',' )')
c--- 98 logical:: format
c   98 format(' (',L20,12x,'[',a,']',' )')
c--- 99 floating point format
c   99 format(' (',f20.4,12x,'[',a,']',' )')

      end

