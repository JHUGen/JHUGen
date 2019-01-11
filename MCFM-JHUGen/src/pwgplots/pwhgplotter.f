      subroutine pwhgplotter(p,pjet,wt0,nd)
      implicit none
      include 'hepevt.h'
      include 'masses.f'
      include 'part.f'
      include 'nproc.f'
      double precision p(12,4),pjet(12,4),ph(4),wt0,pt(1:4),
     1     d1t,d2t,d1a,d2a
      integer switch,nd,ond,j
      logical first
      data first/.true./
      save first,ond
      double precision pjetcom(12,4),wt
      common/pjetcom/pjetcom

c Translate from femto- (MCFM) to pico- (POWHEG)
      wt=wt0/1000

c------------ SET-UP NUMBER OF PARTICLES TO APPEAR IN EVENT RECORD

      if(first) then
c This is trick to communicate the number of jets
c to the histogram booking section
c JC/RKE - what's the purpose of setting these twice?
         if    (nproc.eq.31) then
            nhep=4
            call init_hist_Z
         elseif(nproc.eq.44) then
            nhep=6
            call init_hist_KN
         elseif(nproc.eq.203) then
            nhep=6
            call init_hist_KN
         elseif(nproc.ge.141.and.nproc.le.151) then
            nhep=12
            call init_hist_KN
         elseif(nproc.ge.171.and.nproc.le.177) then
            nhep=8
            call init_hist_ST_sch_dk
         elseif(nproc.ge.231.and.nproc.le.239) then
            nhep=9
            call init_hist_ST_tch_dk
         elseif(nproc.ge.181.and.nproc.le.187) then
            nhep=10
            call init_hist_ST_wt_dk
         else
            stop 'process not implemented in pwhgplotter.f' 
         endif
         ond=0
         first = .false. 
      else
c put partons in the hephep interface
         if    (nproc.eq.31) then
            nhep=4
         elseif(nproc.eq.44) then
            nhep=6
         elseif(nproc.eq.203) then
            nhep=5
         elseif(nproc.ge.141.and.nproc.le.151) then
            nhep=12
         elseif(nproc.ge.171.and.nproc.le.177) then
            nhep=8
         elseif(nproc.ge.231.and.nproc.le.239) then
            nhep=9
         elseif(nproc.ge.181.and.nproc.le.187) then
            nhep=10
         else
            stop 'process not implemented in pwhgplotter.f' 
         endif

c--- additional parton present in real corrections
         if(part.eq.'real'.and.nd.eq.0) then
            nhep=nhep+1
         endif
         
c------------ EVENT RECORD INITIALIZATION 

c----spuriously set everyone to a gluon?
         do j=1,nhep
            idhep(j)=21
            isthep(j)=1
         enddo
         isthep(1)=-1
         isthep(2)=-1

c--- JC/RKE: should leptons be properly defined for all the processes below?

c------------ PROCESS-SPECIFIC FILLING OF EVENT RECORD

c------ Z
         if(nproc.eq.31) then
            idhep(3)=11
            idhep(4)=-11
            do j=1,nhep
               phep(1:4,j)=p(j,:)
            enddo
         endif
         
c------ W/Z + jets
         if(nproc.eq.11.or.nproc.eq.22.or.nproc.eq.44) then
            idhep(3)=11
            idhep(4)=-11
            do j=1,nhep
               phep(1:4,j)=p(j,:)
            enddo
            do j=1,12
               pjetcom(j,1:4)=pjet(j,1:4)
            enddo
         endif
         
c------ H + jet
         if(nproc.eq.203) then
            idhep(3)=25
c            ph=p(3,:)+p(4,:)
c            write(*,*) 
c     1           sqrt(abs(ph(4)**2-ph(1)**2-ph(2)**2-ph(3)**2))
c we want the undecayed Higgs
            do j=1,2
               phep(1:4,j)=p(j,:)
            enddo
            phep(1:4,3)=p(3,:)+p(4,:)
            do j=4,nhep
               phep(1:4,j)=p(j+1,:)
            enddo
         endif

c------ t-tbar            
         if(nproc.ge.141.and.nproc.le.151) then
            phep(1:4,3)=p(3,:)+p(4,:)+p(5,:)
            phep(1:4,4)=p(6,:)+p(7,:)+p(8,:)
            phep(1:4,5)=p(3,:)+p(4,:)
            phep(1:4,6)=p(7,:)+p(8,:)
            phep(1:4,7)=p(4,:)
            phep(1:4,8)=p(3,:)
            phep(1:4,9)=p(7,:)
            phep(1:4,10)=p(8,:)
            phep(1:4,11)=p(5,:)
            phep(1:4,12)=p(6,:)
            if(nhep.eq.13) then
                phep(1:4,13)=p(9,:)
c is 3 a full top?
                pt=phep(1:4,3)
                d1t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
                pt=pt+phep(1:4,13)
                d2t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
                pt=phep(1:4,4)
                d1a=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
                pt=pt+phep(1:4,13)
                d2a=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
c                write(*,*) d1t,d2t,d1a,d2a
                if(d2t.lt.d1t.and.d2a.ge.d2t) then
c                   write(*,*) ' t with radiation'
                   phep(1:4,3)=phep(1:4,3)+phep(1:4,13)
                elseif(d2a.lt.d1a.and.d2t.ge.d2a) then
c                   write(*,*) ' at with radiation'
                   phep(1:4,4)=phep(1:4,4)+phep(1:4,13)
                endif
c                read(*,*)
            endif
            idhep(3)=6
            idhep(4)=-6
            idhep(5)=24
            idhep(6)=-24
            idhep(11)=5
            idhep(12)=-5
         endif

c------ s-channel single top            
         if(nproc.ge.171.and.nproc.le.177) then
            phep(1:4,3)=p(3,:)+p(4,:)+p(5,:)
            phep(1:4,4)=p(3,:)+p(4,:)
            phep(1:4,5)=p(3,:)
            phep(1:4,6)=p(4,:)
            phep(1:4,7)=p(5,:)
            phep(1:4,8)=p(6,:)
            if(nhep.eq.9) then
                phep(1:4,9)=p(7,:)
c is 3 a full top?
                pt=phep(1:4,3)
                d1t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
                pt=pt+phep(1:4,9)
                d2t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
c                write(*,*) d1t,d2t
                if(d2t.lt.d1t) then
c                   write(*,*) ' t with radiation'
                   phep(1:4,3)=phep(1:4,3)+phep(1:4,9)
                endif
c                read(*,*)
            endif
            idhep(3)=6
            idhep(4)=24
            idhep(7)=5
            idhep(8)=-5
         endif
         
c------ t-channel single top            
         if(nproc.ge.231.and.nproc.le.239) then
            phep(1:4,3)=p(3,:)+p(4,:)+p(5,:)
            phep(1:4,4)=p(3,:)+p(4,:)
            phep(1:4,5)=p(3,:)
            phep(1:4,6)=p(4,:)
            phep(1:4,7)=p(5,:)
            phep(1:4,8)=p(6,:)
            phep(1:4,9)=p(7,:)
            if(nhep.eq.10) then
                phep(1:4,10)=p(8,:)
c is 3 a full top?
                pt=phep(1:4,3)
                d1t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
                pt=pt+phep(1:4,10)
                d2t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
c                write(*,*) d1t,d2t
                if(d2t.lt.d1t) then
c                   write(*,*) ' t with radiation'
                   phep(1:4,3)=phep(1:4,3)+phep(1:4,10)
                endif
c                read(*,*)
            endif
            idhep(3)=6
            idhep(4)=24
            idhep(7)=5
            idhep(8)=-5
        endif
        
c------ W-t-single top            
         if(nproc.ge.181.and.nproc.le.187) then
            phep(1:4,3)=p(5,:)+p(6,:)+p(7,:)
            phep(1:4,4)=p(5,:)+p(6,:)
            phep(1:4,5)=p(3,:)+p(4,:)
            phep(1:4,6)=p(3,:)
            phep(1:4,7)=p(4,:)
            phep(1:4,8)=p(5,:)
            phep(1:4,9)=p(6,:)
            phep(1:4,10)=p(7,:)
            if(nhep.eq.11) then
                phep(1:4,11)=p(8,:)
c is 3 a full top?
                pt=phep(1:4,3)
                d1t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
                pt=pt+phep(1:4,11)
                d2t=abs(pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2-mt**2)/mt**2
c                write(*,*) d1t,d2t
                if(d2t.lt.d1t) then
c                   write(*,*) ' t with radiation'
                   phep(1:4,3)=phep(1:4,3)+phep(1:4,11)
                endif
c                read(*,*)
            endif
            idhep(3)=6
            idhep(4)=24
            idhep(5)=-24
            idhep(10)=5
         endif
         
c------------ FILL HISTOGRAMS

         if(part.ne.'real') then
            call analysis_wrap(wt)
            call pwhgaccumup
         else
            if(nd.eq.0.and.ond.ne.0) then
               call pwhgaccumup
               call analysis_wrap(wt)
            else
               call analysis_wrap(wt)
            endif
            ond=nd
         endif
      endif

      end

   
      subroutine analysis_wrap(wt)
      implicit none
      double precision wt
      integer nproc
      common/nproc/nproc
      
         if(nproc.eq.44) then
            write(6,*) 'No analysis routine written for nproc=',nproc
            stop
         elseif(nproc.eq.203) then
            write(6,*) 'No analysis routine written for nproc=',nproc
            stop
         elseif(nproc.ge.31) then
            call analysis_Z(wt)
         elseif(nproc.ge.141.and.nproc.le.151) then
            call analysis_KN(wt)
         elseif(nproc.ge.171.and.nproc.le.177) then
            call analysis_ST_sch_dk(wt)
         elseif(nproc.ge.231.and.nproc.le.239) then
            call analysis_ST_tch_dk(wt)
         elseif(nproc.ge.181.and.nproc.le.186) then
            call analysis_ST_wt_dk(wt)
         else
            write(6,*) 'No analysis routine written for nproc=',nproc
            stop
         endif
      end   



      subroutine pwhghistofin(itno,itmx)
      implicit none
      include 'part.f'
      integer itno,itmx
      character *255 runname,outfile 
      integer iun, iostat
      common/runname/runname
      parameter (iun=99)

      if(part.eq.'real') then
c     in this case a call to accumup is missin
         call pwhgaccumup
      endif
c finalize histograms
      call pwhgsetout
c mcfm delivered weights divided by n; should now multiply
c by n to get proper normalization. Furthermore, if itno#0
c we should divide by itno; otherwise we should divide by itmx
      if(itno.gt.0) then
         call pwhgfixmcfm(itno)
      else
         call pwhgfixmcfm(itmx)
      endif
c print histograms
      
      outfile = trim(runname)//'pwg.gdat'
      open(unit=iun,file=outfile,status='unknown',iostat=iostat) 
      call pwhgtopout
      close(iun) 
      
      end

      subroutine pwhgfixmcfm(it)
c analysis, leaving the yhistarr1 and errhistarr1 unchanged.
      implicit none
      integer it
      include 'pwhg_bookhist-new.h'
      integer j,k
      real *8 xxx,sum,sumsq
      do j=1,jhist
         xxx=ient1(j)
         xxx=xxx/it
         do k=0,nbins(j)+1
            yhistarr2(k,j)=yhistarr2(k,j)*xxx
            errhistarr2(k,j)=errhistarr2(k,j)*xxx
         enddo
      enddo
      end
