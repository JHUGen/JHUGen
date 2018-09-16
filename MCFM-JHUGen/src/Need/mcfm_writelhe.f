c--- Routine to write a generic LHE file for a given MCFM process
      subroutine mcfm_writelhe(pin,xmsq,xfac)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'plabel.f'
      include 'eventbuffer.f'
      include 'facscale.f'
      include 'qcdcouple.f'
      include 'maxwt.f'
      include 'montecarlorpp.f'
      include 'hepeup.f'
      include 'heprup.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      integer:: ilomomenta,plabeltoPDG,ip_parent(10),id_parent(10),ic,
     & i,j,k,jj,kk,j1,j2,j3,j4,ip,nu,ilen
      real(dp):: p(mxpart,4),xfac,p_parent(10,4),pin(mxpart,4),
     & xmsq(-nf:nf,-nf:nf),xx,mm
      character*255 runname,outputstring
      common/ilomomenta/ilomomenta
      common/runname/runname
      logical, save :: first=.true.

c--- work out flavour to use for initial state
c--- (randomly, based on weights passed in xmsq)
      call mcfm_getflavour(xmsq,jj,kk)
      if (jj == 0) then
        idup(1)=21
      else
        idup(1)=jj
      endif
      if (kk == 0) then
        idup(2)=21
      else
        idup(2)=kk
      endif
      istup(1)=-1
      istup(2)=-1

c--- fill simple particle labels for final state
      do i=3,ilomomenta
        idup(i)=plabeltoPDG(plabel(i))
        istup(i)=1
      enddo

c--- copy array pin to new array p
      p(:,:)=pin(:,:)

c--- get parent combinations: at the end of this loop there will be ip of them
      call getparents(ip_parent,id_parent)
      ip=1
      do while (ip_parent(ip) > 0)
        j=ip_parent(ip)
        idup(ilomomenta+ip)=id_parent(ip)
        if     (j < 100) then
          j1=j/10
          j2=mod(j,10)
          if (j2 == 0) j2=10   ! special code: 0 -> 10
          mothup(1,ilomomenta+ip)=1
          mothup(2,ilomomenta+ip)=2
          istup(ilomomenta+ip)=+2
          mothup(:,j1)=ilomomenta+ip
          mothup(:,j2)=ilomomenta+ip
          call fixmasses(p,id_parent(ip),idup(j1),idup(j2),j1,j2)
          p_parent(ip,:)=p(j1,:)+p(j2,:)
        elseif (j < 1000) then
c--- three-particle plots
          j2=(j-j1*100)/10
          j3=mod(j,10)
          if (j3 == 0) j3=10   ! special code: 0 -> 10
          write(6,*) 'Unfinished mcfm_writelhe: j=',j
          stop
        elseif (j < 10000) then
          j1=j/1000
          j2=(j-j1*1000)/100
          j3=(j-j1*1000-j2*100)/10
          j4=mod(j,10)
          if (j4 == 0) j4=10   ! special code: 0 -> 10
          write(6,*) 'Unfinished mcfm_writelhe: j=',j
          stop
        else
          write(6,*) 'Error in ip_parent: ',j
          stop
        endif
        ip=ip+1
      enddo
      ip=ip-1 ! over-counted

c--- now handle color
      if (nqcdjets == 0) then
c---   for no jets in the final state, straightforward
        if     (jj > 0) then
          icolup(1,1)=501
          icolup(2,1)=0
          icolup(1,2)=0
          icolup(2,2)=501
        elseif (jj < 0) then
          icolup(1,1)=0
          icolup(2,1)=501
          icolup(1,2)=501
          icolup(2,2)=0
        else
          icolup(1,1)=501
          icolup(2,1)=502
          icolup(1,2)=502
          icolup(2,2)=501
        endif
      else
        write(6,*) 'LHE output not yet available for nqcdjets=',nqcdjets
        stop
      endif

c--- number of entries to write
      nup=ilomomenta+ip

c--- fill momenta
      do i=1,nup
        do nu=1,4
          if     (i <= 2) then
            xx=-p(i,nu)
          elseif (i <= ilomomenta) then
            xx=p(i,nu)
          else
            xx=p_parent(i-ilomomenta,nu)
          endif
          pup(nu,i)=xx
        enddo
        if     (i <= 2) then
          mm=zip
        elseif (i <= ilomomenta) then
          mm=sqrt(max(zip,p(i,4)**2-p(i,1)**2-p(i,2)**2-p(i,3)**2))
          if (mm < 1.e-12_dp) mm=zip
        else
          mm=sqrt(max(zip,
     &     +p_parent(i-ilomomenta,4)**2-p_parent(i-ilomomenta,1)**2
     &     -p_parent(i-ilomomenta,2)**2-p_parent(i-ilomomenta,3)**2))
        endif
        pup(5,i)=mm
      enddo

c--- Event weight must be multiplied by xfac to account for
c--- events with negative weight or events with weights that exceed wtmax
      xwgtup=xfac
      xmaxup(1)=wtmax

c--- Miscellaneous info
      idprup=10000+nproc
      scalup=facscale
      aqedup=-one
      aqcdup=as

c--- junk entries
      vtimup(:)=zip
      spinup(:)=9._dp

c--- on the first pass, open unit 84 with the right file name and write header
      if(first) then
         first=.false.
         ilen=len_trim(runname)
         outputstring=runname(1:ilen)//'.lhe'
         open(unit=84,file=outputstring,status='unknown')
         call init_lhe_events(84)
      endif

c--- increment event counter and write out entry to unit 84
      numstored=numstored+1
      call lhefwritev(84)

      return
      end


c--- Routine to specify parent particles
      subroutine getparents(ip_parent,id_parent)
      implicit none
      include 'types.f'

      include 'nwz.f'
      include 'kprocess.f'
      include 'montecarlorpp.f'
      integer:: j,ip_parent(10),id_parent(10)

      ip_parent(:)=0

      if ( (kcase==kHZZ_4l) .or. (kcase==kHZZ_tb)
     & .or.(kcase==kHZZint) .or. (kcase==kggZZ4l)
     & .or.(kcase==kHZZHpi) .or. (kcase==kZZlept)
     & .or.(kcase==kggZZbx)) then
        ip_parent= (/ 34, 56, (0,j=1,8) /)
        id_parent= (/ z0_pdg, z0_pdg, (0,j=1,8) /)
      endif

      if ( (kcase==kHWW_4l) .or. (kcase==kHWW_tb)
     & .or.(kcase==kHWWint) .or. (kcase==kggWW4l)
     & .or.(kcase==kHWWHpi) .or. (kcase==kWWqqbr)
     & .or.(kcase==kggWWbx)) then
        ip_parent= (/ 34, 56, (0,j=1,8) /)
        id_parent= (/ wp_pdg, wm_pdg, (0,j=1,8) /)
      endif

      if (kcase==kW_only) then
        ip_parent= (/ 34, (0,j=1,9) /)
        id_parent= (/ wp_pdg*nwz, (0,j=1,9) /)
      endif

      if (kcase==kZ_only) then
        ip_parent= (/ 34, (0,j=1,9) /)
        id_parent= (/ z0_pdg, (0,j=1,9) /)
      endif

      if (kcase==kggfus0) then
        ip_parent= (/ 34, (0,j=1,9) /)
        id_parent= (/ h_pdg, (0,j=1,9) /)
      endif

      if ( (kcase==kWWqqbr) .or. (kcase==kWWnpol)
     & .or.(kcase==kWWqqdk)) then
        ip_parent= (/ 34, 56, (0,j=1,8) /)
        id_parent= (/ wp_pdg, wm_pdg, (0,j=1,8) /)
      endif

      if (kcase==kWZbbar) then
        ip_parent= (/ 34, 56, (0,j=1,8) /)
        id_parent= (/ wp_pdg*nwz, z0_pdg, (0,j=1,8) /)
      endif

      return
      end


c--- Routine to transform massless momenta to massive momenta,
c--- satisfying constraint that parent particle momentum is unchanged
c--- This is implemented by using "Rodrigo" transformation
c---
c---  p -> momentum array that will be changed
c---  idparent -> PDG code of parent particle
c---  id1,id2 -> PDG codes of daughter particles
c---  i1,id2 -> PDG codes of daughter particles
c---  j1,j2 -> entries in momentum array for daughters
      subroutine fixmasses(p,idparent,id1,id2,j1,j2)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'montecarlorpp.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),dot,m,beta,p1(4),p2(4)
      integer:: idparent,id1,id2,j1,j2,imass

      if     (idparent == z0_pdg) then
c--- Z parent
        if     ((abs(id1)==el_pdg) .and. (abs(id2)==el_pdg)) then
          m=mel   ! electrons
        elseif ((abs(id1)==ml_pdg) .and. (abs(id2)==ml_pdg)) then
          m=mmu   ! muons
        elseif ((abs(id1)==tl_pdg) .and. (abs(id2)==tl_pdg)) then
          m=mtau  ! taus
        elseif ((abs(id1)==bq_pdg) .and. (abs(id2)==bq_pdg)) then
          m=mb    ! b-quarks
        elseif ((abs(id1)==nel_pdg) .and. (abs(id2)==nel_pdg))then
          return  ! massless neutrinos -> nothing to do
        elseif ((abs(id1)==nml_pdg) .and. (abs(id2)==nml_pdg))then
          return  ! massless neutrinos -> nothing to do
        elseif ((abs(id1)==ntl_pdg) .and. (abs(id2)==ntl_pdg))then
          return  ! massless neutrinos -> nothing to do
        else
          write(6,*) 'Unexpected Z decay products: id1,id2=',id1,id2
          stop
        endif
        beta=sqrt(max(zip,one-4._dp*m**2/(2._dp*dot(p,j1,j2))))
        p1(:)=(one+beta)/2._dp*p(j1,:)+(one-beta)/2._dp*p(j2,:)
        p2(:)=(one+beta)/2._dp*p(j2,:)+(one-beta)/2._dp*p(j1,:)
        p(j1,:)=p1(:)
        p(j2,:)=p2(:)
      elseif (abs(idparent) == wp_pdg) then
c--- W parent
        if     ((abs(id1)==el_pdg) .and. (abs(id2)==nel_pdg)) then
          m=mel   ! electron
          imass=1
        elseif ((abs(id1)==ml_pdg) .and. (abs(id2)==nml_pdg)) then
          m=mmu   ! muon
          imass=1
        elseif ((abs(id1)==tl_pdg) .and. (abs(id2)==ntl_pdg)) then
          m=mtau  ! tau
          imass=1
        elseif ((abs(id1)==nel_pdg) .and. (abs(id2)==el_pdg)) then
          m=mel   ! electron
          imass=2
        elseif ((abs(id1)==nml_pdg) .and. (abs(id2)==ml_pdg)) then
          m=mmu   ! muon
          imass=2
        elseif ((abs(id1)==ntl_pdg) .and. (abs(id2)==tl_pdg)) then
          m=mtau  ! tau
          imass=2
        else
          write(6,*) 'Unexpected W decay products: id1,id2=',id1,id2
          stop
        endif
        beta=m**2/(2._dp*dot(p,j1,j2))
        if (imass == 1) then ! massive particle is j1
          p2(:)=(one-beta)*p(j2,:)
          p1(:)=p(j1,:)+beta*p(j2,:)
        endif
        if (imass == 2) then ! massive particle is j2
          p1(:)=(one-beta)*p(j1,:)
          p2(:)=p(j2,:)+beta*p(j1,:)
        endif
          p(j1,:)=p1(:)
          p(j2,:)=p2(:)
      else
        write(6,*) 'Unexpected parent in fixmasses: idparent=',idparent
        stop
      endif

      return
      end


c--- Routine to convert MCFM plabel to PDG number
      function plabeltoPDG(ch)
       implicit none
      include 'types.f'
      integer:: plabeltoPDG

      include 'montecarlorpp.f'
      character*2 ch

      if     (ch == 'el') then
        plabeltoPDG=el_pdg
      elseif (ch == 'ea') then
        plabeltoPDG=ea_pdg
      elseif (ch == 'ml') then
        plabeltoPDG=ml_pdg
      elseif (ch == 'ma') then
        plabeltoPDG=ma_pdg
      elseif (ch == 'tl') then
        plabeltoPDG=tl_pdg
      elseif (ch == 'ta') then
        plabeltoPDG=ta_pdg
      elseif (ch == 'bq') then
        plabeltoPDG=bq_pdg
      elseif (ch == 'ba') then
        plabeltoPDG=bb_pdg
      elseif (ch == 'nl') then
        plabeltoPDG=nel_pdg
      elseif (ch == 'na') then
        plabeltoPDG=nea_pdg
      elseif (ch == 'nm') then
        plabeltoPDG=nml_pdg
      elseif (ch == 'bm') then
        plabeltoPDG=nma_pdg
      elseif (ch == 'nt') then
        plabeltoPDG=ntl_pdg
      elseif (ch == 'bt') then
        plabeltoPDG=nta_pdg
      elseif (ch == 'pp') then
        plabeltoPDG=g_pdg
      elseif (ch == 'ga') then
        plabeltoPDG=g0_pdg
      else
        write(6,*) 'Unknown plabel in plabeltoPDG: ',ch
        stop
      endif

      return
      end


      subroutine mcfm_getflavour(msq,j,k)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),msqsum,ran2,ptr

c--- total of absolute weights
      msqsum=zip
      do j=-nf,nf
      do k=-nf,nf
        msqsum=msqsum+abs(msq(j,k))
      enddo
      enddo

c--- random weight between 0 and this total
      ptr=msqsum*ran2()

c--- recover corresponding j,k
      msqsum=zip
      do j=-nf,nf
      do k=-nf,nf
        msqsum=msqsum+abs(msq(j,k))
        if (msqsum >= ptr) return
      enddo
      enddo

      return
      end


