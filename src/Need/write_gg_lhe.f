!=== C.Williams
!=== this routine writes out events in LHE format for
!=== gg initiated processes
      subroutine write_gg_lhe(p,xfac)
      implicit none
      include 'constants.f'
      include 'heprup.f'
      include 'hepeup.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'process.f'
      include 'maxwt.f'
      include 'eventbuffer.f'
      include 'runstring.f'
      include 'nproc.f'
      double precision p(mxpart,4),xfac
      integer j,k,i,nu,it
      double precision dot
      logical first,store_wt
      common/store_wt/store_wt
      data first /.true./
      save first
      character*200 outputstring
      integer iruns

c      store_wt=.true.
c      if(store_wt.eqv..false.) then
c         xwgtup=+1
c      else
c      endif

c--- increment event counter
       numstored=numstored+1

c--- Event weight must be multiplied by xfac to account for
c--- events with negative weight or events with weights that exceed wtmax
       xwgtup=wtmax*xfac
       xmaxup(1)=wtmax


      if(first) then
         first=.false.
         iruns=LEN_TRIM(runstring)
         outputstring=runstring(1:iruns)//'.lhe'
         open(unit=83,file=outputstring,status='unknown')
         call init_lhe_events(83)
      endif

!===== Check we are doing an allowed case
      if ((case.eq.'HZZ_4l') .or. (case.eq.'HZZ_tb')
     &.or.(case.eq.'HZZint') .or. (case.eq.'ggZZ4l')
     &.or.(case.eq.'HZZH+i')) then
         nup=npart+4
         idprup=10000+nproc
         scalup=facscale
         aqedup=-1d0
         aqcdup=as
      else
         write(6,*) 'LHE output not supported for case=',case
         write(6,*) 'Check write_gg_lhe.f for details'
         stop
      endif


      idup(1)=21
      idup(2)=21
!======Z^0
         idup(3)=23
!=======mu-
         idup(4)=11
!========mu+
         idup(5)=-11
!======Z^0
         idup(6)=23
!=======mu-
         idup(7)=13
!========mu+
         idup(8)=-13

         istup(1)=-1
         istup(2)=-1

         istup(3)=2
         istup(6)=2

         istup(4)=1
         istup(5)=1
         istup(7)=1
         istup(8)=1

         mothup(1,3)=1
         mothup(2,3)=2
         mothup(1,4)=3
         mothup(2,4)=3
         mothup(1,5)=3
         mothup(2,5)=3

         mothup(1,6)=1
         mothup(2,6)=2
         mothup(1,7)=3
         mothup(2,7)=3
         mothup(1,8)=3
         mothup(2,8)=3


         icolup(1,1)=501
         icolup(2,1)=511
         icolup(1,2)=511
         icolup(2,2)=501

         it=1
         do i=1,nup
            if((i.ne.3).and.(i.ne.6)) then
               do nu=1,4
                  pup(nu,i)=p(it,nu)
                  if(i.lt.3) pup(nu,i)=-pup(nu,i)
               enddo
               it=it+1
               pup(5,i)=0d0
            elseif(i.eq.3) then
               do nu=1,4
                  pup(nu,i)=p(3,nu)+p(4,nu)
               enddo
               pup(5,i)=dsqrt(dot(p,3,4)*2d0)
            elseif(i.eq.6) then
               do nu=1,4
                  pup(nu,i)=p(5,nu)+p(6,nu)
               enddo
               pup(5,i)=dsqrt(dot(p,5,6)*2d0)
            endif
         enddo

         vtimup(:)=0d0
         spinup(:)=9d0

      call lhefwritev(83)

      return
      end


      subroutine init_lhe_events(iu)
      implicit none
      include 'constants.f'
      include 'heprup.f'
      include 'hepeup.f'
      include 'maxwt.f'
      include 'xs_store_info.f'
      include 'energy.f'
      integer nproc
      common/nproc/nproc
      integer iu,ih1,ih2
      common/density/ih1,ih2

      idbmup(1)=isign(2212,ih1)
      idbmup(2)=isign(2212,ih2)
      ebmup(1)=sqrts/2d0
      ebmup(2)=sqrts/2d0

      pdfgup(1)=-1
      pdfgup(2)=-1

      pdfsup(1)=-1
      pdfsup(2)=-1

      idwtup=+3
      nprup=1

      if(xs_store_val.gt.0d0) then
         xsecup(:)=xs_store_val
         xerrup(:)=xs_err_store_val
         xmaxup(:)=wtmax
      else
         xsecup(:)=1d0
         xerrup(:)=1d0
         xmaxup(:)=1d0
      endif

      lprup=nproc+10000

      call lhefwritehdr(iu)
      return
      end

c--- Adapted from a routine of the same name in POWHEG-BOX, P. Nason et al
c--- Original routine available at http://powhegbox.mib.infn.it/
      subroutine lhefwritehdr(nlf)
      implicit none
      include 'heprup.f'
      include 'hepeup.f'
      include 'codeversion.f'
      include 'xs_store_info.f'
      integer nlf
      integer ipr,iran,n1ran,n2ran
      integer idum
      common/ranno/idum
      save/ranno/

      write(nlf,'(a)') '<LesHouchesEvents version="1.0">'
      write(nlf,'(a)') '<!--'
      write(nlf,'(a)') 'file generated with MCFM version '//
     &   codeversion
      write(nlf,'(a)') 'Input file input.DAT contained:'
      call writeinfo(nlf,'# ',xs_store_val,xs_err_store_val,0)
      write(nlf,'(a)') 'End of input.DAT content'
      write(nlf,*) 'Random number generator initialized with: ',
     &  idum
      write(nlf,'(a)') '-->'
      write(nlf,'(a)') '<init>'
      write(nlf,110) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
      do ipr=1,nprup
         write(nlf,120) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &        lprup(ipr)
      enddo

      write(nlf,'(a)') '</init>'

      return

 110  format(1p,2(1x,i8),2(1x,e12.5),6(1x,i6))
 120  format(1p,3(1x,e12.5),1x,i6)

      end


c--- Adapted from a routine of the same name in POWHEG-BOX, P. Nason et al
c--- Original routine available at http://powhegbox.mib.infn.it/
      subroutine lhefwritev(nlf)
      implicit none
      include 'heprup.f'
      include 'hepeup.f'
      integer nlf
      integer i,j

      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
      enddo

      write(nlf,'(a)')'</event>'

      return

 210  format(1p,2(1x,i6),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)

      end


      subroutine lhefwritefooter(nlf)
      implicit none
      integer nlf

      write(nlf,'(a)') '</LesHouchesEvents>'

      return
      end

