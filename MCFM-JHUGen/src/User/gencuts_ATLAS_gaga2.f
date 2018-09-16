!==== C.W Jan 16, ATLAS cuts for mgg searches, see for instance
!==== ATLAS-CONF-2015-08
      function gencuts_ATLAS_gaga2(p)
      implicit none
      include 'types.f'
      logical gencuts_ATLAS_gaga2
      include 'constants.f'
      include 'mxpart.f'
      include 'first.f'
      include 'mpicommon.f'
      real(dp) :: p(mxpart,4),mgg,ptgaga,pth,pts,ygh,ygs,ygg
      common/gaga_observables/mgg,ptgaga,pth,pts,ygh,ygs,ygg
!$omp threadprivate(/gaga_observables/)
      integer i
      real(dp):: ptfrac1,ptfrac2
      real(dp):: pt,pttwo,yrap,yraptwo
!-----observables for BSM interest
      real(dp):: gaga_highm_ptgaga(0:3),gaga_highm_ptgh(0:3)
      real(dp):: gaga_highm_ptgs(0:3),gaga_highm_ygg(0:3)
      real(dp):: gaga_highm_yh(0:3),gaga_highm_ys(0:3),ayi,pt3,pt4
      integer itag
      common/gaga_highmgg_obs/gaga_highm_ptgaga,gaga_highm_ptgh,
     &     gaga_highm_ptgs,gaga_highm_ygg,gaga_highm_yh,gaga_highm_ys
!$omp threadprivate(/gaga_highmgg_obs/)

      
!=====initalize out of bounds
      gaga_highm_ptgaga(:)=-1._dp
      gaga_highm_ptgh(:)=-1._dp
      gaga_highm_ptgs(:)=-1._dp
      gaga_highm_ygg(:)=-10000._dp
      gaga_highm_yh(:)=-10000._dp
      gaga_highm_ys(:)=-10000._dp
      
      gencuts_ATLAS_gaga2=.false.

      ptfrac1=0.4_dp
      ptfrac2=0.3_dp
!=====basic photon cuts from input file, this routine calculates
!==== mgaga and checks that

!=====pt(ga,hard)/mgaga > ptcut1, and pt(ga,soft)/mgaga > ptcut2

!=====initalize histo observables
      pts=-1._dp
      pth=-1._dp
      mgg=-1._dp
      ptgaga=-1._dp
      ygh=-1000._dp
      ygs=-1000._dp
      ygg=-1000._dp
      
      
      if(first) then
         first=.false.
!$omp master
         if (rank == 0) then
         write(6,*)
         write(6,*) '*********** Additional photon cuts applied *********'
         write(6,*) '*                                                  *'
         write(6,99) '*             pt(hard)/m34 > ',ptfrac1,'              *'
         write(6,99) '*             pt(soft)/m34 > ',ptfrac2,'              *'
         write(6,*) '*                                                  *'
         write(6,*) '******** ATLAS rapidity conditions applied *********'
         write(6,*) '*                                                  *'
         write(6,99) '*            |eta(gamma) | < ',2.37_dp,'              *'
         write(6,*) '*        ATLAS crack (1.37, 1.52) excluded         *'
         write(6,*) '****************************************************'
         endif
!$omp end master
      endif
      
      mgg=(p(3,4)+p(4,4))**2
      do i=1,3
         mgg=mgg-(p(3,i)+p(4,i))**2
      enddo
      mgg=sqrt(mgg)

!-----rapidity and crack measurements
      do i=3,4
         ayi=abs(yrap(i,p))
         if (ayi > 2.37_dp) goto 101
         if ((ayi > 1.37_dp) .and. (ayi < 1.52_dp)) goto 101
      enddo
      
!===== routine also calculates quantities for histogram 
      pt3=pt(3,p)
      pt4=pt(4,p)
      if (pt3 > pt4) then
         pth=pt3
         pts=pt4
         ygh=yrap(3,p)
         ygs=yrap(4,p)
         
      else
         pts=pt3
         pth=pt4
         ygh=yrap(4,p)
         ygs=yrap(3,p)
      endif

      if ((pth/mgg < ptfrac1) .or. (pts/mgg < ptfrac2)) then
         goto 101
      endif

!======calculate stuff for histos
      ptgaga=pttwo(3,4,p)
      ygg=yraptwo(3,4,p) 

      itag=0
!========now organize for mgg histos
      if((mgg >= 200._dp).and.(mgg  <= 700._dp)) then
         itag=1
      elseif((mgg >  700._dp).and.(mgg <=800._dp)) then
         itag=2
      elseif((mgg > 800._dp)) then
         itag=3
      endif

      gaga_highm_ptgaga(0)=ptgaga
      gaga_highm_ptgaga(itag)=ptgaga
      
      gaga_highm_ptgh(0)=pth
      gaga_highm_ptgh(itag)=pth
      gaga_highm_ptgs(0)=pth
      gaga_highm_ptgs(itag)=pts

      gaga_highm_ygg(0)=ygg
      gaga_highm_ygg(itag)=ygg
      
      gaga_highm_yh(0)=ygh
      gaga_highm_yh(itag)=ygh

      gaga_highm_ys(0)=ygs
      gaga_highm_ys(itag)=ygs

      
      return
!=====zero histos if fail
 101  continue
      gencuts_ATLAS_gaga2=.true.
      pts=-1._dp
      pth=-1._dp
      mgg=-1._dp
      ptgaga=-1._dp
      ygh=-1000._dp
      ygs=-1000._dp
      ygg=-1000._dp

      return 

   99 format(1x,a29,f6.2,a17)
      
      end
      
