      subroutine smartzero
      implicit none
      include 'types.f'
c--- zero out all entries in the temporary histograms used for
c--- binning the weights in the real contribution
      
      integer:: maxd,maxhisto
      parameter(maxd=40,maxhisto=100) ! maximum number of dipoles and histograms
      include 'incsmarthisto.f'
            
      ibin(:,:)=0             ! which bin entry should be filled
      xbinwgt(:,:)=0._dp        ! ... with this weight
      icont(:)=0              ! counter for number of entries in ibin
      
      return
      end
      
      
      subroutine smartbook(N,TIT,DEL,XMIN,XMAX)
      implicit none
      include 'types.f'
      
      integer:: N
      real(dp):: DEL,XMIN,XMAX
      character*(*) TIT
      integer:: maxd,maxhisto
      parameter(maxd=40,maxhisto=100) ! maximum number of dipoles and histograms
      include 'incsmarthisto.f'
      
      binmin(N)=XMIN
      bindel(N)=DEL
      
      return
      end
      
      
      subroutine smartfill(N,X,Y)
      implicit none
      include 'types.f'
      
      integer:: N,I
      real(dp):: X,Y
      integer:: maxd,maxhisto
      parameter(maxd=40,maxhisto=100) ! maximum number of dipoles and histograms
      include 'incsmarthisto.f'
      
      I=INT((X-binmin(N))/bindel(N)+1)
      
      icont(N)=icont(N)+1     ! increase counter by one
      ibin(icont(N),N)=I      ! record bin number for weight
      xbinwgt(icont(N),N)=Y/bindel(N)   ! record actual weight

      return
      end
      
      
      subroutine smartadd(wgt)
      implicit none
      include 'types.f'
c--- add temporary histograms to the cumulative ones
      
      integer:: I,L,M
      include 'histo.f'
      real(dp):: wgt,tmpvar
      integer:: nplotmax
      common/nplotmax/nplotmax
      integer:: maxd
      parameter(maxd=40) ! maximum number of dipoles
      include 'incsmarthisto.f'
ccccc!$omp threadprivate(/nplotmax/)
      do I=1,nplotmax
      
      if (icont(I) == 0) cycle ! skip histogram if no entries
!$omp atomic
        IENT(I)=IENT(I)+1
        do L=1,icont(I)
c--- skip if out of bounds
          if ((ibin(L,I) < 1) .or. (ibin(L,I) > NBIN(I))) cycle
c---    now look for other entries in the same bin before adding
          do M=L+1,icont(I)
          if (ibin(M,I) == ibin(L,I)) then
            ibin(M,I)=0
            xbinwgt(L,I)=xbinwgt(L,I)+xbinwgt(M,I)
          endif
          enddo
          tmpvar=xbinwgt(L,I)**2*HDEL(I)
!$omp atomic
          HIST(I,ibin(L,I))=HIST(I,ibin(L,I))+xbinwgt(L,I)
!$omp atomic
          HIST(maxhisto+I,ibin(L,I))=HIST(maxhisto+I,ibin(L,I))+tmpvar
        enddo
      enddo
      
      return
      end
      
      
