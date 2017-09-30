      subroutine bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot) 
      implicit none
      include 'nplot.f'
      include 'part.f'
      include 'outputflags.f'
      include 'vegas_common.f'
      integer n
      character*(*) titlex
      character*3 llplot
      character*4 tag
      double precision var,wt,wt2,xmin,xmax,dx
      logical, save :: threadfirst(500)=.true.
!$omp threadprivate(threadfirst)

      if     (tag .eq. 'book') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call mbook(n,titlex,dx,xmin,xmax)
c--- also book the errors now (in maxhisto+n,2*maxhisto+n)
          call mbook(maxhisto+n,titlex,dx,xmin,xmax)
          call mbook(2*maxhisto+n,titlex,dx,xmin,xmax)
          if ( (part .eq. 'real') .or. (part .eq. 'tota')
     &     .or.(part .eq. 'todk')) then
c            call tmpmbook(n,titlex,dx,xmin,xmax)
            call smartbook(n,titlex,dx,xmin,xmax)
        endif
        else
c--- DSW histograms - call hbook booking routine
          call dswhbook(n,titlex,dx,xmin,xmax)
        endif
      elseif (tag .eq. 'plot') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
c--- also book the errors now (in maxhisto+n); fill temp histos for real
          if ((part .eq. 'lord') .or. (part .eq. 'virt')
     &   .or. (part .eq. 'frag')) then
            call mfill(n,var,wt)
            call mfill(maxhisto+n,var,wt2)
          else
            if (threadfirst(n)) then
c              call tmpmbook(n,titlex,dx,xmin,xmax)
              call smartbook(n,titlex,dx,xmin,xmax)
              threadfirst(n)=.false.
            endif
c            call tmpmfill(n,var,wt)
            call smartfill(n,var,wt)
        endif
        else
c--- DSW histograms - call hbook filling routine
c--- note that we divide by # of iterations here since it is only
c--- handled at the end in the default MCFM histograms
          call dswhfill(n,var,wt/dfloat(itmx))
        endif
        linlog(n)=llplot
        titlearray(n)=titlex
      endif

      return
      end

      subroutine ebookplot(n,tag,var,wt) 
      implicit none
      include 'PDFerrors.f'
      include 'outputflags.f'
      integer n
      double precision var,wt
      character tag*4

      if (PDFerrors .eqv. .false.) return

      if (tag.eq.'book') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call ebook(n)
        endif
      elseif (tag .eq. 'plot') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call efill(n,var,wt)
        endif
      endif

      return
      end

