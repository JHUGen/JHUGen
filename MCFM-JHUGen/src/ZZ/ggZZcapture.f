      subroutine ggZZcapture(label,h34,h56,j1,j2,j3,j4,j5,j6,
     & amp0,amp2,amp4)
      implicit none
      include 'types.f'
c--- This is called from within the coefficient-computing routines;
c--- captures m^0, m^2 and m^4 parts of the coefficient and stores
c--- in a common block
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ZZclabels.f'
      include 'ZZdlabels.f'
      include 'first.f'
      integer:: imt0,imt2,imt4,imp,ipp,h34,h56,j1,j2,j3,j4,j5,j6,
     & htag,ltag,itag
      complex(dp):: res(2,4,10,3),amp0,amp2,amp4
      character*(*) label
      parameter(imt0=1,imt2=2,imt4=3)
      parameter(ipp=1,imp=2)
      common/ggZZcaptureres/res
!$omp threadprivate(/ggZZcaptureres/)
      
      if (first) then
        res(:,:,:,:)=czip
        first=.false.
      endif
      
      htag=-1
     
c      write(6,*) 'label=',label
 
c--- Canonical permutation:  1 2 3 4 5 6      
      if (   (j1 == 1) .and. (j2 == 2) .and. (j3 == 3) 
     & .and. (j4 == 4) .and. (j5 == 5) .and. (j6 == 6) ) then
        if ((trim(label)=='1x2').and.(h34==1).and.(h56==1)) then
          htag=imp
          ltag=3
          itag=c1_2
        endif
        if ((trim(label)=='1x2pp').and.(h34==1).and.(h56==1)) then
          htag=ipp
          ltag=3
          itag=c1_2
        endif
        if ((trim(label)=='1x34').and.(h34==1).and.(h56==1)) then
          htag=imp
          ltag=3
          itag=c1_34
        endif
        if ((trim(label)=='1x34pp').and.(h34==1).and.(h56==1))then
          htag=ipp
          ltag=3
          itag=c1_34
        endif
        if ((trim(label)=='12x34').and.(h34==1).and.(h56==1)) then
          htag=imp
          ltag=3
          itag=c12_34
        endif
        if ((trim(label)=='12x34pp').and.(h34==1).and.(h56==1))
     &    then
          htag=ipp
          ltag=3
          itag=c12_34
        endif
        if ((trim(label)=='1x34x2mp').and.(h34==1).and.(h56==1))
     &    then
          htag=imp
          ltag=4
          itag=d1_34_2
        endif
        if ((trim(label)=='1x34x2pp').and.(h34==1).and.(h56==1))
     &    then
          htag=ipp
          ltag=4
          itag=d1_34_2
        endif
        if ((trim(label)=='2x1x34pp').and.(h34==1).and.(h56==1))
     &    then
          htag=ipp
          ltag=4
          itag=d2_1_34
        endif
        if ((trim(label)=='2x1x34').and.(h34==1).and.(h56==1))then
          htag=imp
          ltag=4
          itag=d2_1_34
        endif
        if ((trim(label)=='d62x1x34').and.(h34==1).and.(h56==1))
     &    then
          htag=imp
          ltag=4
          itag=d6_2_1_34
        endif
        if ((label(1:5)=='bubmp').and.(h34==1).and.(h56==1))
     &    then
          htag=imp
          ltag=2
          itag=ichar(label(6:6))-ichar('0')
        endif
        if ((label(1:5)=='bubpp').and.(h34==1).and.(h56==1))
     &    then
          htag=ipp
          ltag=2
          itag=ichar(label(6:6))-ichar('0')
        endif
      endif

c--- Permute 1 and 2:  2 1 3 4 5 6       
      if (   (j1 == 2) .and. (j2 == 1) .and. (j3 == 3) 
     & .and. (j4 == 4) .and. (j5 == 5) .and. (j6 == 6) ) then
        if ((trim(label)=='1x34').and.(h34==1).and.(h56==1)) then
          htag=imp
          ltag=3
          itag=c2_34
          if (amp0 .ne. czip) amp0=conjg(amp0)
          if (amp2 .ne. czip) amp2=conjg(amp2)
          if (amp4 .ne. czip) amp4=conjg(amp4)
        endif
        if ((trim(label)=='1x34pp').and.(h34==1).and.(h56==1))then
          htag=ipp
          ltag=3
          itag=c2_34
          if (amp0 .ne. czip) amp0=conjg(amp0)
          if (amp2 .ne. czip) amp2=conjg(amp2)
          if (amp4 .ne. czip) amp4=conjg(amp4)
        endif
        if ((trim(label)=='2x1x34').and.(h34==1).and.(h56==1))then
          htag=imp
          ltag=4
          itag=d1_2_34
          if (amp0 .ne. czip) amp0=conjg(amp0)
          if (amp2 .ne. czip) amp2=conjg(amp2)
          if (amp4 .ne. czip) amp4=conjg(amp4)
        endif
        if ((trim(label)=='d62x1x34').and.(h34==1).and.(h56==1))
     &    then
          htag=imp
          ltag=4
          itag=d6_1_2_34
          if (amp0 .ne. czip) amp0=conjg(amp0)
          if (amp2 .ne. czip) amp2=conjg(amp2)
          if (amp4 .ne. czip) amp4=conjg(amp4)
        endif
      endif
       
c--- Permute 3,4 and 5,6:  1 2 5 6 3 4
      if (   (j1 == 1) .and. (j2 == 2) .and. (j3 == 5) 
     & .and. (j4 == 6) .and. (j5 == 3) .and. (j6 == 4) ) then
        if ((trim(label)=='1x34').and.(h34==1).and.(h56==1)) then
          htag=imp
          ltag=3
          itag=c1_56
        endif
        if ((trim(label)=='1x34pp').and.(h34==1).and.(h56==1))then
          htag=ipp
          ltag=3
          itag=c1_56
        endif
        if ((trim(label)=='2x1x34pp').and.(h34==1).and.(h56==1))
     &    then
          htag=ipp
          ltag=4
          itag=d1_2_34
        endif
      endif
       
c--- Permute 1,3,4 and 2,5,6:  2 1 3 4 5 6       
      if (   (j1 == 2) .and. (j2 == 1) .and. (j3 == 5) 
     & .and. (j4 == 6) .and. (j5 == 3) .and. (j6 == 4) ) then
        if ((trim(label)=='1x34').and.(h34==1).and.(h56==1)) then
          htag=imp
          ltag=3
          itag=c2_56
          if (amp0 .ne. czip) amp0=conjg(amp0)
          if (amp2 .ne. czip) amp2=conjg(amp2)
          if (amp4 .ne. czip) amp4=conjg(amp4)
        endif
        if ((trim(label)=='1x34pp').and.(h34==1).and.(h56==1))then
          htag=ipp
          ltag=3
          itag=c2_56
          if (amp0 .ne. czip) amp0=conjg(amp0)
          if (amp2 .ne. czip) amp2=conjg(amp2)
          if (amp4 .ne. czip) amp4=conjg(amp4)
        endif
      endif
       
      if (htag < 0) then
c        write(6,*) 'ggZZcapture: no match found'
c        stop
        return
      endif

c      write(6,*) 'label',label
c      write(6,*) 'htag=',htag
c      write(6,*) 'ltag=',ltag
c      write(6,*) 'itag=',itag
c      pause 'got one'

c--- assign results      
      res(htag,ltag,itag,imt0)=amp0
      res(htag,ltag,itag,imt2)=amp2
      res(htag,ltag,itag,imt4)=amp4
      
      if (abs(res(htag,ltag,itag,imt0)) < 1d-8)
     &    res(htag,ltag,itag,imt0)=czip
      if (abs(res(htag,ltag,itag,imt2)) < 1d-8)
     &    res(htag,ltag,itag,imt2)=czip
      if (abs(res(htag,ltag,itag,imt4)) < 1d-8)
     &    res(htag,ltag,itag,imt4)=czip

      return
      end
      
