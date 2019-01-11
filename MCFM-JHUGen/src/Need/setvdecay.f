      subroutine setvdecay(idv,vcharge)
c--- Set up particle labels and couplings for vector boson
c--- decays according to contents of vdecayid common block
c---
c--- Vector boson is p(idv) with charge vcharge
c---
c--- If vdecayid = FALSE, returns without doing anything
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      include 'plabel.f'
      include 'vdecayid.f'
      include 'zcouple.f'
      include 'nqcdjets.f'
      integer idv,vcharge
      character*2 vdecay,decay1,decay2
      character*15 decaystring
      double precision decayq,decayl,decayr

c--- return if decay not specified in input file
      if (vdecayid .eqv. .false.) return
      
      if     (idv .eq. 34) then
        vdecay=v34id
      elseif (idv .eq. 56) then
        vdecay=v56id
      else
        write(6,*) 'Subroutine setvdecay called improperly;'
        write(6,*) 'only idv=34, 56 allowed but idv=',idv
        stop
      endif

c--- allowed string for specifying vector boson decays,
c---  (particle, anti-particle)

c--- Z decays
      if (vcharge .eq. 0) then
      if     ((vdecay .eq. 'el') .or. (vdecay .eq. 'EL')) then
        decaystring='(e-, e+)       '
        decay1='el'
        decay2='ea'
        decayq=-1d0
        decayl=le
        decayr=re
      elseif ((vdecay .eq. 'mu') .or. (vdecay .eq. 'MU')
     &   .or. (vdecay .eq. 'ml') .or. (vdecay .eq. 'ML')) then
        decaystring='(mu-, mu+)     '
        decay1='ml'
        decay2='ma'
        decayq=-1d0
        decayl=le
        decayr=re
      elseif ((vdecay .eq. 'tl') .or. (vdecay .eq. 'TL')) then
        decaystring='(tau-, tau+)   '
        decay1='tl'
        decay2='ta'
        decayq=-1d0
        decayl=le
        decayr=re
      elseif ((vdecay .eq. 'nu') .or. (vdecay .eq. 'NU')
     &   .or. (vdecay .eq. 'nl') .or. (vdecay .eq. 'NL')) then
        decaystring='(nu, nubar) x 3'
        decay1='nl'
        decay2='na'
        decayq=0d0
        decayl=ln*dsqrt(3d0)
        decayr=rn*dsqrt(3d0)
      elseif ((vdecay .eq. 'bq') .or. (vdecay .eq. 'BQ')) then
        decaystring='(b, b-bar)     '
        decay1='bq'
        decay2='ba'
        nqcdjets=nqcdjets+2
        decayq=Q(5)*dsqrt(xn)
        decayl=l(5)*dsqrt(xn)
        decayr=r(5)*dsqrt(xn)
      else
        write(6,*) 'Decay string not recognized, vdecay=',vdecay
        stop
      endif
      endif

c--- W+ decays      
      if(vcharge .eq. +1) then
      if     ((vdecay .eq. 'en') .or. (vdecay .eq. 'EN')) then
        decaystring='(ve, e+)       '
        decay1='nl'
        decay2='ea'
      elseif ((vdecay .eq. 'mn') .or. (vdecay .eq. 'MN')) then
        decaystring='(vmu, mu+)     '
        decay1='nm'
        decay2='ma'
      elseif ((vdecay .eq. 'tn') .or. (vdecay .eq. 'TN')) then
        decaystring='(vtau, tau+)   '
        decay1='nt'
        decay2='ta'
      else
        write(6,*) 'Decay string not recognized, vdecay=',vdecay
        stop
      endif
      endif
      
c--- W- decays      
      if(vcharge .eq. -1) then
      if     ((vdecay .eq. 'en') .or. (vdecay .eq. 'EN')) then
        decaystring='(e-, ve~)      '
        decay1='el'
        decay2='na'
      elseif ((vdecay .eq. 'mn') .or. (vdecay .eq. 'MN')) then
        decaystring='(mu-, vmu~)    '
        decay1='ml'
        decay2='bm'
      elseif ((vdecay .eq. 'tn') .or. (vdecay .eq. 'TN')) then
        decaystring='(tau-, vtau~)  '
        decay1='tl'
        decay2='bt'
      else
        write(6,*) 'Decay string not recognized, vdecay=',vdecay
        stop
      endif
      endif
      
      if (idv .eq. 34) then
        plabel(3)=decay1
        plabel(4)=decay2
        if (vcharge .eq. 0) then
          q1=decayq
          l1=decayl
          r1=decayr
        endif
      endif

      if (idv .eq. 56) then
        plabel(5)=decay1
        plabel(6)=decay2
        if (vcharge .eq. 0) then
          q2=decayq
          l2=decayl
          r2=decayr
        endif
      endif
      
      write(6,*)
      write(6,*) '****************************************************'
      write(6,99) idv,decaystring
      write(6,*) '****************************************************'
      
      
      return
      
   99 format('*      Vector boson decay ',i2,' -> ',a15,'     *')
      end
      
