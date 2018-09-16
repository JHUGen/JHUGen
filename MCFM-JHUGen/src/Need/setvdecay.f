      subroutine setvdecay(idv,vcharge)
      implicit none
      include 'types.f'
c--- Set up particle labels and couplings for vector boson
c--- decays according to contents of vdecayid common block
c---
c--- Vector boson is p(idv) with charge vcharge
c---
c--- If vdecayid = FALSE, returns without doing anything
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'plabel.f'
      include 'vdecayid.f'
      include 'zcouple.f'
      include 'nqcdjets.f'
      integer:: idv,vcharge
      character*2 vdecay,decay1,decay2
      character*15 decaystring
      real(dp):: decayq,decayl,decayr

c--- return if decay not specified in input file
      if (vdecayid .eqv. .false.) return
      
      if     (idv == 34) then
        vdecay=v34id
      elseif (idv == 56) then
        vdecay=v56id
      else
        write(6,*) 'Subroutine setvdecay called improperly;'
        write(6,*) 'only idv=34, 56 allowed but idv=',idv
        stop
      endif

c--- allowed string for specifying vector boson decays,
c---  (particle, anti-particle)

c--- Z decays
      if (vcharge == 0) then
      if     ((vdecay == 'el') .or. (vdecay == 'EL')) then
        decaystring='(e-, e+)       '
        decay1='el'
        decay2='ea'
        decayq=-1._dp
        decayl=le
        decayr=re
      elseif ((vdecay == 'mu') .or. (vdecay == 'MU')
     &   .or. (vdecay == 'ml') .or. (vdecay == 'ML')) then
        decaystring='(mu-, mu+)     '
        decay1='ml'
        decay2='ma'
        decayq=-1._dp
        decayl=le
        decayr=re
      elseif ((vdecay == 'tl') .or. (vdecay == 'TL')) then
        decaystring='(tau-, tau+)   '
        decay1='tl'
        decay2='ta'
        decayq=-1._dp
        decayl=le
        decayr=re
      elseif ((vdecay == 'nu') .or. (vdecay == 'NU')
     &   .or. (vdecay == 'nl') .or. (vdecay == 'NL')) then
        decaystring='(nu, nubar) x 3'
        decay1='nl'
        decay2='na'
        decayq=0._dp
        decayl=ln*sqrt(3._dp)
        decayr=rn*sqrt(3._dp)
      elseif ((vdecay == 'bq') .or. (vdecay == 'BQ')) then
        decaystring='(b, b-bar)     '
        decay1='bq'
        decay2='ba'
        nqcdjets=nqcdjets+2
        decayq=Q(5)*sqrt(xn)
        decayl=l(5)*sqrt(xn)
        decayr=r(5)*sqrt(xn)
      else
        write(6,*) 'Decay string not recognized, vdecay=',vdecay
        stop
      endif
      endif

c--- W+ decays      
      if(vcharge == +1) then
      if     ((vdecay == 'en') .or. (vdecay == 'EN')) then
        decaystring='(ve, e+)       '
        decay1='nl'
        decay2='ea'
      elseif ((vdecay == 'mn') .or. (vdecay == 'MN')) then
        decaystring='(vmu, mu+)     '
        decay1='nm'
        decay2='ma'
      elseif ((vdecay == 'tn') .or. (vdecay == 'TN')) then
        decaystring='(vtau, tau+)   '
        decay1='nt'
        decay2='ta'
      else
        write(6,*) 'Decay string not recognized, vdecay=',vdecay
        stop
      endif
      endif
      
c--- W- decays      
      if(vcharge == -1) then
      if     ((vdecay == 'en') .or. (vdecay == 'EN')) then
        decaystring='(e-, ve~)      '
        decay1='el'
        decay2='na'
      elseif ((vdecay == 'mn') .or. (vdecay == 'MN')) then
        decaystring='(mu-, vmu~)    '
        decay1='ml'
        decay2='bm'
      elseif ((vdecay == 'tn') .or. (vdecay == 'TN')) then
        decaystring='(tau-, vtau~)  '
        decay1='tl'
        decay2='bt'
      else
        write(6,*) 'Decay string not recognized, vdecay=',vdecay
        stop
      endif
      endif
      
      if (idv == 34) then
        plabel(3)=decay1
        plabel(4)=decay2
        if (vcharge == 0) then
          q1=decayq
          l1=decayl
          r1=decayr
        endif
      endif

      if (idv == 56) then
        plabel(5)=decay1
        plabel(6)=decay2
        if (vcharge == 0) then
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
      
