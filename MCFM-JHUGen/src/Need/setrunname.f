      subroutine setrunname(scalestart,fscalestart) 
      implicit none
      include 'types.f'
      include 'flags.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'jetcuts.f'
      include 'werkdir.f'
      include 'pdlabel.f'
      include 'kpart.f'
      include 'vdecayid.f'
      include 'runstring.f'
      include 'dynamicscale.f'
      include 'taucut.f'
      real(dp):: scalestart,fscalestart
      integer:: nlength,lenocc
      character*255 outlabel1,runname,outlabeltmp
      character*3 strmh,getstr,strpt
      character*5 strtaucut
      character*9 strscale
      character*15 part,kpartstring
      character*6 case,kcasestring
      common/runname/runname
      common/nlength/nlength

c      if (abs(scalestart-fscalestart) < 1._dp) then
c--- if the scales are the same, use this scale as the label
c        strscale=getstr(int(scalestart))
c       else
c--- .... otherwise, use the percentage of (muR/muF)
c        strscale=getstr(int(scalestart/fscalestart*100._dp))
c      endif
      strscale=getstr(int(scalestart))//'__'
     &       //getstr(int(fscalestart))//'_'

      if (dynamicscale) then
        write(strscale,'(F4.2,"_",F4.2)') scalestart,fscalestart
      endif

c--- convert kpart and kcase to strings
      part=kpartstring(kpart)
      case=kcasestring(kcase)

      if     ( (kcase==kWHbbar)
     &    .or. (kcase==kZHbbar)
     &    .or. (kcase==kqq_tth)
     &    .or. (kcase==ktottth)
     &    .or. (kcase==kHWW_4l)
     &    .or. (kcase==kHWW_tb)
     &    .or. (kcase==kHWWint)
     &    .or. (kcase==kHWWHpi)
     &    .or. (kcase==kggWW4l)
     &    .or. (kcase==kggVV4l)
     &    .or. (kcase==kHZZ_4l)
     &    .or. (kcase==kHZZ_tb)
     &    .or. (kcase==kHZZint)
     &    .or. (kcase==kHZZHpi)
     &    .or. (kcase==kggZZ4l)
     &    .or. (kcase==kggfus0)
     &    .or. (kcase==kggfus1)
     &    .or. (kcase==kggfus2)
     &    .or. (kcase==kggfus3) ) then
        strmh=getstr(int(hmass))
        outlabel1=case//'_'//trim(part)//'_'//pdlabel//'_'//strscale//
     &   '_'//strmh
      elseif (  (kcase==kH_1jet) ) then
        strmh=getstr(int(hmass))
        strpt=getstr(int(ptjetmin))
        outlabel1=case//'_'//trim(part)//'_'//pdlabel//'_'//strscale//
     &   '_'//strmh//'_pt'//strpt(1:2)      
      elseif ( (kcase==kW_2jet)
     &    .or. (kcase==kZ_2jet) ) then
        if     (Gflag .eqv. .false.) then
          outlabel1=case//'_'//trim(part)//'_'//pdlabel//'_'//strscale//'_qrk'
        elseif (Qflag .eqv. .false.) then
          outlabel1=case//'_'//trim(part)//'_'//pdlabel//'_'//strscale//'_glu'
        else
          outlabel1=case//'_'//trim(part)//'_'//pdlabel//'_'//strscale
        endif
      else
        outlabel1=case//'_'//trim(part)//'_'//pdlabel//'_'//strscale
      endif
      
      nlength=lenocc(outlabel1)
      if (usescet) then
        write(strtaucut,'(ES5.0E1)') taucut
        runname=outlabel1(1:nlength)//'_'//strtaucut//'_'//runstring
      else
        if (vdecayid) then
         runname=outlabel1(1:nlength)//'_'//v34id//v56id//'_'//runstring
        else
         runname=outlabel1(1:nlength)//'_'//runstring
        endif
      endif
      nlength=lenocc(runname)

c--- add working directory, if necessary 
      if (werkdir .ne. '') then
        outlabeltmp=runname
        runname=werkdir(1:lenocc(werkdir))//'/'//outlabeltmp
        nlength=nlength+1+lenocc(werkdir)
      endif
      
      return
      end


      character*3 function getstr(no)
c returns a string of length 3 from an integer::
      integer:: no,i1,i2,i3,izero
      
      izero=ichar('0')

      i1=abs(no)/100
      i2=(abs(no)-i1*100)/10
      i3=abs(no)-i1*100-i2*10

      if    (i1==0.and.i2==0) then
        if (no < 0) then
        getstr='-'//char(i3+izero)//'_'
        else
        getstr=char(i3+izero)//'__'
        endif
      elseif(i1==0) then
        getstr=char(i2+izero)//char(i3+izero)//'_'
      else
        getstr=char(i1+izero)//char(i2+izero)//char(i3+izero)
      endif
      
      return
      end

