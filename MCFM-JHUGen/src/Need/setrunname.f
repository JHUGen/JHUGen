      subroutine setrunname(scalestart,fscalestart) 
      implicit none
      include 'flags.f'
      include 'masses.f'
      include 'process.f'
      include 'jetcuts.f'
      include 'werkdir.f'
      include 'pdlabel.f'
      include 'part.f'
      include 'vdecayid.f'
      include 'runstring.f'
      double precision scalestart,fscalestart
      integer nlength,lenocc
      character*255 outlabel1,runname,outlabeltmp
      character*3 strmh,getstr,strpt
      character*7 strscale
      common/runname/runname
      common/nlength/nlength

c      if (abs(scalestart-fscalestart) .lt. 1d0) then
c--- if the scales are the same, use this scale as the label
c        strscale=getstr(int(scalestart))
c       else
c--- .... otherwise, use the percentage of (muR/muF)
c        strscale=getstr(int(scalestart/fscalestart*100d0))
c      endif
      strscale=getstr(int(scalestart))//'_'//getstr(int(fscalestart))

      if     ( (case .eq. 'WHbbar')
     .    .or. (case .eq. 'ZHbbar')
     .    .or. (case .eq. 'qq_tth')
     .    .or. (case .eq. 'tottth')
     .    .or. (case .eq. 'HWW_4l')
     .    .or. (case .eq. 'HWW_tb')
     .    .or. (case .eq. 'HWWint')
     .    .or. (case .eq. 'HWWH+i')
     .    .or. (case .eq. 'ggWW4l')
     .    .or. (case .eq. 'ggVV4l')
     .    .or. (case .eq. 'HZZ_4l')
     .    .or. (case .eq. 'HZZ_tb')
     .    .or. (case .eq. 'HZZint')
     .    .or. (case .eq. 'HZZH+i')
     .    .or. (case .eq. 'ggZZ4l')
     .    .or. (case .eq. 'ggfus0')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'ggfus2')
     .    .or. (case .eq. 'ggfus3') ) then
        strmh=getstr(int(hmass))
        outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//
     .   '_'//strmh
      elseif (  (case .eq. 'H_1jet') ) then
        strmh=getstr(int(hmass))
        strpt=getstr(int(ptjetmin))
        outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//
     .   '_'//strmh//'_pt'//strpt(1:2)      
      elseif ( (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet') ) then
        if     (Gflag .eqv. .false.) then
          outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//'_qrk'
        elseif (Qflag .eqv. .false.) then
          outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//'_glu'
        else
          outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale
        endif
      else
        outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale
      endif
      
      nlength=lenocc(outlabel1)
      if (vdecayid) then
        runname=outlabel1(1:nlength)//'_'//v34id//v56id//'_'//runstring
      else
        runname=outlabel1(1:nlength)//'_'//runstring
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
c returns a string of length 3 from an integer
      integer no,i1,i2,i3,zero
      
      zero=ichar('0')

      i1=abs(no)/100
      i2=(abs(no)-i1*100)/10
      i3=abs(no)-i1*100-i2*10

      if    (i1.eq.0.and.i2.eq.0) then
        if (no .lt. 0) then
        getstr='-'//char(i3+zero)//'_'
        else
        getstr=char(i3+zero)//'__'
        endif
      elseif(i1.eq.0) then
        getstr=char(i2+zero)//char(i3+zero)//'_'
      else
        getstr=char(i1+zero)//char(i2+zero)//char(i3+zero)
      endif
      
      return
      end

