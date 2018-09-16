************************************************************************
*  Function that converts a 2-character particle label to an integer   *
*  corresponding to the PDG/STDHEP representation of that particle     *
************************************************************************
      integer function jetlabel_to_stdhep(label)
      character*2 label
      
      if     (label .eq. 'nl') then
        jetlabel_to_stdhep=+12
      elseif (label .eq. 'na') then
        jetlabel_to_stdhep=-12
      elseif (label .eq. 'el') then
        jetlabel_to_stdhep=+11
      elseif (label .eq. 'ea') then
        jetlabel_to_stdhep=-11
      elseif (label .eq. 'ml') then
        jetlabel_to_stdhep=+13
      elseif (label .eq. 'ma') then
        jetlabel_to_stdhep=-13
      elseif (label .eq. 'tl') then
        jetlabel_to_stdhep=+15
      elseif (label .eq. 'ta') then
        jetlabel_to_stdhep=-15
      elseif (label .eq. 'bq') then
        jetlabel_to_stdhep=+5
      elseif (label .eq. 'ba') then
        jetlabel_to_stdhep=-5
      elseif (label .eq. 'pp') then
        jetlabel_to_stdhep=21      ! Use gluon label for any jet
      elseif (label .eq. 'ga') then
        jetlabel_to_stdhep=22
      elseif (label .eq. 'ig') then
c--- this is a dummy value: no particle id available      
        jetlabel_to_stdhep=0
      else
        write(6,*) 'Conversion routine jetlabel_to_stdhep.f failed'
        write(6,*) 'I do not know what a [',label,'] label means'
        stop         
      endif
      
      return
      end
