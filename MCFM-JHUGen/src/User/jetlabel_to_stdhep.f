************************************************************************
*  Function that converts a 2-character particle label to an integer::   *
*  corresponding to the PDG/STDHEP representation of that particle     *
************************************************************************
      function jetlabel_to_stdhep(label)
       implicit none
      include 'types.f'
      integer:: jetlabel_to_stdhep
      character*2 label
      
      if     (label == 'nl') then
        jetlabel_to_stdhep=+12
      elseif (label == 'na') then
        jetlabel_to_stdhep=-12
      elseif (label == 'el') then
        jetlabel_to_stdhep=+11
      elseif (label == 'ea') then
        jetlabel_to_stdhep=-11
      elseif (label == 'ml') then
        jetlabel_to_stdhep=+13
      elseif (label == 'ma') then
        jetlabel_to_stdhep=-13
      elseif (label == 'tl') then
        jetlabel_to_stdhep=+15
      elseif (label == 'ta') then
        jetlabel_to_stdhep=-15
      elseif (label == 'bq') then
        jetlabel_to_stdhep=+5
      elseif (label == 'ba') then
        jetlabel_to_stdhep=-5
      elseif (label == 'pp') then
        jetlabel_to_stdhep=21      ! Use gluon label for any jet
      elseif (label == 'ga') then
        jetlabel_to_stdhep=22
      elseif (label == 'ig') then
c--- this is a dummy value: no particle id available      
        jetlabel_to_stdhep=0
      else
        write(6,*) 'Conversion routine jetlabel_to_stdhep.f failed'
        write(6,*) 'I do not know what a [',label,'] label means'
        stop         
      endif
      
      return
      end
