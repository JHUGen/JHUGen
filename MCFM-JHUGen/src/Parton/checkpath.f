      character*72 function checkpath(inpstr)
c--- take an input string (path to PDF file) and strip the leading
c--- directory if the parameter "vanillafiles" is true
      implicit none
      include 'vanillafiles.f'
      character inpstr*(*)
      
      if (vanillafiles) then
        checkpath=inpstr(9:)
      else
        checkpath=inpstr
      endif
      
      return
      end
      
