      subroutine checkversion(inpunit,filename)
      implicit none
      include 'types.f'
************************************************************************
*  Checks that the version of MCFM specified in the next line of unit  *
*  "inpunit" agrees with the version number of the code                *
************************************************************************
      
      include 'codeversion.f'
      integer:: inpunit,j,dat
      character*6 fileversion
      character*(*) filename

      read(inpunit,*) fileversion

      if (fileversion .ne. codeversion) then
        dat=18
        do j=1,20
          if (filename(j:j) == 'D') dat=j
        enddo        
        write(6,*)
        write(6,*) 'Sorry, the version of this input file does not'
        write(6,*) 'match with the code version number. Please refer'
        write(6,*) 'to the documentation and adjust accordingly.'
        write(6,*)
        write(6,*) '     Filename: ',filename(1:dat+2)
        write(6,*) ' File version: ',fileversion
        write(6,*) ' Code version: ',codeversion
        write(6,*)
        stop
      endif
      
      return
      end
      
