      subroutine determinefilenames(inputfile,workdir)
      implicit none
      include 'types.f'
      
      character*72 workdir
      character*72 inputfile,getinput
      integer:: nargs,lenocc,lenarg
c--- work out the name of the input file and return it


      nargs=iargc()
      if (nargs >= 1) then
        call getarg(1,inputfile)
      else
        inputfile='input.DAT'
      endif
      
      lenarg=lenocc(inputfile)

      if ((lenarg<4).or.(inputfile(lenarg-3:lenarg).ne.'.DAT')) then
        workdir=inputfile
c--- truncate if the directory / is included
        if (workdir(lenarg:lenarg) == '/') then
           workdir(lenarg:lenarg)=' '
           lenarg=lenarg-1
        endif
        if (nargs >= 2) then
          call getarg(2,getinput)
          inputfile=workdir(1:lenarg)//'/'//getinput
        else  
          inputfile=workdir(1:lenarg)//'/input.DAT'
        endif
      else
        workdir=''
      endif
      return
      end      


c      write(6,*) '* Using input file named ',inputfile
