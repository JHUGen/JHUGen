MODULE ModCommandLine
implicit none

integer, parameter :: maxlogicals=20, maxintegers=20, maxreal8s=20, maxcomplex8s=150, maxstrings=10

interface ReadCommandLineArgument
    module procedure ReadCommandLineArgument_logical, ReadCommandLineArgument_integer, ReadCommandLineArgument_real8,&
                     ReadCommandLineArgument_complex8, ReadCommandLineArgument_string
end interface

interface savevalue
  module procedure :: savevalue_logical, savevalue_integer, savevalue_real8, savevalue_complex8, savevalue_string
end interface savevalue

type SaveValues
  character(len=100) :: logicalnames(1:maxlogicals), integernames(1:maxintegers), real8names(1:maxreal8s), complex8names(1:maxcomplex8s), stringnames(1:maxstrings)
  integer :: nlogicals, nintegers, nreal8s, ncomplex8s, nstrings
  logical :: logicals(1:maxlogicals)
  integer :: integers(1:maxintegers)
  real(8) :: real8s(1:maxreal8s)
  complex(8) :: complex8s(1:maxcomplex8s)
  character(len=100) :: strings(1:maxstrings)
contains
  procedure :: savevalue_logical
  procedure :: savevalue_integer
  procedure :: savevalue_real8
  procedure :: savevalue_complex8
  procedure :: savevalue_string
  procedure :: sort => SortSaveValues
  procedure :: WriteToFile => WriteSaveValuesToFile
  procedure :: ReadFromFile => ReadSaveValuesFromFile
end type SaveValues

contains

function ForceOneWord(string)
character(len=*), intent(in) :: string
character(len=len(string) + 100) :: ForceOneWord
character(len=*), parameter :: onewordchars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._"
integer :: i

ForceOneWord = string
do while (.true.)
  do i=1,len(trim(ForceOneWord))
    if (ForceOneWord(i:i) == "/") then
      ForceOneWord = ForceOneWord(:i-1) // "_SLASH_" // ForceOneWord(i+1:)
      exit
    endif
    if (index(onewordchars, ForceOneWord(i:i)).eq.0) then
      print *, "Don't know what to do with character '"//ForceOneWord(i:i)//"'"
      stop 1
    endif
  enddo
  exit
enddo

end function ForceOneWord

function SaveValuesConstructor() result(self)
implicit none
type(SaveValues) :: self
  self%nlogicals = 0
  self%nintegers = 0
  self%nreal8s = 0
  self%ncomplex8s = 0
  self%nstrings = 0

  self%logicalnames(:) = "NONE"
  self%integernames(:) = "NONE"
  self%real8names(:) = "NONE"
  self%complex8names(:) = "NONE"
  self%stringnames(:) = "NONE"
  self%strings(:) = "NONE"
end function SaveValuesConstructor

subroutine savevalue_logical(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*), intent(in) :: name
character(len=len(name)+50) :: name_
logical :: value
  name_ = ForceOneWord(name)
  if (self%nlogicals .ge. maxlogicals) call Error("Too many logicals")
  if (len(trim(name_)) .gt. 100) call Error("Parameter name is too long: "//name_)
  self%nlogicals = self%nlogicals + 1
  self%logicalnames(self%nlogicals) = name_
  self%logicals(self%nlogicals) = value
end subroutine savevalue_logical

subroutine savevalue_integer(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*), intent(in) :: name
character(len=len(name)+50) :: name_
integer :: value
  name_ = ForceOneWord(name)
  if (self%nintegers .ge. maxintegers) call Error("Too many integers")
  if (len(trim(name_)) .gt. 100) call Error("Parameter name is too long: "//name_)
  self%nintegers = self%nintegers + 1
  self%integernames(self%nintegers) = name_
  self%integers(self%nintegers) = value
end subroutine savevalue_integer

subroutine savevalue_real8(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*), intent(in) :: name
character(len=len(name)+50) :: name_
real(8) :: value
  name_ = ForceOneWord(name)
  if (self%nreal8s .ge. maxreal8s) call Error("Too many real8s")
  if (len(trim(name_)) .gt. 100) call Error("Parameter name is too long: "//name_)
  self%nreal8s = self%nreal8s + 1
  self%real8names(self%nreal8s) = name_
  self%real8s(self%nreal8s) = value
end subroutine savevalue_real8

subroutine savevalue_complex8(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*), intent(in) :: name
character(len=len(name)+50) :: name_
complex(8) :: value
  name_ = ForceOneWord(name)
  if (self%ncomplex8s .ge. maxcomplex8s) call Error("Too many complex8s")
  if (len(trim(name_)) .gt. 100) call Error("Parameter name is too long: "//name_)
  self%ncomplex8s = self%ncomplex8s + 1
  self%complex8names(self%ncomplex8s) = name_
  self%complex8s(self%ncomplex8s) = value
end subroutine savevalue_complex8

subroutine savevalue_string(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*), intent(in) :: name, value
character(len=len(name)+50) :: name_
character(len=len(value)+50) :: value_
  name_ = ForceOneWord(name)
  value_ = ForceOneWord(value)
  if (self%nstrings .ge. maxstrings) call Error("Too many strings")
  if (len(trim(name_)) .gt. 100) call Error("Parameter name is too long: "//name_)
  if (len(trim(value_)) .gt. 100) call Error("Parameter value is too long: "//value_)
  self%nstrings = self%nstrings + 1
  self%stringnames(self%nstrings) = name_
  self%strings(self%nstrings) = value_
end subroutine savevalue_string

subroutine SortSaveValues(self)
use modmisc
implicit none
class(SaveValues) :: self
  call BubleSort(self%nlogicals, self%logicalnames(1:self%nlogicals), self%logicals(1:self%nlogicals))
  call BubleSort(self%nintegers, self%integernames(1:self%nintegers), self%integers(1:self%nintegers))
  call BubleSort(self%nreal8s, self%real8names(1:self%nreal8s), self%real8s(1:self%nreal8s))
  call BubleSort(self%ncomplex8s, self%complex8names(1:self%ncomplex8s), self%complex8s(1:self%ncomplex8s))
  call BubleSort(self%nstrings, self%stringnames(1:self%nstrings), self%strings(1:self%nstrings))
end subroutine SortSaveValues

subroutine WriteSaveValuesToFile(self, filename)
use ModParameters
implicit none
class(SaveValues) :: self
character(len=*) :: filename
call self%sort()
open(unit=io_TmpFile, file=filename,form='formatted',status='replace')
write(io_TmpFile,fmt=*) self%nlogicals, self%logicalnames, self%logicals
write(io_TmpFile,fmt=*) self%nintegers, self%integernames, self%integers
write(io_TmpFile,fmt=*) self%nreal8s, self%real8names, self%real8s
write(io_TmpFile,fmt=*) self%ncomplex8s, self%complex8names, self%complex8s
write(io_TmpFile,fmt=*) self%nstrings, self%stringnames, self%strings
close(unit=io_TmpFile)
end subroutine WriteSaveValuesToFile

subroutine ReadSaveValuesFromFile(self, filename)
use ModParameters
implicit none
class(SaveValues) :: self
character(len=*) :: filename
open(unit=io_TmpFile, file=filename,form='formatted',status='old')
read(io_TmpFile,fmt=*) self%nlogicals, self%logicalnames, self%logicals
read(io_TmpFile,fmt=*) self%nintegers, self%integernames, self%integers
read(io_TmpFile,fmt=*) self%nreal8s, self%real8names, self%real8s
read(io_TmpFile,fmt=*) self%ncomplex8s, self%complex8names, self%complex8s
read(io_TmpFile,fmt=*) self%nstrings, self%stringnames, self%strings
close(unit=io_TmpFile)
call self%sort()
end subroutine ReadSaveValuesFromFile

subroutine CompareSaveValues(new, old)
implicit none
class(SaveValues) :: new, old
integer :: i
call new%sort()
call old%sort()

do i=1, max(new%nlogicals, old%nlogicals)
  if (new%logicalnames(i) .gt. old%logicalnames(i) .or. i.gt.old%nlogicals) then
    print *, trim(new%logicalnames(i)), " is set to ", new%logicals(i), ", but was previously not set.  Can't do ReadCSmax."
    stop 1
  else if (new%logicalnames(i) .lt. old%logicalnames(i) .or. i.gt.new%nlogicals) then
    print *, trim(old%logicalnames(i)), " was previously set to ", old%logicals(i), ", but is not set now.  Can't do ReadCSmax"
    stop 1
  else if (new%logicals(i) .neqv. old%logicals(i)) then
    print *, trim(old%logicalnames(i)), " was previously set to ", old%logicals(i), ", but is now set to ", new%logicals(i), ".  Can't do ReadCSmax"
    stop 1
  endif
end do

do i=1, max(new%nintegers, old%nintegers)
  if (new%integernames(i) .gt. old%integernames(i) .or. i.gt.old%nintegers) then
    print *, trim(new%integernames(i)), " is set to ", new%integers(i), ", but was previously not set.  Can't do ReadCSmax."
    stop 1
  else if (new%integernames(i) .lt. old%integernames(i) .or. i.gt.new%nintegers) then
    print *, trim(old%integernames(i)), " was previously set to ", old%integers(i), ", but is not set now.  Can't do ReadCSmax"
    stop 1
  else if (new%integers(i) .ne. old%integers(i)) then
    print *, trim(old%integernames(i)), " was previously set to ", old%integers(i), ", but is now set to ", new%integers(i), ".  Can't do ReadCSmax"
    stop 1
  endif
end do

do i=1, max(new%nreal8s, old%nreal8s)
  if (new%real8names(i) .gt. old%real8names(i) .or. i.gt.old%nreal8s) then
    print *, trim(new%real8names(i)), " is set to ", new%real8s(i), ", but was previously not set.  Can't do ReadCSmax."
    stop 1
  else if (new%real8names(i) .lt. old%real8names(i) .or. i.gt.new%nreal8s) then
    print *, trim(old%real8names(i)), " was previously set to ", old%real8s(i), ", but is not set now.  Can't do ReadCSmax"
    stop 1
  else if (new%real8s(i) .ne. old%real8s(i)) then
    print *, trim(old%real8names(i)), " was previously set to ", old%real8s(i), ", but is now set to ", new%real8s(i), ".  Can't do ReadCSmax"
    stop 1
  endif
end do

do i=1, max(new%ncomplex8s, old%ncomplex8s)
  if (new%complex8names(i) .gt. old%complex8names(i) .or. i.gt.old%ncomplex8s) then
    print *, trim(new%complex8names(i)), " is set to ", new%complex8s(i), ", but was previously not set.  Can't do ReadCSmax."
    stop 1
  else if (new%complex8names(i) .lt. old%complex8names(i) .or. i.gt.new%ncomplex8s) then
    print *, new%complex8names(1:new%ncomplex8s)
    print *, old%complex8names(1:old%ncomplex8s)
    print *, trim(old%complex8names(i)), " was previously set to ", old%complex8s(i), ", but is not set now.  Can't do ReadCSmax"
    stop 1
  else if (new%complex8s(i) .ne. old%complex8s(i)) then
    print *, trim(old%complex8names(i)), " was previously set to ", old%complex8s(i), ", but is now set to ", new%complex8s(i), ".  Can't do ReadCSmax"
    stop 1
  endif
end do

do i=1, max(new%nstrings, old%nstrings)
  if (new%stringnames(i) .gt. old%stringnames(i) .or. i.gt.old%nstrings) then
    print *, trim(new%stringnames(i)), " is set to ", trim(new%strings(i)), ", but was previously not set.  Can't do ReadCSmax."
    stop 1
  else if (new%stringnames(i) .lt. old%stringnames(i) .or. i.gt.new%nstrings) then
    print *, trim(old%stringnames(i)), " was previously set to ", trim(old%strings(i)), ", but is not set now.  Can't do ReadCSmax"
    stop 1
  else if (new%strings(i) .ne. old%strings(i)) then
    print *, trim(old%stringnames(i)), " was previously set to ", trim(old%strings(i)), ", but is now set to ", trim(new%strings(i)), ".  Can't do ReadCSmax"
    stop 1
  endif
end do

end subroutine CompareSaveValues

!ReadCommandLineArgument is overloaded.  Pass the type needed as "dest"
!success is set to true if the argument passed matches argumentname, otherwise it's left alone
!same for success2, success3, success4, success5, and success6 (optional, can be used for other things, see main.F90)
!checkdestchange (optional) modifies how success[2-6] are set by checking if the destination value has changed.
!SetLastArgument (optional) is set to true if the argument matches, otherwise it's set to false
!for examples of all of them see main.F90

subroutine ReadCommandLineArgument_logical(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, checkdestchange, tosave)
implicit none
character(len=*) :: argument, argumentname
logical, intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
logical, optional, intent(in) :: checkdestchange
type(SaveValues), optional :: tosave
character(len=*), parameter :: numbers = "0123456789"
logical :: successval
logical :: dest_store
integer :: length, temp_int

    successval = .true.
    dest_store = dest

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( trim(argument).eq.trim(argumentname) ) then
        dest=.true.
        if (present(checkdestchange)) then
           successval = (.not. checkdestchange .or. dest .neqv. dest_store)
        endif
        if (present(SetLastArgument)) SetLastArgument=successval
        if (present(success2)) success2=success2 .or. successval
        if (present(success3)) success3=success3 .or. successval
        if (present(success4)) success4=success4 .or. successval
        if (present(success5)) success5=success5 .or. successval
        if (present(success6)) success6=success6 .or. successval
        if (successval .and. present(tosave)) call tosave%savevalue_logical(argumentname, dest)
        success = .true.
    elseif( trim(argument).eq."No"//trim(argumentname) ) then
        dest=.false.
        if (present(checkdestchange)) then
           successval = (.not. checkdestchange .or. dest .neqv. dest_store)
        endif
        if (present(SetLastArgument)) SetLastArgument=successval
        if (present(success2)) success2=success2 .or. successval
        if (present(success3)) success3=success3 .or. successval
        if (present(success4)) success4=success4 .or. successval
        if (present(success5)) success5=success5 .or. successval
        if (present(success6)) success6=success6 .or. successval
        if (successval .and. present(tosave)) call tosave%savevalue_logical(argumentname, dest)
        success = .true.
    elseif( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        if( Index(numbers, argument(length+2:length+2)) .ne. 0 ) then
            read(argument(length+2:len(argument)), *) temp_int
            dest = (temp_int.ne.0)
            if (present(checkdestchange)) then
               successval = (.not. checkdestchange .or. dest .neqv. dest_store)
            endif
            if (present(SetLastArgument)) SetLastArgument=successval
            if (present(success2)) success2=success2 .or. successval
            if (present(success3)) success3=success3 .or. successval
            if (present(success4)) success4=success4 .or. successval
            if (present(success5)) success5=success5 .or. successval
            if (present(success6)) success6=success6 .or. successval
            if (successval .and. present(tosave)) call tosave%savevalue_logical(argumentname, dest)
        else
            read(argument(length+2:len(argument)), *) dest
            if (present(checkdestchange)) then
               successval = (.not. checkdestchange .or. dest .neqv. dest_store)
            endif
            if (present(SetLastArgument)) SetLastArgument=successval
            if (present(success2)) success2=success2 .or. successval
            if (present(success3)) success3=success3 .or. successval
            if (present(success4)) success4=success4 .or. successval
            if (present(success5)) success5=success5 .or. successval
            if (present(success6)) success6=success6 .or. successval
            if (successval .and. present(tosave)) call tosave%savevalue_logical(argumentname, dest)
        endif
        success = .true.
    endif

end subroutine ReadCommandLineArgument_logical


subroutine ReadCommandLineArgument_integer(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, checkdestchange, multiply, tosave)
implicit none
character(len=*) :: argument, argumentname
integer, intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
logical, optional, intent(in) :: checkdestchange
type(SaveValues), optional :: tosave
integer, optional, intent(in) :: multiply
logical :: successval
integer :: length, dest_store

    successval = .true.
    dest_store = dest

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        read(argument(length+2:len(argument)), *) dest
        if (present(multiply)) dest = dest*multiply
        ! Checkdestchange after all multiplications are done!
        if (present(checkdestchange)) then
           successval = (.not. checkdestchange .or. dest.ne.dest_store)
        endif
        if (present(SetLastArgument)) SetLastArgument=successval
        if (present(success2)) success2=success2 .or. successval
        if (present(success3)) success3=success3 .or. successval
        if (present(success4)) success4=success4 .or. successval
        if (present(success5)) success5=success5 .or. successval
        if (present(success6)) success6=success6 .or. successval
        if (successval .and. present(tosave)) call tosave%savevalue_integer(argumentname, dest)
        success = .true.
    endif

end subroutine ReadCommandLineArgument_integer


subroutine ReadCommandLineArgument_real8(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, checkdestchange, multiply, tosave)
implicit none
character(len=*) :: argument, argumentname
real(8), intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
logical, optional, intent(in) :: checkdestchange
type(SaveValues), optional :: tosave
real(8), optional, intent(in) :: multiply
logical :: successval
integer :: length
real(8) :: dest_store

    successval = .true.
    dest_store = dest

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        read(argument(length+2:len(argument)), *) dest
        if (present(multiply)) dest = dest*multiply
        ! Checkdestchange after all multiplications are done!
        if (present(checkdestchange)) then
           successval = (.not. checkdestchange .or. dest.ne.dest_store)
        endif
        if (present(SetLastArgument)) SetLastArgument=successval
        if (present(success2)) success2=success2 .or. successval
        if (present(success3)) success3=success3 .or. successval
        if (present(success4)) success4=success4 .or. successval
        if (present(success5)) success5=success5 .or. successval
        if (present(success6)) success6=success6 .or. successval
        if (successval .and. present(tosave)) call tosave%savevalue_real8(argumentname, dest)
        success = .true.
    endif

end subroutine ReadCommandLineArgument_real8


subroutine ReadCommandLineArgument_complex8(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, checkdestchange, multiply, multiplyreal, tosave)
implicit none
character(len=*) :: argument, argumentname
complex(8), intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
logical, optional, intent(in) :: checkdestchange
type(SaveValues), optional :: tosave
complex(8), optional, intent(in) :: multiply
real(8), optional, intent(in) :: multiplyreal
logical :: successval
integer :: length
real(8) :: re, im
complex(8) :: dest_store

    successval = .true.
    dest_store = dest

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        if( Index(argument(length+2:len(trim(argument))),",").eq.0 &
       .and. Index(argument(length+2:len(trim(argument)))," ").eq.0 ) then
            print *, "Argument ", argumentname, " is complex."
            print *, "Example syntax for complex arguments:"
            print *, "      ", argumentname, "=1,2"
            print *, "   or ", argumentname, "=1.0d0,2.0d0"
            print *, "for 1+2i"
            stop 1
        endif
        read(argument(length+2:len(argument)), *) re, im
        dest = dcmplx(re, im)
        if (present(multiply)) dest = dest*multiply
        if (present(multiplyreal)) dest = dest*multiplyreal
        ! Checkdestchange after all multiplications are done!
        if (present(checkdestchange)) then
           !print *, argumentname,": checkdestchange, dest_store, dest",checkdestchange,dest_store,dest
           successval = (.not. checkdestchange .or. dest.ne.dest_store)
        endif
        if (present(SetLastArgument)) SetLastArgument=successval
        if (present(success2)) success2=success2 .or. successval
        if (present(success3)) success3=success3 .or. successval
        if (present(success4)) success4=success4 .or. successval
        if (present(success5)) success5=success5 .or. successval
        if (present(success6)) success6=success6 .or. successval
        if (successval .and. present(tosave)) call tosave%savevalue_complex8(argumentname, dest)
        success = .true.
    endif

end subroutine ReadCommandLineArgument_complex8


subroutine ReadCommandLineArgument_string(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, checkdestchange, tosave)
implicit none
character(len=*) :: argument, argumentname
character(len=*), intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
logical, optional, intent(in) :: checkdestchange
type(SaveValues), optional :: tosave
logical :: successval
integer :: length
character(len=len(dest)) :: dest_store

    successval = .true.
    dest_store = dest

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        if( len(dest).lt.len(trim(argument))-(length+1) ) then
            print "(A,A,A,I4,A)", "Argument ", argument, " is too long!  Maximum allowed length is ", len(dest), " characters."
            stop 1
        endif
        dest = argument(length+2:len(argument))
        if (present(checkdestchange)) then
           successval = (.not. checkdestchange .or. trim(dest)==trim(dest_store))
        endif
        if (present(SetLastArgument)) SetLastArgument=successval
        if (present(success2)) success2=success2 .or. successval
        if (present(success3)) success3=success3 .or. successval
        if (present(success4)) success4=success4 .or. successval
        if (present(success5)) success5=success5 .or. successval
        if (present(success6)) success6=success6 .or. successval
        if (successval .and. present(tosave)) call tosave%savevalue_string(argumentname, dest)
        success = .true.
    endif

end subroutine ReadCommandLineArgument_string

END MODULE
