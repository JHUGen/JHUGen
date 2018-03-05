MODULE ModCommandLine
implicit none

interface ReadCommandLineArgument
    module procedure ReadCommandLineArgument_logical, ReadCommandLineArgument_integer, ReadCommandLineArgument_real8,&
                     ReadCommandLineArgument_complex8, ReadCommandLineArgument_string
end interface

interface savevalue
  module procedure :: savevalue_logical, savevalue_integer, savevalue_real8, savevalue_complex8, savevalue_string
end interface savevalue

type SaveValues
  character(len=100) :: logicalnames(1:100), integernames(1:100), real8names(1:100), complex8names(1:100), stringnames(1:100)
  integer :: nlogicals, nintegers, nreal8s, ncomplex8s, nstrings
  logical :: logicals(1:100)
  integer :: integers(1:100)
  real(8) :: real8s(1:100)
  complex(8) :: complex8s(1:100)
  character(len=100) :: strings(1:100)
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

interface SaveValues
  procedure :: SaveValuesConstructor
end interface SaveValues

contains

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
character(len=*) :: name
logical :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%nlogicals = self%nlogicals + 1
  self%logicalnames(self%nlogicals) = name
  self%logicals(self%nlogicals) = value
end subroutine savevalue_logical

subroutine savevalue_integer(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
integer :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%nintegers = self%nintegers + 1
  self%integernames(self%nintegers) = name
  self%integers(self%nintegers) = value
end subroutine savevalue_integer

subroutine savevalue_real8(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
real(8) :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%nreal8s = self%nreal8s + 1
  self%real8names(self%nreal8s) = name
  self%real8s(self%nreal8s) = value
end subroutine savevalue_real8

subroutine savevalue_complex8(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
complex(8) :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%ncomplex8s = self%ncomplex8s + 1
  self%complex8names(self%ncomplex8s) = name
  self%complex8s(self%ncomplex8s) = value
end subroutine savevalue_complex8

subroutine savevalue_string(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
character(len=*) :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  if (len(value) .gt. 100) call Error("Parameter value is too long: "//value)
  self%nstrings = self%nstrings + 1
  self%stringnames(self%nstrings) = name
  self%strings(self%nstrings) = value
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
!SetLastArgument (optional) is set to true if the argument matches, otherwise it's set to false
!for examples of all of them see main.F90

subroutine ReadCommandLineArgument_logical(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, tosave)
implicit none
character(len=*) :: argument, argumentname
logical, intent(inout) :: dest
logical, intent(inout) :: success
integer :: length
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
type(SaveValues), optional :: tosave
integer :: temp_int
character(len=*), parameter :: numbers = "0123456789"

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( trim(argument).eq.trim(argumentname) ) then
        dest=.true.
        success=.true.
        if (present(SetLastArgument)) SetLastArgument=.true.
        if (present(success2)) success2=.true.
        if (present(success3)) success3=.true.
        if (present(success4)) success4=.true.
        if (present(success5)) success5=.true.
        if (present(success6)) success6 = .true.
        if (present(tosave)) call tosave%savevalue_logical(argumentname, dest)
    elseif( trim(argument).eq."No"//trim(argumentname) ) then
        dest=.false.
        success=.true.
        if (present(SetLastArgument)) SetLastArgument=.true.
        if (present(success2)) success2=.true.
        if (present(success3)) success3=.true.
        if (present(success4)) success4=.true.
        if (present(success5)) success5=.true.
        if (present(success6)) success6 = .true.
        if (present(tosave)) call tosave%savevalue_logical(argumentname, dest)
    elseif( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        if( Index(numbers, argument(length+2:length+2)) .ne. 0 ) then
            read(argument(length+2:len(argument)), *) temp_int
            dest = (temp_int.ne.0)
            success=.true.
            if (present(SetLastArgument)) SetLastArgument=.true.
            if (present(success2)) success2=.true.
            if (present(success3)) success3=.true.
            if (present(success4)) success4=.true.
            if (present(success5)) success5=.true.
            if (present(success6)) success6 = .true.
            if (present(tosave)) call tosave%savevalue_logical(argumentname, dest)
        else
            read(argument(length+2:len(argument)), *) dest
            success=.true.
            if (present(SetLastArgument)) SetLastArgument=.true.
            if (present(success2)) success2=.true.
            if (present(success3)) success3=.true.
            if (present(success4)) success4=.true.
            if (present(success5)) success5=.true.
            if (present(success6)) success6 = .true.
            if (present(tosave)) call tosave%savevalue_logical(argumentname, dest)
        endif
    endif

end subroutine ReadCommandLineArgument_logical


subroutine ReadCommandLineArgument_integer(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, multiply, tosave)
implicit none
character(len=*) :: argument, argumentname
integer, intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
type(SaveValues), optional :: tosave
integer, optional, intent(in) :: multiply
integer :: length

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        read(argument(length+2:len(argument)), *) dest
        if (present(multiply)) dest = dest*multiply
        success=.true.
        if (present(SetLastArgument)) SetLastArgument=.true.
        if (present(success2)) success2=.true.
        if (present(success3)) success3=.true.
        if (present(success4)) success4=.true.
        if (present(success5)) success5=.true.
        if (present(success6)) success6 = .true.
        if (present(tosave)) call tosave%savevalue_integer(argumentname, dest)
    endif

end subroutine ReadCommandLineArgument_integer


subroutine ReadCommandLineArgument_real8(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, multiply, tosave)
implicit none
character(len=*) :: argument, argumentname
real(8), intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
type(SaveValues), optional :: tosave
real(8), optional, intent(in) :: multiply
integer :: length

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        read(argument(length+2:len(argument)), *) dest
        if (present(multiply)) dest = dest*multiply
        success=.true.
        if (present(SetLastArgument)) SetLastArgument=.true.
        if (present(success2)) success2=.true.
        if (present(success3)) success3=.true.
        if (present(success4)) success4=.true.
        if (present(success5)) success5=.true.
        if (present(success6)) success6 = .true.
        if (present(tosave)) call tosave%savevalue_real8(argumentname, dest)
    endif

end subroutine ReadCommandLineArgument_real8


subroutine ReadCommandLineArgument_complex8(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, multiply, multiplyreal, tosave)
implicit none
character(len=*) :: argument, argumentname
complex(8), intent(inout) :: dest
real(8) :: re, im
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
type(SaveValues), optional :: tosave
complex(8), optional, intent(in) :: multiply
real(8), optional, intent(in) :: multiplyreal
integer :: length

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
        success=.true.
        if (present(SetLastArgument)) SetLastArgument=.true.
        if (present(success2)) success2=.true.
        if (present(success3)) success3=.true.
        if (present(success4)) success4=.true.
        if (present(success5)) success5=.true.
        if (present(success6)) success6 = .true.
        if (present(tosave)) call tosave%savevalue_complex8(argumentname, dest)
    endif

end subroutine ReadCommandLineArgument_complex8


subroutine ReadCommandLineArgument_string(argument, argumentname, success, dest, SetLastArgument, success2, success3, success4, success5, success6, tosave)
implicit none
character(len=*) :: argument, argumentname
character(len=*), intent(inout) :: dest
logical, intent(inout) :: success
logical, optional, intent(inout) :: SetLastArgument, success2, success3, success4, success5, success6
type(SaveValues), optional :: tosave
integer :: length

    if (present(SetLastArgument)) SetLastArgument=.false.

    length=len(trim(argumentname))

    if( argument(1:length+1) .eq. trim(argumentname)//"=" ) then
        if( len(dest).lt.len(trim(argument))-(length+1) ) then
            print "(A,A,A,I4,A)", "Argument ", argument, " is too long!  Maximum allowed length is ", len(dest), " characters."
            stop 1
        endif
        dest = argument(length+2:len(argument))
        success=.true.
        if (present(SetLastArgument)) SetLastArgument=.true.
        if (present(success2)) success2=.true.
        if (present(success3)) success3=.true.
        if (present(success4)) success4=.true.
        if (present(success5)) success5=.true.
        if (present(success6)) success6 = .true.
        if (present(tosave)) call tosave%savevalue_string(argumentname, dest)
    endif

end subroutine ReadCommandLineArgument_string

END MODULE
