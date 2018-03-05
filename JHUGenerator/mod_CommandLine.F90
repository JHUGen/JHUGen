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
  integer :: logicalupto, integerupto, real8upto, complex8upto, stringupto
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
end type SaveValues

interface SaveValues
  procedure :: SaveValuesConstructor
end interface SaveValues

contains

function SaveValuesConstructor() result(self)
implicit none
type(SaveValues) :: self
  self%logicalupto = 0
  self%integerupto = 0
  self%real8upto = 0
  self%complex8upto = 0
  self%stringupto = 0

  self%logicalnames(:) = ""
  self%integernames(:) = ""
  self%real8names(:) = ""
  self%complex8names(:) = ""
  self%stringnames(:) = ""
end function SaveValuesConstructor

subroutine savevalue_logical(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
logical :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%logicalnames(self%logicalupto) = name
  self%logicals(self%logicalupto) = value
  self%logicalupto = self%logicalupto + 1
end subroutine savevalue_logical

subroutine savevalue_integer(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
integer :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%integernames(self%integerupto) = name
  self%integers(self%integerupto) = value
  self%integerupto = self%integerupto + 1
end subroutine savevalue_integer

subroutine savevalue_real8(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
real(8) :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%real8names(self%real8upto) = name
  self%real8s(self%real8upto) = value
  self%real8upto = self%real8upto + 1
end subroutine savevalue_real8

subroutine savevalue_complex8(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
complex(8) :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  self%complex8names(self%complex8upto) = name
  self%complex8s(self%complex8upto) = value
  self%complex8upto = self%complex8upto + 1
end subroutine savevalue_complex8

subroutine savevalue_string(self, name, value)
use ModMisc
implicit none
class(SaveValues) :: self
character(len=*) :: name
character(len=*) :: value
  if (len(name) .gt. 100) call Error("Parameter name is too long: "//name)
  if (len(value) .gt. 100) call Error("Parameter value is too long: "//value)
  self%stringnames(self%stringupto) = name
  self%strings(self%stringupto) = value
  self%stringupto = self%stringupto + 1
end subroutine savevalue_string

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
