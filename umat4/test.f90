module modTypeA
  implicit none
 
  type, public :: typeA
    real(kind = 8)    :: dblA
    integer(kind = 4) :: intA
    real(kind = 8), allocatable :: arrA(:)
  contains
    procedure, public :: printA
    procedure, public :: initA
  end type typeA

  interface typeA
    module procedure Constructor
  end interface

contains

  subroutine initA(this, intVal, dblVal)
    class(typeA) :: this
    integer(kind = 4) :: intVal
    real(kind = 8) :: dblVal

    this%dblA = dblVal
    this%intA = intVal

    allocate(this%arrA(this%intA))
    this%arrA = this%dblA

  end subroutine


  function Constructor(intVal, dblVal) result(this)
    integer(kind = 4), intent(in) :: intVal
    real(kind = 8),    intent(in) :: dblVal
    class(typeA),      pointer :: this

    write(*, *) "Constructor is called!"
    allocate(this)

    this%intA = intVal
    this%dblA = dblVal
    allocate(this%arrA(this%intA))
    this%arrA = this%dblA

  end function

  subroutine printA(this)
    class(typeA), intent(in) :: this

    write(*, *) "value of array A"
    write(*, *) this%intA, this%arrA

  end subroutine

end module


module data
  use modTypeA
  class(typeA), pointer :: pA1 => null()
  class(typeA), pointer :: pA2 => null()

  contains
    subroutine init1(intA, dblA)
      integer(kind = 4) :: intA
      real(kind = 8)    :: dblA

      pA1 => typeA(intA, dblA)
    end subroutine


    subroutine init2(dblA)
      real(kind = 8)    :: dblA

      pA2 => pA1
      pA2%intA = pA2%intA + 1
      pA2%dblA = dblA + 2.0d0
      pA2%arrA = pA2%dblA

    end subroutine

end module data


program test
  use modTypeA
  use data, only : pA1, pA2, init1, init2
    
  implicit none

  integer(kind = 4)     :: testInt
  type(typeA),  target  :: theA,theB

  !theA = typeA(4, 1.12d0)
  call init1(4, 1.12d0)
  call init2(1.53d0)

  pA1%intA = pA1%intA + 2
  call pA1%printA
  call pA2%printA
  
end program
