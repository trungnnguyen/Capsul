module typeThermoExpa
  use utils, only : RKIND, IKIND, LKIND
  implicit none

  integer(kind = IKIND), public, parameter :: kIsoThExpa    = 1
  integer(kind = IKIND), public, parameter :: kAnisoThExpa  = 2

  type, abstract, public :: ThermoExpa
    real(kind = RKIND)  :: fTempRef
    real(kind = RKIND)  :: fAlpha(3)
    real(kind = RKIND), allocatable :: fMatProp(:, :)

    logical(kind = LKIND) :: fIsCoefTempDep
  contains
    procedure, public :: GetDefGrdTh
    procedure(UpdExpaCoefIntf), public, deferred :: UpdExpaCoef
  end type ThermoExpa

  abstract interface
    subroutine UpdExpaCoefIntf(this, tempCur)
      import ThermoExpa
      class(ThermoExpa) :: this
      real(kind = 8), intent(in) :: tempCur
    end subroutine UpdExpaCoefIntf
  end interface

contains

  function GetDefGrdTh(this, tempCur) result(defGrdTh)
    class(ThermoExpa), intent(inout) :: this
    real(kind = RKIND), intent(in) :: tempCur
    real(kind = RKIND) :: defGrdTh(3, 3)

    integer(kind = IKIND) :: i

    defGrdTh = 0.0d0
    if (this%fIsCoefTempDep) call this%UpdExpaCoef(tempCur)

    do i = 1, 3
      defGrdTh(i, i) = 1 + this%fAlpha(i)*(tempCur - this%fTempRef)
    end do

  end function GetDefGrdTh

end module typeThermoExpa


module typeThermoExpaIso
  use typeThermoExpa
  use utils, only : RKIND, IKIND

  implicit none


  type, extends(ThermoExpa), public :: ThermoExpaIso
  contains
    procedure, public :: UpdExpaCoef
    final             :: Destructor
  end type ThermoExpaIso

  interface ThermoExpaIso
    module procedure Constructor
  end interface

contains

  function Constructor(fID, polyOrder, tempRef, tempInit) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    integer(kind = IKIND), intent(in) :: polyOrder
    real(kind = RKIND),    intent(in) :: tempRef
    real(kind = RKIND),    intent(in) :: tempInit
    type(ThermoExpaIso), pointer      :: this

    character(len = 256) :: line
    integer(kind = IKIND) :: i 


    allocate(this)

    this%fTempRef = tempRef
    if (polyOrder == 0) then
      this%fIsCoefTempDep = .false.
    else
      this%fIsCoefTempDep = .true.
    end if

    allocate(this%fMatProp(polyOrder+1, 1))
    do i = 1, polyOrder + 1
      call ReadLine(fID, line)
      read(line, *) this%fMatProp(i, 1)
    end do

    call this%UpdExpaCoef(tempInit)

  end function Constructor


  subroutine Destructor(this)
    type(ThermoExpaIso) :: this

    if (allocated(this%fMatProp)) deallocate(this%fMatProp)

  end subroutine


  subroutine UpdExpaCoef(this, tempCur)
    use algebra, only : Polynomial
    class(ThermoExpaIso) :: this
    real(kind = RKIND), intent(in) :: tempCur

    this%fAlpha(1) = Polynomial(this%fMatProp(:, 1), tempCur)
    this%fAlpha(2) = this%fAlpha(1)
    this%fAlpha(3) = this%fAlpha(1)

  end subroutine UpdExpaCoef


end module typeThermoExpaIso



module typeThermoExpaAniso
  use typeThermoExpa
  use utils, only : RKIND, IKIND

  implicit none

  type, extends(ThermoExpa), public :: ThermoExpaAniso
  contains
    procedure, public :: UpdExpaCoef
    final             :: Destructor
  end type ThermoExpaAniso

  interface ThermoExpaAniso
    module procedure Constructor
  end interface

contains

  function Constructor(fID, polyOrder, tempRef, tempInit) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    integer(kind = IKIND), intent(in) :: polyOrder
    real(kind = RKIND),    intent(in) :: tempRef
    real(kind = RKIND),    intent(in) :: tempInit
    type(ThermoExpaAniso), pointer    :: this

    character(len = 256) :: line
    integer(kind = IKIND) :: i, j


    allocate(this)

    this%fTempRef = tempRef

    if (polyOrder == 0) then
      this%fIsCoefTempDep = .false.
    else
      this%fIsCoefTempDep = .true.
    end if

    allocate(this%fMatProp(polyOrder+1, 3))
    do i = 1, polyOrder + 1
      call ReadLine(fID, line)
      read(line, *) (this%fMatProp(i, j), j = 1, 3)
    end do

    call this%UpdExpaCoef(tempInit)

  end function Constructor


  subroutine Destructor(this)
    type(ThermoExpaAniso) :: this

    if (allocated(this%fMatProp)) deallocate(this%fMatProp)

  end subroutine


  subroutine UpdExpaCoef(this, tempCur)
    use algebra, only : Polynomial
    class(ThermoExpaAniso) :: this
    real(kind = RKIND), intent(in) :: tempCur

    this%fAlpha(1) = Polynomial(this%fMatProp(:, 1), tempCur)
    this%fAlpha(2) = Polynomial(this%fMatProp(:, 2), tempCur)
    this%fAlpha(3) = Polynomial(this%fMatProp(:, 3), tempCur)

  end subroutine UpdExpaCoef

      

end module typeThermoExpaAniso

