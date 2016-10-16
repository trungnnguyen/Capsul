module typeCrysElas
  use utils, only : RKIND, IKIND, LKIND
  implicit none

  integer(kind = RKIND), public, parameter :: kIso   = 1
  integer(kind = RKIND), public, parameter :: kCubic = 2
  integer(kind = RKIND), public, parameter :: kAniso = 3

  type, abstract, public :: CrysElas
    real(kind = RKIND) :: fModuli(3, 3, 3, 3)
    real(kind = RKIND) :: fC11
    real(kind = RKIND) :: fC12
    real(kind = RKIND) :: fC44
    real(kind = RKIND) :: fC33
    real(kind = RKIND) :: fC13
    real(kind = RKIND) :: fC66

    real(kind = RKIND), allocatable :: fMatProp(:, :)

    logical(kind = LKIND) :: fIsTempDep
  contains
    ! public subroutines and functions
    procedure, public :: ElasticModuli

    ! protected subroutines and functions
    procedure, public :: FormModuli
    procedure(UpdModuliIntf), public, deferred :: UpdModuli

  end type CrysElas

  abstract interface
    subroutine UpdModuliIntf(this, temp)
      import CrysElas
      integer(kind = 4), parameter :: RKIND = 8
      class(CrysElas)   :: this
      real(kind = RKIND), intent(in) :: temp
    end subroutine UpdModuliIntf
  end interface

contains

  function ElasticModuli(this, tempCur) result(moduli)
    class(CrysElas),    intent(in) :: this
    real(kind = RKIND), intent(in) :: tempCur
    real(kind = RKIND) :: moduli(3, 3, 3, 3)


    if (this%fIsTempDep) call this%UpdModuli(tempCur)

    moduli = this%fModuli

  end function ElasticModuli



  subroutine FormModuli(this)
    class(CrysElas) :: this

    this%fModuli = 0.0d0

    this%fModuli(1, 1, 1, 1) = this%fC11
    this%fModuli(2, 2, 2, 2) = this%fC11
    this%fModuli(3, 3, 3, 3) = this%fC33

    this%fModuli(1, 3, 1, 3) = this%fC44
    this%fModuli(3, 1, 1, 3) = this%fC44
    this%fModuli(1, 3, 3, 1) = this%fC44
    this%fModuli(3, 1, 3, 1) = this%fC44

    this%fModuli(2, 3, 2, 3) = this%fC44
    this%fModuli(3, 2, 2, 3) = this%fC44
    this%fModuli(2, 3, 3, 2) = this%fC44
    this%fModuli(3, 2, 3, 2) = this%fC44

    this%fModuli(1, 2, 1, 2) = this%fC66
    this%fModuli(2, 1, 1, 2) = this%fC66
    this%fModuli(1, 2, 2, 1) = this%fC66
    this%fModuli(2, 1, 2, 1) = this%fC66

    this%fModuli(1, 1, 2, 2) = this%fC12
    this%fModuli(2, 2, 1, 1) = this%fC12

    this%fModuli(1, 1, 3, 3) = this%fC13
    this%fModuli(3, 3, 1, 1) = this%fC13

    this%fModuli(2, 2, 3, 3) = this%fC13
    this%fModuli(3, 3, 2, 2) = this%fC13

  end subroutine


end module typeCrysElas



module typeCrysElasIso
  use typeCrysElas
  use utils, only : RKIND, IKIND
  implicit none

  type, extends(CrysElas), public :: CrysElasIso
    real(kind = RKIND) :: fYoung
    real(kind = RKIND) :: fPoisson
  contains
    procedure, public :: UpdModuli
    final :: Destructor
  end type

  interface CrysElasIso
    module procedure Constructor
  end interface

contains

  function Constructor(fID, tempDep, polyOrder, tempInit) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    logical(kind = LKIND), intent(in) :: tempDep
    integer(kind = IKIND), intent(in) :: polyOrder
    real(kind = RKIND),    intent(in) :: tempInit
    type(CrysElasIso), pointer :: this
    
    character(len = 256) :: line
    integer(kind = IKIND) :: i, j


    allocate(this)

    this%fIsTempDep = tempDep

    if (.not. this%fIsTempDep .and. polyOrder > 0) then
      write(*, *) "Temperature-depended elastic properties specified for a temperature-independed problem!" 
      write(*, *) "Only the constant term of the polynomial are used"
    end if

   
    allocate(this%fMatProp(polyOrder+1, 2))
    do i = 1, polyOrder+1
      call ReadLine(fID, line)
      read(line, *) (this%fMatProp(i, j), j = 1, 2) 
    end do

    call this%UpdModuli(tempInit)

  end function Constructor



  subroutine Destructor(this)
    type(CrysElasIso) :: this

    if (allocated(this%fMatProp)) deallocate(this%fMatProp)

  end subroutine Destructor



  subroutine UpdModuli(this, temp)
    use algebra, only : Polynomial

    class(CrysElasIso) :: this
    real(kind = RKIND), intent(in) :: temp
    

    if (this%fIsTempDep) then
      this%fYoung   = Polynomial(this%fMatProp(:, 1), temp)
      this%fPoisson = Polynomial(this%fMatProp(:, 2), temp)
    else
      this%fYoung   = this%fMatProp(1, 1)
      this%fPoisson = this%fMatProp(1, 2) 
    end if

    this%fC44 = 0.5d0*this%fYoung/(1.0 + this%fPoisson)                       ! mu
    this%fC12 = 2.0d0*this%fC44*this%fPoisson/(1.0d0 - 2.0d0*this%fPoisson)   ! lamda
    this%fC11 = this%fC12 + 2.0d0*this%fC44
    this%fC33 = this%fC11
    this%fC13 = this%fC12
    this%fC66 = this%fC44

    call this%FormModuli()

  end subroutine


end module typeCrysElasIso



module typeCrysElasCubic
  use typeCrysElas
  use utils, only : RKIND, IKIND
  implicit none

  type, extends(CrysElas), public :: CrysElasCubic
  contains
    procedure, public :: UpdModuli
    final             :: Destructor
  end type

  interface CrysElasCubic
    module procedure Constructor
  end interface

contains

  function Constructor(fID, tempDep, polyOrder, temp) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    logical(kind = LKIND), intent(in) :: tempDep
    integer(kind = IKIND), intent(in) :: polyOrder
    real(kind = RKIND),    intent(in) :: temp
    type(CrysElasCubic), pointer      :: this
    
    character(len = 256)  :: line
    integer(kind = IKIND) :: i, j


    allocate(this)

    this%fIsTempDep = tempDep
    
    if (.not. this%fIsTempDep .and. polyOrder > 0) then
      write(*, *) "Temperature-dependent elastic properties specified for a temperature-independent problem!" 
      write(*, *) "Only the constant term of the polynomial are used"
    end if

    allocate(this%fMatProp(polyOrder+1, 3))
    do i = 1, polyOrder+1
      call ReadLine(fID, line)
      read(line, *) (this%fMatProp(i, j), j = 1, 3)
    end do


    call this%UpdModuli(temp)
    
  end function Constructor


  subroutine Destructor(this)
    type(CrysElasCubic) :: this

    if (allocated(this%fMatProp)) deallocate(this%fMatProp)

  end subroutine Destructor



  subroutine UpdModuli(this, temp)
    use algebra, only : Polynomial

    class(CrysElasCubic) :: this
    real(kind = RKIND), intent(in) :: temp
    
    if (this%fIsTempDep) then
      this%fC11 = Polynomial(this%fMatProp(:, 1), temp)
      this%fC12 = Polynomial(this%fMatProp(:, 2), temp)
      this%fC44 = Polynomial(this%fMatProp(:, 3), temp)
    else
      this%fC11 = this%fMatProp(1, 1)
      this%fC12 = this%fMatProp(1, 2)
      this%fC44 = this%fMatProp(1, 3)
    end if

    this%fC33 = this%fC11
    this%fC13 = this%fC12
    this%fC66 = this%fC44

    call this%FormModuli()

  end subroutine

end module typeCrysElasCubic



module typeCrysElasAniso
  use typeCrysElas
  use utils, only : RKIND, IKIND
  implicit none

  type, extends(CrysElas), public :: CrysElasAniso
  contains
    procedure, public :: UpdModuli
    final :: Destructor
  end type

  interface CrysElasAniso
    module procedure Constructor
  end interface

contains

  function Constructor(fID, tempDep, polyOrder, temp) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    logical(kind = LKIND), intent(in) :: tempDep
    integer(kind = IKIND), intent(in) :: polyOrder
    real(kind = RKIND),    intent(in) :: temp
    type(CrysElasAniso),   pointer    :: this

    character(len = 256) :: line
    integer(kind = IKIND) :: i, j

    allocate(this)

    this%fIsTempDep = tempDep

    if (.not. this%fIsTempDep .and. polyOrder > 0) then
      write(*, *) "Temperature-depended elastic properties specified for a temperature-independed problem!" 
      write(*, *) "Only the constant term of the polynomial are used"
    end if

    allocate(this%fMatProp(polyOrder+1, 6))
    do i = 1, polyOrder+1
      call ReadLine(fID, line)
      read(line, *) (this%fMatProp(i, j), j = 1, 6)
    end do
    
    call this%UpdModuli(temp)
    
  end function Constructor


  subroutine Destructor(this)
    type(CrysElasAniso) :: this

    if (allocated(this%fMatProp)) deallocate(this%fMatProp)

  end subroutine Destructor



  subroutine UpdModuli(this, temp)
    use algebra, only : Polynomial

    class(CrysElasAniso) :: this
    real(kind = RKIND), intent(in) :: temp
   

    if (this%fIsTempDep) then 
      this%fC11 = Polynomial(this%fMatProp(:, 1), temp)
      this%fC12 = Polynomial(this%fMatProp(:, 2), temp)
      this%fC44 = Polynomial(this%fMatProp(:, 3), temp)
      this%fC33 = Polynomial(this%fMatProp(:, 4), temp)
      this%fC13 = Polynomial(this%fMatProp(:, 5), temp)
      this%fC66 = Polynomial(this%fMatProp(:, 6), temp)
    else
      this%fC11 = this%fMatProp(1, 1)
      this%fC12 = this%fMatProp(1, 2)
      this%fC44 = this%fMatProp(1, 3)                

      this%fC33 = this%fMatProp(1, 4)
      this%fC13 = this%fMatProp(1, 5)
      this%fC66 = this%fMatProp(1, 6)
    end if

    call this%FormModuli()

  end subroutine



end module typeCrysElasAniso

