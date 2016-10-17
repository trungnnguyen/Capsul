  module modSlipGeom
    use utils, only : RKIND, IKIND
    implicit none
    private

    !********* Declaration of base class SlipGeom
    type, abstract, public :: SlipGeom
      integer(kind = IKIND) :: fNumSlip
      real(kind = RKIND), allocatable:: fSlipVecS(:, :)
      real(kind = RKIND), allocatable:: fSlipVecM(:, :)
      real(kind = RKIND), allocatable:: fSchmidt(:, :, :)

    contains
      procedure, public, non_overridable :: GetNumSlip
      procedure, public, non_overridable :: GetSchmidt

      procedure, private :: InitSlipGeom
      procedure, private :: DestroySlipGeom
      procedure, private :: SetSchmidt

      procedure(SetSlipVectorsInterface), public, deferred :: SetSlipVectors
    end type SlipGeom

    abstract interface
      subroutine SetSlipVectorsInterface(this)
        import SlipGeom
        class(SlipGeom) :: this
      end subroutine SetSlipVectorsInterface
    end interface
    !**********End of declaration of class SlipGeom


    !**********Declaration of derived class FCCGeom
    type, extends(SlipGeom), public :: FCCGeom
    contains
      procedure, public :: SetSlipVectors => SetSlipVecsFCC
      final :: FCCGeomDestructor
    end type FCCGeom

    interface FCCGeom
      module procedure FCCGeomConstructor
    end interface
    !********** End of declaration of class FCCGeom

   
    !********Declaration of derived class BCCGeom
    type, extends(SlipGeom), public :: BCCGeom
    contains
      procedure, public :: SetSlipVectors => SetSlipVecsBCC
      final :: BCCGeomDestructor
    end type BCCGeom

    interface BCCGeom
      module procedure BCCGeomConstructor
    end interface
    !******* End of declaration of class BCCGeom

    type, extends(SlipGeom), public :: UserGeom
      integer(kind = IKIND) :: fID  
    contains
      procedure, public :: SetSlipVectors => SetSlipVecsUser
      final :: UserGeomDestructor
    end type UserGeom

    interface UserGeom
      module procedure UserGeomConstructor
    end interface

  contains
          
    subroutine InitSlipGeom(this, numSlip)
      class(SlipGeom), intent(inout) :: this
      integer(kind = IKIND), intent(in) :: NumSlip


      write(*, *) "InitSlipGeom is called"
      this%fNumSlip = numSlip

      allocate(this%fSlipVecS(3, this%fNumSlip))
      allocate(this%fSlipVecM(3, this%fNumSlip))
      allocate(this%fSchmidt(3, 3, this%fNumSlip))

      call this%SetSlipVectors()
      call this%SetSchmidt()

    end subroutine InitSlipGeom


    subroutine DestroySlipGeom(this)
      class(SlipGeom) :: this

      if (allocated(this%fSchmidt)) then
        deallocate(this%fSchmidt)
      end if

    end subroutine DestroySlipGeom


    function GetNumSlip(this) result(numSlip) 
      integer(kind = IKIND) :: numSlip
      class(SlipGeom), intent(in) :: this

      numSlip = this%fNumSlip
    
    end function GetNumSlip


    subroutine SetSchmidt(this)
      use algebra, only : normalize, out_product
      class(SlipGeom) :: this
      integer :: i


      do i = 1, this%fNumSlip
        call normalize(this%fSlipVecS(:, i))
        call normalize(this%fSlipVecM(:, i))
        this%fSchmidt(:, :, i) = out_product(this%fSlipVecS(:, i), this%fSlipVecM(:, i))
      end do

    end subroutine SetSchmidt


    function GetSchmidt(this, iSlip) result(schmidt)
      class(SlipGeom), intent(in) :: this
      integer(kind = IKIND), intent(in) :: iSlip
      real(kind = RKIND) :: schmidt(3, 3)

      schmidt = this%fSchmidt(:, :, iSlip)

    end function GetSchmidt


    function FCCGeomConstructor(numSlip) result(this)
      type(FCCGeom), pointer :: this
      integer(kind = IKIND), intent(in) :: NumSlip

      allocate(this)

      call InitSlipGeom(this, numSlip)

    end function FCCGeomConstructor


    subroutine FCCGeomDestructor(this)
      type(FCCGeom) :: this

      call DestroySlipGeom(this)

    end subroutine FCCGeomDestructor


    subroutine SetSlipVecsFCC(this)
      class(FCCGeom) :: this

      this%fSlipVecS(:,  1) = (/ 1, -1,  0/)
      this%fSlipVecS(:,  2) = (/ 0, -1,  1/)
      this%fSlipVecS(:,  3) = (/-1,  0,  1/)
      this%fSlipVecS(:,  4) = (/ 1,  1,  0/)
      this%fSlipVecS(:,  5) = (/ 0,  1,  1/)
      this%fSlipVecS(:,  6) = (/-1,  0,  1/)
      this%fSlipVecS(:,  7) = (/ 1, -1,  0/)
      this%fSlipVecS(:,  8) = (/ 0,  1,  1/)
      this%fSlipVecS(:,  9) = (/ 1,  0,  1/)
      this%fSlipVecS(:, 10) = (/ 1,  1,  0/)
      this%fSlipVecS(:, 11) = (/ 0, -1,  1/)
      this%fSlipVecS(:, 12) = (/ 1,  0,  1/)

      this%fSlipVecM(:,  1) = (/ 1,  1,  1/)
      this%fSlipVecM(:,  2) = (/ 1,  1,  1/)
      this%fSlipVecM(:,  3) = (/ 1,  1,  1/)
      this%fSlipVecM(:,  4) = (/ 1, -1,  1/)
      this%fSlipVecM(:,  5) = (/ 1, -1,  1/)
      this%fSlipVecM(:,  6) = (/ 1, -1,  1/)
      this%fSlipVecM(:,  7) = (/-1, -1,  1/)
      this%fSlipVecM(:,  8) = (/-1, -1,  1/)
      this%fSlipVecM(:,  9) = (/-1, -1,  1/)
      this%fSlipVecM(:, 10) = (/-1,  1,  1/)
      this%fSlipVecM(:, 11) = (/-1,  1,  1/)
      this%fSlipVecM(:, 12) = (/-1,  1,  1/)

    end subroutine



    function BCCGeomConstructor(numSlip) result(this)
      integer(kind = IKIND), intent(in) :: NumSlip
      type(BCCGeom), pointer :: this

      allocate(this)
      call InitSlipGeom(this, numSlip)

    end function BCCGeomConstructor

    subroutine BCCGeomDestructor(this)
      type(BCCGeom) :: this

      call DestroySlipGeom(this)

    end subroutine BCCGeomDestructor

    subroutine SetSlipVecsBCC(this)
      class(BCCGeom) :: this


    end subroutine


    function UserGeomConstructor(numSlip, fID) result(this)
      type(UserGeom), pointer :: this
      integer(kind = IKIND), intent(in) :: numSlip
      integer(kind = IKIND), intent(in) :: fID

      allocate(this)

      this%fID = fID
      call InitSlipGeom(this, numSlip)

    end function UserGeomConstructor


    subroutine UserGeomDestructor(this)
      type(UserGeom) :: this

      call DestroySlipGeom(this)

    end subroutine UserGeomDestructor


    subroutine SetSlipVecsUser(this)
      use utils, only : ReadLine
      class(UserGeom) :: this

      character(len = 256) :: line
      integer :: i

      do i = 1, this%fNumSlip
        call ReadLine(this%fID, line)
        read(line, *) this%fSlipVecM(1:3, i), this%fSlipVecS(1:3, i)
      end do
        

    end subroutine
    
  end module modSlipGeom



  program main
    use modSlipGeom

    class(SlipGeom), pointer :: theSlipGeom
    type(FCCGeom),  target :: theFCCGeom
    type(BCCGeom),  target :: theBCCGeom
    type(UserGeom), target :: theUserGeom
    real(kind = 8) :: schmidt(3, 3)
    
    integer :: GeomType, numSlip, fID, i

    GeomType = 1
    fID = 100

    open(unit=fID, file = "test.in", status="old") 

    if (GeomType == 1) then
      theSlipGeom => FCCGeom(12)
    else if (GeomType == 2) then
      theSlipGeom => BCCGeom(12)
    else if (GeomType == 3) then
      theSlipGeom => UserGeom(2, fID)
    end if

    numSlip = theSlipGeom%GetNumSlip()

    do i = 1, numSlip
      schmidt = theSlipGeom%GetSchmidt(i)
      write(*, *) "Slip ", i
      write(*, *) schmidt
    end do
    
    close(unit=fID)

  end program main

    
