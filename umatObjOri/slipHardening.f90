  module typeSlipHard
    use utils, only : RKIND, IKIND, LKIND
    implicit none
  
  
    integer(kind = RKIND), public, parameter :: kAssaNeed = 1
    integer(kind = RKIND), public, parameter :: kVoceTome = 2
    integer(kind = RKIND), public, parameter :: kVoceKock = 3
  
    type, abstract, public :: SlipHard
      integer(kind = IKIND)           :: fNumSlipSys
      logical(kind = LKIND)           :: fIsTempDep
      real(kind = RKIND), allocatable :: fMatProp(:, :)
      real(kind = RKIND), allocatable :: fQab(:, :) 
    contains

      ! subs and funs used by external call. Main purpose of this class.
      procedure(InitialCRSSIntf),     public, deferred :: InitialCRSS
      procedure,  public :: TauCritRate
      procedure,  public :: DTauCritDGamma

      ! subs and funs used by this class and its derived classes.
      procedure,  public :: InitHardQnts
      procedure,  public :: Destroy
      procedure(HardMatrixIntf),   public, deferred :: HardMatrix
      procedure(DHardMatrixIntf),  public, deferred :: DHardMatrix
  
    end type SlipHard
  
    abstract interface
          
      function InitialCRSSIntf(this) result(tauCrit0)
        import SlipHard
        integer, parameter :: RKIND = 8
        class(SlipHard),    intent(in) :: this
        real(kind = RKIND) :: tauCrit0(this%fNumSlipSys)
      end function InitialCRSSIntf

      function HardMatrixIntf(this, gamma, gammaDot, tauCrit, temp) result(HardMatx)
        import SlipHard
        integer, parameter :: RKIND = 8
        class(SlipHard),    intent(in) :: this
        real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
        real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
        real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
        real(kind = RKIND), intent(in) :: temp
        real(kind = RKIND)             :: HardMatx(this%fNumSlipSys, this%fNumSlipSys)
      end function HardMatrixIntf

      function DHardMatrixIntf(this, gamma, gammaDot, tauCrit, temp) result(DHardMatx)
        import SlipHard
        integer, parameter :: RKIND = 8
        class(SlipHard),    intent(in) :: this
        real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
        real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
        real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
        real(kind = RKIND), intent(in) :: temp
        real(kind = RKIND)             :: DHardMatx(this%fNumSlipSys, this%fNumSlipSys)
      end function DHardMatrixIntf

    end interface
  
  contains
  
    subroutine InitHardQnts(this, fID, numSlipSet, SlipSet, kNumMatProp, kNumMatAux)
      use utils, only : ReadLine
  
      class(SlipHard) :: this
      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: numSlipSet
      integer(kind = IKIND), intent(in) :: SlipSet(this%fNumSlipSys)
      integer(kind = IKIND), intent(in) :: kNumMatProp
      integer(kind = IKIND), intent(in) :: kNumMatAux
  
      real(kind = RKIND) :: auxMatProp(numSlipSet, kNumMatProp+kNumMatAux)
      character(len = 256)  :: line
      integer(kind = IKIND) :: i, j
  
  
      allocate(this%fMatProp(this%fNumSlipSys, kNumMatProp+kNumMatAux))
  
      ! Read material properties for hardening in each slip set
      do i = 1, numSlipSet
        call ReadLine(fID, line)
        read(line, *) (auxMatProp(i, j), j = 1, kNumMatProp)
      end do
  
      do i = 1, this%fNumSlipSys
        this%fMatProp(i, 1:kNumMatProp) = auxMatProp(slipSet(i), 1:kNumMatProp)
      end do
  
  
    end subroutine InitHardQnts
  
  
    subroutine Destroy(this)
      class(SlipHard) :: this
  
      if (allocated(this%fMatProp))    deallocate(this%fMatProp)
      if (allocated(this%fQab))        deallocate(this%fQab)
  
    end subroutine Destroy
  
           
  
    function TauCritRate(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHard),    intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: TauCritRate(this%fNumSlipSys)
  
      real(kind = RKIND) :: hh(this%fNumSlipSys, this%fNumSlipSys)
      real(kind = RKIND) :: absGammaDot(this%fNumSlipSys)
      integer(kind = IKIND) :: i

      hh = this%HardMatrix(gamma, gammaDot, tauCrit, temp)
      absGammaDot = abs(gammaDot)
      do i = 1, this%fNumSlipSys
        TauCritRate(i) = dot_product(hh(i, :), absGammaDot)
      end do
  
    end function TauCritRate
  

    function DTauCritDGamma(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHard),    intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: DTauCritDGamma(this%fNumSlipSys, this%fNumSlipSys)
  
      real(kind = RKIND) :: gammaDotSign(this%fNumSlipSys)
      real(kind = RKIND) :: absGammaDot(this%fNumSlipSys)
      real(kind = RKIND) :: hh(this%fNumSlipSys, this%fNumSlipSys)
      real(kind = RKIND) :: gg(this%fNumSlipSys, this%fNumSlipSys)
      integer(kind = RKIND) :: ii, jj

      gammaDotSign = sign(1.0d0, gammaDot)
      absGammaDot  = abs(gammaDot)
      hh = this%HardMatrix( gamma, gammaDot, tauCrit, temp)
      gg = this%DHardMatrix(gamma, gammaDot, tauCrit, temp)

      do ii = 1, this%fNumSlipSys
        do jj = 1, this%fNumSlipSys
          DTauCritDGamma(ii, jj) = hh(ii, jj)*gammaDotSign(jj) + gg(ii, jj)*absGammaDot(jj)
        end do
      end do

    end function DTauCritDGamma
  

  end module typeSlipHard
  
  
  module typeSlipHardAN
    use typeSlipHard
    use utils, only : RKIND, IKIND
    implicit none
  
  
    type, extends(SlipHard), public :: SlipHardAN
    contains
      procedure, public :: InitialCRSS
      procedure, public :: HardMatrix
      procedure, public :: DHardMatrix
      final             :: Destructor
    end type
  
    interface SlipHardAN
      module procedure Constructor
    end interface
  
  contains
  
    function Constructor(fID, isTempDep, numSlipSys, numSlipSet, slipSet) result(this)
      use utils, only : ReadLine
  
      integer(kind = IKIND), intent(in) :: fID
      logical(kind = LKIND), intent(in) :: isTempDep
      integer(kind = IKIND), intent(in) :: numSlipSys
      integer(kind = IKIND), intent(in) :: numSlipSet
      integer(kind = IKIND), intent(in) :: slipSet(numSlipSys)
      type(SlipHardAN), pointer         :: this
  
      integer(kind = IKIND), parameter :: kNumMatProp = 3
      integer(kind = IKIND), parameter :: kNumMatAux  = 1
      real(kind = RKIND) :: auxQab(numSlipSet, numSlipSet)

      character(len = 256) :: line
      integer :: i, j
  

      allocate(this)

      this%fIsTempDep  = isTempDep
  
      this%fNumSlipSys = numSlipSys
      allocate(this%fQab(this%fNumSlipSys, this%fNumSlipSys)) 
      do i = 1, numSlipSet
        call ReadLine(fID, line) 
        read(line, *) (auxQab(i, j), j = 1, numSlipSet)
      end do

      do i = 1, this%fNumSlipSys
        do j = 1, this%fNumSlipSys
          this%fQab(i, j) = auxQab(slipSet(i), slipSet(j))
        end do
      end do
  
      call this%InitHardQnts(fID, numSlipSet, slipSet, kNumMatProp, kNumMatAux)

      ! h0/(taus - tau0) 
      this%fMatProp(:, 4) = this%fMatProp(:, 3)/(this%fMatProp(:, 2) - this%fMatProp(:, 1))
  
    end function Constructor
  
  
    subroutine Destructor(this)
      type(SlipHardAN) :: this
  
      call this%Destroy()
  
    end subroutine Destructor


    function InitialCRSS(this) result(tauCrit0)
      class(SlipHardAN),  intent(in) :: this
      real(kind = RKIND) :: tauCrit0(this%fNumSlipSys)

      tauCrit0 = this%fMatProp(:, 1)

    end function InitialCRSS


    function HardMatrix(this, gamma, gammaDot, tauCrit, temp)
      use algebra, only : sech

      class(SlipHardAN),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: HardMatrix(this%fNumSlipSys, this%fNumSlipSys)

      real(kind = RKIND) :: h0(this%fNumSlipSys)
      real(kind = RKIND) :: aux(this%fNumSlipSys)
      real(kind = RKIND) :: hii
      integer(kind = IKIND) :: ii, jj

      h0  = this%fMatProp(:, 3)
      aux = sech(this%fMatProp(:, 4) * sum(gamma))
      do ii = 1, this%fNumSlipSys
        hii = h0(ii)*aux(ii)*aux(ii) 
        do jj = 1, this%fNumSlipSys
          HardMatrix(ii, jj) = this%fQab(ii, jj)*hii
        end do
      end do

    end function HardMatrix


    function DHardMatrix(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHardAN),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: DHardMatrix(this%fNumSlipSys, this%fNumSlipSys)


      DHardMatrix = 0.0d0

    end function DHardMatrix

  end module typeSlipHardAN
  
  
  
  module typeSlipHardVT
    use typeSlipHard
    use utils, only : RKIND, IKIND
    implicit none
  
  
    type, extends(SlipHard), public :: SlipHardVT
    contains
      procedure, public :: InitialCRSS
      procedure, public :: HardMatrix
      procedure, public :: DHardMatrix
      final             :: Destructor
    end type
  
    interface SlipHardVT
      module procedure Constructor
    end interface
  
  contains
  
    function Constructor(fID, isTempDep, numSlipSys, numSlipSet, slipSet) result(this)
      integer(kind = IKIND), intent(in) :: fID
      logical(kind = LKIND), intent(in) :: isTempDep
      integer(kind = IKIND), intent(in) :: numSlipSys
      integer(kind = IKIND), intent(in) :: numSlipSet
      integer(kind = IKIND), intent(in) :: slipSet(numSlipSys)
      type(SlipHardVT), pointer         :: this
  
      integer(kind = IKIND), parameter :: kNumMatProp = 4
      integer(kind = IKIND), parameter :: kNumMatAux  = 1
      real(kind = RKIND) :: auxQab(numSlipSet, numSlipSet)

      character(len = 256) :: line
      integer :: i, j
  
      allocate(this)
      this%fIsTempDep  = isTempDep

      this%fNumSlipSys = numSlipSys
      allocate(this%fQab(this%fNumSlipSys, this%fNumSlipSys)) 
      do i = 1, numSlipSet
        call ReadLine(fID, line) 
        read(line, *) (auxQab(i, j), j = 1, numSlipSet)
      end do
  
      do i = 1, this%fNumSlipSys
        do j = 1, this%fNumSlipSys
          this%fQab(i, j) = auxQab(slipSet(i), slipSet(j))
        end do
      end do
  
  
      call this%InitHardQnts(fID, numSlipSet, slipSet, kNumMatProp, kNumMatAux)

      this%fMatProp(:, 2) = this%fMatProp(:, 2) - this%fMatProp(:, 1)
      this%fMatProp(:, 5) = this%fMatProp(:, 3)/this%fMatProp(:, 2)
  
    end function Constructor
  
  
    subroutine Destructor(this)
      type(SlipHardVT) :: this
  
      call this%Destroy()
  
    end subroutine Destructor


    function InitialCRSS(this) result(tauCrit0)
      class(SlipHardVT),  intent(in) :: this
      real(kind = RKIND) :: tauCrit0(this%fNumSlipSys)

      tauCrit0 = this%fMatProp(:, 1)

    end function InitialCRSS


    function HardMatrix(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHardVT),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: HardMatrix(this%fNumSlipSys, this%fNumSlipSys)

      real(kind = RKIND) :: taus(this%fNumSlipSys)
      real(kind = RKIND) :: h1(this%fNumSlipSys)
      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: aux(this%fNumSlipSys)
      real(kind = RKIND) :: hii

      integer(kind = IKIND) :: ii, jj


      taus = this%fMatProp(:, 2)
      h1   = this%fMatProp(:, 4)

      gammaTot = sum(gamma)
      aux = dexp(this%fMatProp(:, 5)*(-gammaTot))
      do ii = 1, this%fNumSlipSys
        hii = h1(ii)*(1.0d0 - aux(ii)) + (taus(ii) + h1(ii)*gammaTot) * this%fMatProp(ii, 5)*aux(ii)
        do jj = 1, this%fNumSlipSys
          HardMatrix(ii, jj) = this%fQab(ii, jj)*hii
        end do
      end do

    end function HardMatrix


    function DHardMatrix(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHardVT),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: DHardMatrix(this%fNumSlipSys, this%fNumSlipSys)


      DHardMatrix = 0.0d0

    end function DHardMatrix


  end module typeSlipHardVT
  
  
  module typeSlipHardVK
    use typeSlipHard
    use utils, only : RKIND, IKIND
    implicit none
  
  
    type, extends(SlipHard), public :: SlipHardVK
    contains
      procedure, public :: InitialCRSS
      procedure, public :: HardMatrix
      procedure, public :: DHardMatrix
      final             :: Destructor
    end type
  
    interface SlipHardVK
      module procedure Constructor
    end interface
  
  contains
  
    function Constructor(fID, isTempDep, numSlipSys, numSlipSet, slipSet) result(this)
      integer(kind = IKIND), intent(in) :: fID
      logical(kind = LKIND), intent(in) :: isTempDep
      integer(kind = IKIND), intent(in) :: numSlipSys
      integer(kind = IKIND), intent(in) :: numSlipSet
      integer(kind = IKIND), intent(in) :: slipSet(numSlipSys)
      type(SlipHardVK), pointer         :: this
  
      integer(kind = IKIND), parameter  :: dummyNumSlipSet = 1
      integer(kind = IKIND) :: dummySlipSet(dummyNumSlipSet)

      integer(kind = IKIND), parameter :: kNumMatProp = 5
      integer(kind = IKIND), parameter :: kNumMatAux  = 1
      integer :: i, j
  
      allocate(this)

      this%fIsTempDep  = isTempDep
  
      this%fNumSlipSys = numSlipSys
      allocate(this%fQab(this%fNumSlipSys, this%fNumSlipSys)) 
      this%fQab = 1.0d0
  
      dummySlipSet(1:dummyNumSlipSet) = 1
      call this%InitHardQnts(fID, dummyNumSlipSet, dummySlipSet, kNumMatProp, kNumMatAux)
  
    end function Constructor
  
  
    subroutine Destructor(this)
      type(SlipHardVK) :: this
  
      call this%Destroy()
  
    end subroutine Destructor


    function InitialCRSS(this) result(tauCrit0)
      class(SlipHardVK),  intent(in) :: this
      real(kind = RKIND) :: tauCrit0(this%fNumSlipSys)

      tauCrit0 = this%fMatProp(1, 1)

    end function InitialCRSS


    function HardMatrix(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHardVK),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: HardMatrix(this%fNumSlipSys, this%fNumSlipSys)

      real(kind = RKIND) :: tau0  
      real(kind = RKIND) :: taus0  
      real(kind = RKIND) :: h0  
      real(kind = RKIND) :: gammas0  
      real(kind = RKIND) :: ma

      real(kind = RKIND) :: gammaDotSum
      real(kind = RKIND) :: taus
      real(kind = RKIND) :: aux
      real(kind = RKIND) :: mm
      integer(kind = IKIND) :: ii


      tau0    = this%fMatProp(1,  1)
      taus0   = this%fMatProp(1,  2)
      h0      = this%fMatProp(1,  3)
      gammas0 = this%fMatProp(1,  4)
      ma      = this%fMatProp(1,  5)

      mm = ma

      gammaDotSum = sum(abs(gammaDot))
      taus = taus0 * (gammaDotSum/gammas0) ** (mm)
      aux = h0/(taus - tau0)
      do ii = 1, this%fNumSlipSys
        HardMatrix(1:this%fNumSlipSys, ii) = aux*(taus - tauCrit(ii))
      end do

    end function HardMatrix


    function DHardMatrix(this, gamma, gammaDot, tauCrit, temp)
      class(SlipHardVK),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: DHardMatrix(this%fNumSlipSys, this%fNumSlipSys)

      real(kind = RKIND) :: tau0  
      real(kind = RKIND) :: taus0  
      real(kind = RKIND) :: h0  
      real(kind = RKIND) :: gammas0  
      real(kind = RKIND) :: ma


      real(kind = RKIND) :: gammaDotSign(this%fNumSlipSys)
      real(kind = RKIND) :: gammaDotSum
      real(kind = RKIND) :: mm, taus
      real(kind = RKIND) :: aux1, aux2, aux3 
      integer(kind = IKIND) :: i, j


      tau0    = this%fMatProp(1, 1)
      taus0   = this%fMatProp(1, 2)
      h0      = this%fMatProp(1, 3)
      gammas0 = this%fMatProp(1, 4)
      ma      = this%fMatProp(1, 5)

!      mm = temp/ma
      mm = ma
      
      gammaDotSign = sign(1.0d0, gammaDot)
      gammaDotSum = sum(abs(gammaDot))
      if (gammaDotSum < 1.0d-12) then
        DHardMatrix = 0.0d0
        return 
      end if

      taus = taus0 * (gammaDotSum/gammas0) ** mm
      aux1 = taus - tau0
      !aux2 = -h0*mm*taus/(aux1 * aux1 * gammaDotSum)
      aux2 = h0*mm*taus0*(gammaDotSum/gammas0)**(mm-1)/(aux1 * aux1 * gammas0)
      do i = 1, this%fNumSlipSys
        aux3 = aux2*(tauCrit(i) - tau0)
        do j = 1, this%fNumSlipSys
          DHardMatrix(i, j) = aux3*gammaDotSign(j)
        end do
      end do

    end function DHardMatrix

  end module typeSlipHardVK
  
  
