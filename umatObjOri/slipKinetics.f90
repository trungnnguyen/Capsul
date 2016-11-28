module typeSlipKine
  use utils, only : RKIND, IKIND, LKIND
  implicit none


  integer(kind = IKIND), public, parameter :: kSKPowLaw1    = 1
  integer(kind = IKIND), public, parameter :: kSKFrostAshby = 2

  type, abstract, public :: SlipKine
    integer(kind = IKIND) :: fNumSlipSys
    logical(kind = LKIND) :: fIsTempDep
  contains
    procedure(SlipRateIntf),          public, deferred :: SlipRate
    procedure(DSlipRateDTauReslIntf), public, deferred :: DSlipRateDTauResl
    procedure(DSlipRateDTauCritIntf), public, deferred :: DSlipRateDTauCrit
    procedure(StressDivergenceStateIntf),  public, deferred :: StressDivergenceState

    procedure, nopass, public :: Power
  end type SlipKine

  abstract interface
    function SlipRateIntf(this, tauResl, tauCrit, temp) result(gammaDot)
      import SlipKine
      integer, parameter :: RKIND = 8
      class(SlipKine),    intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND)             :: gammaDot(this%fNumSlipSys)
    end function SlipRateIntf
    
    function DSlipRateDTauReslIntf(this, tauResl, tauCrit, temp) result(dGammaDotDTauResl)
      import SlipKine
      integer, parameter :: RKIND = 8
      class(SlipKine),    intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND)             :: dGammaDotDTauResl(this%fNumSlipSys)
    end function DSlipRateDTauReslIntf

    function DSlipRateDTauCritIntf(this, tauResl, tauCrit, temp) result(dGammaDotDTauCrit)
      import SlipKine
      integer, parameter :: RKIND = 8
      class(SlipKine),    intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND)             :: dGammaDotDTauCrit(this%fNumSlipSys)
    end function DSlipRateDTauCritIntf

    function StressDivergenceStateIntf(this, tauResl, tauCrit) result(diverged)
      import SlipKine
      integer, parameter :: RKIND = 8
      integer, parameter :: LKIND = 1
      class(SlipKine),    intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      logical(kind = LKIND)             :: diverged
    end function StressDivergenceStateIntf


  end interface
contains

  function Power(x, y) result(rlt)
    real(kind = RKIND), intent(in) :: x
    real(kind = RKIND), intent(in) :: y
    real(kind = RKIND) :: rlt
    real(kind = RKIND), parameter :: MINDBL = 1.0d-15

    if (dabs(x) <= MINDBL) then
      if (dabs(y) <= MINDBL) then
        rlt = 1.0d0
      else if (y > MINDBL) then
        rlt = 0.0d0
      else if (y < -MINDBL) then
        rlt = 1.0d300
      end if
    else
      rlt = y*log10(dabs(x))
      if (rlt > 300.0d0) then
        rlt = 1.0d300
      else
        rlt = 10.0d0**rlt
      end if

      if (x < 0.0d0) rlt = -rlt
    end if

  end function Power

        
end module typeSlipKine



module typeSlipKinePowLaw1
  use typeSlipKine
  use utils, only : RKIND, IKIND, LKIND
  implicit none

  type, extends(SlipKine), public :: SlipKinePowLaw1
    real(kind = RKIND) :: fGammaDot0
    real(kind = RKIND) :: fXm
    real(kind = RKIND) :: fQActive
  contains
    procedure, public :: SlipRate
    procedure, public :: DSlipRateDTauResl
    procedure, public :: DSlipRateDTauCrit
    procedure, public :: StressDivergenceState
  end type

  interface SlipKinePowLaw1
    module procedure  Constructor
  end interface


contains

  function Constructor(fID, isTempDep, numSlipSys) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    logical(kind = LKIND), intent(in) :: isTempDep
    integer(kind = IKIND), intent(in) :: numSlipSys
    type(SlipKinePowLaw1), pointer    :: this

    character(len = 256) :: line

    allocate(this)
    this%fNumSlipSys = numSlipSys
    this%fIsTempDep  = isTempDep

    call ReadLine(fID, line)
    read(line, *) this%fGammaDot0, this%fXm, this%fQActive
  
  end function Constructor



  function SlipRate(this, tauResl, tauCrit, temp) result(gammaDot)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),     intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),     intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),     intent(in) :: temp
    real(kind = RKIND)                 :: gammaDot(this%fNumSlipSys)

    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: tauSign(this%fNumSlipSys)
    real(kind = RKIND) :: absTauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: factor

    real(kind = RKIND) :: aux1, aux2
    integer(kind = IKIND) :: i


    if (this%fIsTempDep) then
      factor = exp(-this%fQActive/temp)
    else
      factor = 1.0d0
    end if

    tauRatio = tauResl/tauCrit
    absTauRatio = abs(tauRatio)


    aux1 = factor*this%fGammaDot0
    aux2 = 1.0d0/this%fXm - 1.0d0
    tauSign = sign(1.0d0, tauRatio)
    do i = 1, this%fNumSlipSys
!      SlipRate(i) = aux1*tauRatio(i)*this%Power(dabs(tauRatio(i)), aux2)
      gammaDot(i) = aux1*tauSign(i)*abs(tauRatio(i))**(1.0d0/this%fXm)
    end do

  end function SlipRate



  function DSlipRateDTauResl(this, tauResl, tauCrit, temp) result(dGammaDotDTauResl)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),     intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),     intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),     intent(in) :: temp
    real(kind = RKIND)                 :: dGammaDotDTauResl(this%fNumSlipSys)

    real(kind = RKIND)    :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND)    :: factor
    real(kind = RKIND)    :: aux1, aux2
    integer(kind = IKIND) :: i


    if (this%fIsTempDep) then
      factor = exp(-this%fQActive/temp)
    else
      factor = 1.0d0
    end if

    tauRatio = tauResl/tauCrit

    aux1 = factor*this%fGammaDot0/this%fXm
    aux2 = 1.0d0/this%fXm - 1.0d0
    do i = 1, this%fNumSlipSys
     ! DSlipRateDTauResl(i) = aux1*this%Power(dabs(tauRatio(i)), aux2)/tauCrit(i)
      dGammaDotDTauResl(i) = aux1*dabs(tauRatio(i))**aux2/tauCrit(i)
    end do

  end function DSlipRateDTauResl



  function DSlipRateDTauCrit(this, tauResl, tauCrit, temp)  result(dGammaDotDTauCrit)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),     intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),     intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),     intent(in) :: temp
    real(kind = RKIND)                 :: dGammaDotDTauCrit(this%fNumSlipSys)

    real(kind = RKIND)    :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND)    :: factor

    real(kind = RKIND)    :: aux1, aux2
    integer(kind = IKIND) :: i

    if (this%fIsTempDep) then
      factor = exp(-this%fQActive/temp)
    else
      factor = 1.0d0
    end if

    tauRatio = tauResl/tauCrit
    
    aux1 = -factor*this%fGammaDot0/this%fXm
    aux2 = 1.0d0/this%fXm - 1.0d0

    do i = 1, this%fNumSlipSys
!      DSlipRateDTauCrit(i) = aux1*this%Power(dabs(tauRatio(i)), aux2)*tauResl(i)
      dGammaDotDTauCrit(i) = aux1*dabs(tauRatio(i))**aux2/tauCrit(i)*tauRatio(i)
    end do

  end function DSlipRateDTauCrit



  function StressDivergenceState(this, tauResl, tauCrit) result(diverged)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),     intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),     intent(in) :: tauCrit(this%fNumSlipSys)
    logical(kind = LKIND)              :: diverged

    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: mmBig


    tauRatio = tauResl/tauCrit
    mmBig    = 10.0d0**(150.0d0*this%fXm)
    diverged = any(dabs(tauRatio) > mmBig)

  end function StressDivergenceState

    
end module typeSlipKinePowLaw1



module typeSlipKineFrostAshby
  use typeSlipKine
  use utils, only : RKIND, IKIND, LKIND
  implicit none

  integer(kind = RKIND), parameter :: kConstTauCritT   = 1
  integer(kind = RKIND), parameter :: kVaryingTauCritT = 2

  type, extends(SlipKine), public :: SlipKineFrostAshby
    real(kind = RKIND)    :: fGammaDot0
    real(kind = RKIND)    :: fEAct         ! Activation Energy
    real(kind = RKIND)    :: fp         
    real(kind = RKIND)    :: fq
    real(kind = RKIND)    :: fKb
    integer(kind = IKIND) :: fTauCritTCode  
    real(kind = RKIND)    :: fYeta
    real(kind = RKIND), allocatable :: fTauCritT(:)    
  contains
    procedure, public :: SlipRate
    procedure, public :: DSlipRateDTauResl
    procedure, public :: DSlipRateDTauCrit
    procedure, public :: StressDivergenceState
    final             :: Destructor
  end type

  interface SlipKineFrostAshby
    module procedure  Constructor
  end interface

contains

  function Constructor(fID, isTempDep, numSlipSys) result(this)
    use utils, only : ReadLine

    integer(kind = IKIND), intent(in) :: fID
    logical(kind = LKIND), intent(in) :: isTempDep
    integer(kind = IKIND), intent(in) :: numSlipSys
    type(SlipKineFrostAshby), pointer    :: this

    real(kind = RKIND)   :: tmpVar
    character(len = 256) :: line


    allocate(this)
    this%fNumSlipSys = numSlipSys
    this%fIsTempDep  = isTempDep

    allocate(this%fTauCritT(this%fNumSlipSys)) 

    call ReadLine(fID, line)
    read(line, *) this%fGammaDot0, this%fEAct, this%fp, this%fq,  this%fKb, this%fTauCritTCode, tmpVar
    if (this%fTauCritTCode == kConstTauCritT) then
       this%fTauCritT = tmpVar
       this%fYeta     = -1.0d0   ! fYeta is not used in this case
    else if (this%fTauCritTCode == kVaryingTauCritT) then
    ! this%fYeta is actually yeta/(1.0+yeta)
       this%fYeta     = tmpVar/(1.0d0 + tmpVar)
       this%fTauCritT = 0.0d0
    else
      write(*, *) "Error! Unknown code for the constantness of thermal part of TauCrit"
      stop
    end if 
  
  end function Constructor



  subroutine Destructor(this)
    type(SlipKineFrostAshby) :: this

    if (allocated(this%fTauCritT)) deallocate(this%fTauCritT)

  end subroutine Destructor



  function SlipRate(this, tauResl, tauCrit, temp)  result(gammaDot)
    class(SlipKineFrostAshby), intent(in) :: this
    real(kind = RKIND),        intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),        intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),        intent(in) :: temp
    real(kind = RKIND)                    :: gammaDot(this%fNumSlipSys)

    real(kind = RKIND) :: tauCritAT(this%fNumSlipSys)
    real(kind = RKIND) :: tauCritT(this%fNumSlipSys)
    real(kind = RKIND) :: tauReslT(this%fNumSlipSys)
    real(kind = RKIND) :: tauSign(this%fNumSlipSys)
    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: dltG
    real(kind = RKIND) :: aux1
    integer(kind = IKIND) :: i



    if (this%fTauCritTCode == kConstTauCritT) then
      tauCritT = this%fTauCritT
    else if (this%fTauCritTCode == kVaryingTauCritT) then
      !tauCritT = tauCrit*(this%fYeta/(1.0d0 + this%fYeta))
      tauCritT = tauCrit*this%fYeta
    end if
    tauCritAT = tauCrit - tauCritT 
    tauReslT  = dabs(tauResl) - tauCritAT 
    tauRatio  = tauReslT/tauCritT

    tauSign = sign(1.0d0, tauResl)

    aux1 = this%fKb*temp
    do i = 1, this%fNumSlipSys
      if (tauReslT(i) <= 0.0d0) then
        gammaDot(i) = 0.0d0
      else if (tauReslT(i) > 0.0d0 .and. tauReslT(i) < tauCritT(i)) then
        dltG = this%fEAct*(1.0d0 - tauRatio(i)**this%fp)**this%fq
        gammaDot(i) = this%fGammaDot0*exp(-dltG/aux1)*tauSign(i)
      else
        write(*, *) "Error in SlipKineFrostAshby::SlipRate"
        stop
      end if
    end do


  end function SlipRate



  function DSlipRateDTauResl(this, tauResl, tauCrit, temp) result(dGammaDotDTauResl)
    class(SlipKineFrostAshby), intent(in) :: this
    real(kind = RKIND),        intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),        intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),        intent(in) :: temp
    real(kind = RKIND)                    :: dGammaDotDTauResl(this%fNumSlipSys)

    real(kind = RKIND) :: tauCritAT(this%fNumSlipSys)
    real(kind = RKIND) :: tauCritT(this%fNumSlipSys)
    real(kind = RKIND) :: tauReslT(this%fNumSlipSys)
    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)

    real(kind = RKIND) :: aux1, aux2, aux3, aux4, aux5, y1, dfdy1, dy1dy2, dy2dx, dfdx
    integer(kind = IKIND) :: i

    if (this%fTauCritTCode == kConstTauCritT) then
      tauCritT = this%fTauCritT
    else if (this%fTauCritTCode == kVaryingTauCritT) then
    ! this%fYeta is actually yeta/(1.0+yeta)
      tauCritT = tauCrit*this%fYeta
    end if
    tauCritAT = tauCrit - tauCritT 
    tauReslT  = dabs(tauResl) - tauCritAT 
    tauRatio  = tauReslT/tauCritT

    aux1 = this%fEAct/(this%fKb*temp)
    aux2 = aux1*this%fp*this%fq
    do i = 1, this%fNumSlipSys
      if (tauReslT(i) <= 0.0d0) then
        dGammaDotDTauResl(i) = 0.0d0
      else if (tauReslT(i) > 0 .and. tauReslT(i) < tauCritT(i)) then
        aux3   = tauRatio(i)**(this%fp - 1.0d0)
        aux4   = 1.0d0 - aux3*tauRatio(i)
        aux5   = aux4**(this%fq - 1.0d0)
        y1     = -aux1*aux5*aux4
        dfDy1  = this%fGammaDot0*exp(y1)
        dy1Dy2 = aux2*aux5*aux3
        dy2dx  = 1.0d0/tauCritT(i) 
        dfdx   = dfdy1*dy1dy2*dy2dx
        dGammaDotDTauResl(i) = dfdx
      else
        write(*, *) "Error in SlipKineFrostAshby::DSlipRateDTauResl"
        stop
      end if
    end do
     
  end function DSlipRateDTauResl




  function DSlipRateDTauCrit(this, tauResl, tauCrit, temp) result(dGammaDotDTauCrit)
    class(SlipKineFrostAshby), intent(in) :: this
    real(kind = RKIND),        intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),        intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),        intent(in) :: temp
    real(kind = RKIND)                    :: dGammaDotDTauCrit(this%fNumSlipSys)

    real(kind = RKIND) :: tauCritAT(this%fNumSlipSys)
    real(kind = RKIND) :: tauCritT(this%fNumSlipSys)
    real(kind = RKIND) :: tauReslT(this%fNumSlipSys)
    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)

    real(kind = RKIND) :: aux1, aux2, aux3, aux4, aux5, y1, dfdy1, dy1dy2, dy2dx, dfdx
    integer(kind = IKIND) :: i

    if (this%fTauCritTCode == kConstTauCritT) then
      tauCritT = this%fTauCritT
      dGammaDotDTauCrit = 0
      return
    else if (this%fTauCritTCode == kVaryingTauCritT) then
    ! this%fYeta is actually yeta/(1.0+yeta)
      tauCritT = tauCrit*this%fYeta
    end if
    tauCritAT = tauCrit - tauCritT 
    tauReslT  = dabs(tauResl) - tauCritAT 
    tauRatio  = tauReslT/tauCritT

    aux1 = this%fEAct/(this%fKb*temp)
    aux2 = aux1*this%fp*this%fq
    do i = 1, this%fNumSlipSys
      if (tauReslT(i) <= 0.0d0) then
        dGammaDotDTauCrit(i) = 0.0d0
      else if (tauReslT(i) > 0 .and. tauReslT(i) < tauCritT(i)) then
        aux3   = tauRatio(i)**(this%fp - 1.0d0)
        aux4   = 1.0d0 - aux3*tauRatio(i)
        aux5   = aux4**(this%fq - 1.0d0)
        y1     = -aux1*aux5*aux4
        dfDy1  = this%fGammaDot0*exp(y1)*sign(1.0d0, tauResl(i))
        dy1Dy2 = aux2*aux5*aux3
        if (this%fTauCritTCode == kConstTauCritT) then
          dy2dx  = -1.0d0/tauCritT(i)
        else if (this%fTauCritTCode == kVaryingTauCritT) then
          dy2dx  = -tauReslT(i)/(tauCritT(i)*tauCritT(i))*this%fYeta
        end if
        dGammaDotDTauCrit(i) = dfdy1*dy1dy2*dy2dx
      else
        write(*, *) "Error in SlipKineFrostAshby::DSlipRateDTauCrit"
        stop
      end if
    end do
    
  end function DSlipRateDTauCrit


 

  function StressDivergenceState(this, tauResl, tauCrit) result(diverged)
    class(SlipKineFrostAshby), intent(in) :: this
    real(kind = RKIND),        intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),        intent(in) :: tauCrit(this%fNumSlipSys)
    logical(kind = LKIND)              :: diverged

    real(kind = RKIND) :: tauCritAT(this%fNumSlipSys)
    real(kind = RKIND) :: tauCritT(this%fNumSlipSys)
    real(kind = RKIND) :: tauReslT(this%fNumSlipSys)


    if (this%fTauCritTCode == kConstTauCritT) then
      tauCritT = this%fTauCritT
    else if (this%fTauCritTCode == kVaryingTauCritT) then
    ! this%fYeta is actually yeta/(1.0+yeta)
      tauCritT = tauCrit*this%fYeta
    end if
    tauCritAT = tauCrit - tauCritT 
    tauReslT  = dabs(tauResl) - tauCritAT 

    diverged = any((tauReslT > tauCritT))

  end function StressDivergenceState

    
end module typeSlipKineFrostAshby


