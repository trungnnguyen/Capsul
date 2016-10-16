module typeSlipKine
  use utils, only : RKIND, IKIND, LKIND
  implicit none


  integer(kind = IKIND), public, parameter :: kPowLaw1 = 1

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
    
    function DSlipRateDTauReslIntf(this, tauResl, tauCrit, temp) result(DGammaDotDTauResl)
      import SlipKine
      integer, parameter :: RKIND = 8
      class(SlipKine),    intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND)             :: DGammaDotDTauResl(this%fNumSlipSys)
    end function DSlipRateDTauReslIntf

    function DSlipRateDTauCritIntf(this, tauResl, tauCrit, temp) result(DGammaDotDTauCrit)
      import SlipKine
      integer, parameter :: RKIND = 8
      class(SlipKine),    intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND)             :: DGammaDotDTauCrit(this%fNumSlipSys)
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



  function SlipRate(this, tauResl, tauCrit, temp)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),    intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),    intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),    intent(in) :: temp
    real(kind = RKIND)                :: SlipRate(this%fNumSlipSys)

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
      SlipRate(i) = aux1*tauSign(i)*abs(tauRatio(i))**(1.0d0/this%fXm)
    end do

  end function SlipRate



  function DSlipRateDTauResl(this, tauResl, tauCrit, temp)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),    intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),    intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),    intent(in) :: temp
    real(kind = RKIND) :: DSlipRateDTauResl(this%fNumSlipSys)

    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: factor

    real(kind = RKIND) :: aux1, aux2
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
      DSlipRateDTauResl(i) = aux1*dabs(tauRatio(i))**aux2/tauCrit(i)
    end do

  end function DSlipRateDTauResl



  function DSlipRateDTauCrit(this, tauResl, tauCrit, temp)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),    intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),    intent(in) :: tauCrit(this%fNumSlipSys)
    real(kind = RKIND),    intent(in) :: temp
    real(kind = RKIND) :: DSlipRateDTauCrit(this%fNumSlipSys)

    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: factor

    real(kind = RKIND) :: aux1, aux2
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
      DSlipRateDTauCrit(i) = aux1*dabs(tauRatio(i)**aux2)*tauResl(i)
    end do


  end function DSlipRateDTauCrit


  function StressDivergenceState(this, tauResl, tauCrit) result(diverged)
    class(SlipKinePowLaw1), intent(in) :: this
    real(kind = RKIND),    intent(in) :: tauResl(this%fNumSlipSys) 
    real(kind = RKIND),    intent(in) :: tauCrit(this%fNumSlipSys)
    logical(kind = LKIND) :: diverged

    real(kind = RKIND) :: tauRatio(this%fNumSlipSys)
    real(kind = RKIND) :: mmBig


    tauRatio = tauResl/tauCrit
    mmBig = 10.0d0**(150.0d0*this%fXm)
    diverged = any(dabs(tauRatio) > mmBig)

  end function StressDivergenceState

    
end module typeSlipKinePowLaw1


