  module typeCrysPlasMat
    use utils, only : RKIND, IKIND, LKIND
    use typeCrysElas,    only : CrysElas
    use typeThermoExpa,  only : ThermoExpa
    use typeSlipGeom,    only : SlipGeom
    use typeSlipKine,    only : SlipKine
    use typeSlipHard,    only : SlipHard
    implicit none
  

    type,  public :: CrysPlasMat

      class(CrysElas),   pointer :: fCrysElas    => null()
      class(ThermoExpa), pointer :: fThermoExpa  => null()
      class(SlipGeom),   pointer :: fSlipGeom    => null()
      class(SlipKine),   pointer :: fSlipKine    => null()
      class(SlipHard),   pointer :: fSlipHard    => null()

      integer(kind = IKIND) :: fNumSlipSys
      real(kind = RKIND)    :: fTempInit
  
      real(kind = RKIND)    :: fCijkl(3, 3, 3, 3)
      real(kind = RKIND), allocatable :: fSchmidt(:, :, :)
  
      logical(kind = LKIND) :: fIsTempDep
  
    contains
      ! universal service interface
      procedure, public  :: GetElasticModuli

      procedure, public  :: GetNumSlipSys 
      procedure, public  :: GetSchmidt

      procedure, public  :: GetSlipRate
      procedure, public  :: GetDSlipRateDTauResl
      procedure, public  :: GetDSlipRateDTauCrit
      procedure, public  :: GetStressDivergenceState

      procedure, public  :: GetInitialCRSS
      procedure, public  :: GetTauCritRate
      procedure, public  :: GetDTauCritDGamma


      procedure, private :: ParseInput
      final              :: Destructor

    end type CrysPlasMat

    interface CrysPlasMat
      module procedure Constructor
    end interface
  
  contains
  
    function Constructor(fID, isTempDep, tempInit) result(this)
      integer(kind = IKIND), intent(in) :: fID
      logical(kind = LKIND), intent(in) :: isTempDep
      real(kind = RKIND),    intent(in) :: tempInit
      type(CrysPlasMat),     pointer    :: this
  
      allocate(this)
      this%fIsTempDep = isTempDep
      this%fTempInit = tempInit

      call this%ParseInput(fID)
      
      allocate(this%fSchmidt(3, 3, this%fNumSlipSys))
      this%fSchmidt = this%fSlipGeom%Schmidt()

      if (.not. this%fIsTempDep) then
        this%fCijkl = this%fCrysElas%ElasticModuli(this%fTempInit)
      end if

    end function Constructor
          


    subroutine Destructor(this)
      type(CrysPlasMat) :: this

      if (allocated(this%fSchmidt))  deallocate(this%fSchmidt)

    end subroutine Destructor



    function GetElasticModuli(this, tempCur) result(Cijkl)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND),  intent(in) :: tempCur
      real(kind = RKIND) :: Cijkl(3, 3, 3, 3)


      if (this%fIsTempDep) then
        Cijkl = this%fCrysElas%ElasticModuli(tempCur)
      else
        Cijkl = this%fCijkl
      end if

    end function GetElasticModuli



    function GetNumSlipSys(this) result(numSlipSys)
      class(CrysPlasMat) :: this
      integer(kind = IKIND) :: numSlipSys

      numSlipSys = this%fNumSlipSys

    end function GetNumSlipSys



    function GetSchmidt(this) result(schmidt)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND) :: schmidt(3, 3, this%fNumSlipSys)
      class(SlipGeom), pointer :: auxSlipGeom

      schmidt = this%fSchmidt

    end function GetSchmidt



    function GetSlipRate(this, tauResl, tauCrit, temp) result(gammaDot)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: gammaDot(this%fNumSlipSys)

      gammaDot = this%fSlipKine%SlipRate(tauResl, tauCrit, temp)

    end function GetSlipRate



    function GetDSlipRateDTauResl(this, tauResl, tauCrit, temp) result(dGammaDotDTauResl)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: dGammaDotDTauResl(this%fNumSlipSys)

      dGammaDotDTauResl = this%fSlipKine%DSlipRateDTauResl(tauResl, tauCrit, temp)

    end function GetDSlipRateDTauResl



    function GetDSlipRateDTauCrit(this, tauResl, tauCrit, temp) result(dGammaDotDTauCrit)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: dGammaDotDTauCrit(this%fNumSlipSys)

      dGammaDotDTauCrit = this%fSlipKine%DSlipRateDTauCrit(tauResl, tauCrit, temp)

    end function GetDSlipRateDTauCrit



    function GetStressDivergenceState(this, tauResl, tauCrit) result(diverged)
      class(CrysPlasMat), intent(in) :: this
      real(kind = RKIND), intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      logical(kind = LKIND) :: diverged

      diverged = this%fSlipKine%StressDivergenceState(tauResl, tauCrit)

    end function GetStressDivergenceState

    


    function GetInitialCRSS(this) result(tauCrit0)
      class(CrysPlasMat) :: this
      real(kind = RKIND) :: tauCrit0(this%fNumSlipSys)

      tauCrit0 = this%fSlipHard%InitialCRSS()

    end function GetInitialCRSS



    function GetTauCritRate(this, gamma, gammaDot, tauCrit, temp) result(tauCritDot)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: tauCritDot(this%fNumSlipSys)

      tauCritDot = this%fSlipHard%TauCritRate(gamma, gammaDot, tauCrit, temp)

    end function GetTauCritRate




    function GetDTauCritDGamma(this, gamma, gammaDot, tauCrit, temp) result(dTauCritDGamma)
      class(CrysPlasMat),  intent(in) :: this
      real(kind = RKIND), intent(in) :: gamma(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: dTauCritDGamma(this%fNumSlipSys, this%fNumSlipSys)

      dTauCritDGamma = this%fSlipHard%DTauCritDGamma(gamma, gammaDot, tauCrit, temp)

    end function GetDTauCritDGamma




    subroutine ParseInput(this, fID)
      use utils, only : ReadLine
      include "crystalElasticity.inc"
      use typeThermoExpa
      use typeThermoExpaIso
      use typeThermoExpaAniso
      include "slipGeometry.inc"
      include "slipKinetics.inc"
      include "slipHardening.inc"

      class(CrysPlasMat)                :: this
      integer(kind = IKIND), intent(in) :: fID

      integer(kind = IKIND) :: crysElasCode
      integer(kind = IKIND) :: thermoExpaCode
      integer(kind = IKIND) :: slipGeomCode
      integer(kind = IKIND) :: slipKineCode
      integer(kind = IKIND) :: slipHardCode
  

      integer(kind = IKIND) :: nSlipSet
      integer(kind = IKIND), allocatable :: theSlipSet(:)

      logical(kind = 4)     :: elasPropTempDep
      integer(kind = IKIND) :: elasPropOrder

      integer(kind = IKIND) :: thermoExpaPropOrder
      real(kind = RKIND)    :: tempRef

      character(len = 256)  :: line
  
  
      ! Crystal elasticity
      call ReadLine(fID, line)
      read(line, *) crysElasCode, elasPropOrder 
      select case (crysElasCode)
      case (kIso)
        this%fCrysElas => CrysElasIso(  fID, this%fIsTempDep, elasPropOrder, this%fTempInit)
      case (kCubic)
        this%fCrysElas => CrysElasCubic(fID, this%fIsTempDep, elasPropOrder, this%fTempInit)
      case (kAniso)
        this%fCrysElas => CrysElasAniso(fID, this%fIsTempDep, elasPropOrder, this%fTempInit)
      case default
        write(*, *) "Unknown Elasticity type ", crysElasCode, " in input file"
        stop
      end select



      ! Thermal  expansion
      call ReadLine(fID, line)
      read(line, *) thermoExpaCode, thermoExpaPropOrder, tempRef
      select case (thermoExpaCode) 
      case (kIsoThExpa)
        this%fThermoExpa => ThermoExpaIso(  fID, thermoExpaPropOrder, tempRef, this%fTempInit)
      case (kAnisoThExpa)
        this%fThermoExpa => ThermoExpaAniso(fID, thermoExpaPropOrder, tempRef, this%fTempInit)
      case default
        write(*, *) "Unknown thermo expansion type ", thermoExpaCode, " in input file"
        stop
      end select
  

      ! Slip Geometroy
      call ReadLine(fID, line)
      read(line, *) slipGeomCode
      select case (slipGeomCode)
      case (kFCC)
        this%fSlipGeom => SlipGeomFCC()
      case (kBCC)
        this%fSlipGeom => SlipGeomBCC()
      case (kHCP)
        this%fSlipGeom => SlipGeomHCP()
      case (kUser)
        this%fSlipGeom => SlipGeomUser(fID)
      case default
        write(*, *) "Unknown Geometry Type ", slipGeomCode, " in input file"
        stop
      end select

      this%fNumSlipSys = this%fSlipGeom%NumSlipSys()
  

      ! Slip Kinetics
      call ReadLine(fID, line)
      read(line, *) slipKineCode
      select case (slipKineCode)
      case (kPowLaw1)
        this%fSlipKine => SlipKinePowLaw1(fID, this%fIsTempDep,  this%fNumSlipSys)
      case default
        write(*, *) "Unknown Slip Kinetics Type ", slipKineCode, " in input file"
        stop
      end select


      ! Slip Hardening
      nSlipSet = this%fSlipGeom%NumSlipSet()
      allocate(theSlipSet(this%fNumSlipSys))
      theSlipSet = this%fSlipGeom%SlipSet()

      call ReadLine(fID, line)
      read(line, *) slipHardCode
      select case (slipHardCode)
      case (kAssaNeed)
        this%fSlipHard => SlipHardAN(fID, this%fIsTempDep, this%fNumSlipSys, nSlipSet, theSlipSet)
      case (kVoceTome)
        this%fSlipHard => SlipHardVT(fID, this%fIsTempDep, this%fNumSlipSys, nSlipSet, theSlipSet)
      case (kVoceKock)
        this%fSlipHard => SlipHardVK(fID, this%fIsTempDep, this%fNumSlipSys, nSlipSet, theSlipSet)
      case default
        write(*, *) "Unknown Slip Hardening Type ", slipHardCode, " in input file"
        stop
      end select

      if (allocated(theSlipSet)) deallocate(theSlipSet)
  

    end subroutine ParseInput

      
  end module typeCrysPlasMat
  
