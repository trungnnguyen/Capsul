  module typeCrysPlasMatPt
    use utils, only : RKIND, IKIND, LKIND
    use typeCrysPlasMat, only : CrysPlasMat

    implicit none

    integer(kind = RKIND), parameter :: kUpdSig = 1
    integer(kind = RKIND), parameter :: kUpdJac = 2

    type, public :: CrysPlasMatPt
      real(kind = RKIND)    :: fStress(6) 
      real(kind = RKIND)    :: fDSigDEps(6, 6) 

      class(CrysPlasMat),  pointer :: fCrysPlasMat
      real(kind = RKIND)    :: fDt
      real(kind = RKIND)    :: fTempOld  
      real(kind = RKIND)    :: fFtotOld(3, 3)    

      real(kind = RKIND)    :: fTempNew  
      real(kind = RKIND)    :: fFtotNew(3, 3)  

      real(kind = RKIND)    :: fFElsOld(3, 3)  
      real(kind = RKIND)    :: fLPlsOld(3, 3)  
      real(kind = RKIND), allocatable :: fTauCritOld(:)  
      real(kind = RKIND), allocatable :: fGammaOld(:)  

      real(kind = RKIND)    :: fFElsNew(3, 3)  
      real(kind = RKIND)    :: fLPlsNew(3, 3)  
      real(kind = RKIND), allocatable :: fTauCritNew(:)  
      real(kind = RKIND), allocatable :: fGammaNew(:)  
  
      real(kind = RKIND), allocatable :: fTauResl(:)
      real(kind = RKIND), allocatable :: fGammaDot(:)
      real(kind = RKIND), allocatable :: fDGamma(:)  
      real(kind = RKIND)    :: fRotMatx(3, 3)
      real(kind = RKIND)    :: fEulerAng(3)


      integer(kind = IKIND) :: fNumSlipSys
      real(kind = RKIND)    :: fCijkl(3, 3, 3, 3)
      real(kind = RKIND), allocatable :: fSchmidt(:, :, :)

      real(kind = RKIND)    :: fNewDt

      logical(kind = LKIND) :: fIsTempDep

    contains
      procedure, public  :: AdvanceStep
      procedure, public  :: InitStateVar
      procedure, public  :: LoadStateVar
      procedure, public  :: SaveStateVar
      procedure, public  :: GetCauchyStress
      procedure, public  :: GetJacob
      procedure, public  :: GetNewDtScaleFactor

      procedure, private :: UpdateStress
      procedure, private :: UpdateJacob

      procedure, private :: CauchyStress
      procedure, private :: GreenStrain
      procedure, private :: ElasticStress
      procedure, private :: ResolvedSlipStress
      procedure, private :: VelocityGradientPlastic
      procedure, private :: HeatGenerationRate

      procedure, private :: InitXX
      procedure, private :: FormLHS
      procedure, private :: FormRHS
      procedure, private :: UpdateXX
      procedure, private :: Correction

      procedure, private :: PreStep
      procedure, private :: SolveStep
      procedure, private :: PostStep

      procedure, private :: FormRotMatx
      final              :: Destructor

    end type CrysPlasMatPt


    public
    interface CrysPlasMatPt
      module procedure Constructor
    end interface

  contains
  
    function Constructor(pCrysPlasMat) result(this)
      use utils, only : UNITMAT
      type(CrysPlasMat),   pointer :: pCrysPlasMat
      type(CrysPlasMatPt), pointer :: this
  

      allocate(this)

      this%fCrysPlasMat => pCrysPlasMat

      this%fNumSlipSys  = this%fCrysPlasMat%GetNumSlipSys()

      allocate(this%fTauCritOld(this%fNumSlipSys))  
      allocate(this%fTauCritNew(this%fNumSlipSys))  
      allocate(this%fGammaOld(this%fNumSlipSys))  
      allocate(this%fGammaNew(this%fNumSlipSys))  
  
      this%fFElsOld    = UNITMAT
      this%fLPlsOld    = 0.0d0
      this%fTauCritOld = this%fCrysPlasMat%GetInitialCRSS()
      this%fGammaOld   = 0.0d0
  
      this%fFElsNew    = this%fFElsOld
      this%fLPlsNew    = this%fLPlsOld
      this%fTauCritNew = this%fTauCritOld
      this%fGammaNew   = this%fGammaOld

      this%fTempOld    = 300
      this%fTempNew    = this%fTempOld
      this%fFtotOld    = UNITMAT
      this%fFtotNew    = this%fFtotOld

      this%fDt         = 0.0d0
      this%fNewDt      = 1.0d0
      this%fIsTempDep  = .false.

      this%fRotMatx    = UNITMAT
      this%fEulerAng   = (/0.0d0, 0.0d0, 0.0d0/)

      allocate(this%fTauResl(this%fNumSlipSys))
      this%fTauResl = 0.0d0

      allocate(this%fGammaDot(this%fNumSlipSys))
      this%fGammaDot = 0.0d0

      allocate(this%fDGamma(this%fNumSlipSys))  
      this%fDGamma  = 0.0d0

      allocate(this%fSchmidt(3, 3, this%fNumSlipSys))
      this%fSchmidt = 0.0d0


    end function Constructor
          


    subroutine Destructor(this)
      type(CrysPlasMatPt) :: this

      if (allocated(this%fTauCritOld))  deallocate(this%fTauCritOld)
      if (allocated(this%fTauCritNew))  deallocate(this%fTauCritNew)
      if (allocated(this%fGammaOld))    deallocate(this%fGammaOld)
      if (allocated(this%fGammaNew))    deallocate(this%fGammaNew)
      if (allocated(this%fTauResl))     deallocate(this%fTauResl)
      if (allocated(this%fGammaDot))    deallocate(this%fGammaDot)
      if (allocated(this%fDGamma))      deallocate(this%fDGamma)
      if (allocated(this%fSchmidt))     deallocate(this%fSchmidt)

    end subroutine Destructor



    subroutine FormRotMatx(this, oriArray)
      use algebra, only : normalize, cross_product
  
      class(CrysPlasMatPt)           :: this
      real(kind = RKIND), intent(in) :: oriArray(6)
  
      real(kind = RKIND) :: auxVec1(3)
      real(kind = RKIND) :: auxVec2(3)
  
      auxVec1 = oriArray(1:3)
      auxVec2 = oriArray(4:6)
  
      call normalize(auxVec1)
      call normalize(auxVec2)
  
      this%fRotMatx(:, 1) = auxVec1
      this%fRotMatx(:, 2) = auxVec2
      this%fRotMatx(:, 3) = cross_product(this%fRotMatx(:, 1), this%fRotMatx(:, 2))
  
    end subroutine FormRotMatx
  


    subroutine InitStateVar(this, oriArray, statev, nstatv)
      class(CrysPlasMatPt)               :: this
      real(kind = RKIND),    intent(in)  :: oriArray(6)   
      integer(kind = IKIND), intent(in)  :: nstatv
      real(kind = RKIND),    intent(out) :: statev(nstatv)

      call this%FormRotMatx(oriArray)
      call this%SaveStateVar(statev, nstatv)

    end subroutine InitStateVar



    subroutine LoadStateVar(this, statev, nstatv)
      class(CrysPlasMatPt)              :: this
      integer(kind = IKIND), intent(in) :: nstatv
      real(kind = RKIND),    intent(in) :: statev(nstatv)

      integer(kind = IKIND) :: arrLen, idxLo, idxHi


      arrLen = 9
      idxLo = 1
      idxHi = idxLo + arrLen - 1
      this%fFElsOld = reshape(statev(idxLo : idxHi), (/3, 3/))  
      this%fFElsNew = this%fFElsOld
  
    
      arrLen = this%fNumSlipSys
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      this%fGammaOld = statev(idxLo : idxHi)
      this%fGammaNew = this%fGammaOld
  
      arrLen = this%fNumSlipSys
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      this%fTauCritOld = statev(idxLo : idxHi)
      this%fTauCritNew = this%fTauCritOld
  
      arrLen = 9
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      this%fLPlsOld = reshape(statev(idxLo : idxHi), (/3, 3/))
      this%fLPlsNew = this%fLPlsOld

      arrLen = 9
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      this%fRotMatx = reshape(statev(idxLo: idxHi), (/3, 3/))
  

      arrLen = 3
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      this%fEulerAng = statev(idxLo: idxHi)

    end subroutine LoadStateVar


  
    subroutine SaveStateVar(this, statev, nstatv)
      class(crysPlasMatPt),  intent(in)  :: this
      integer(kind = IKIND), intent(in)  :: nstatv
      real(kind = RKIND),    intent(out) :: statev(nstatv)
  
      integer(kind = IKIND) :: idxLo, idxHi, arrLen
  
      arrLen = 9
      idxLo = 1
      idxHi = idxLo + arrLen - 1
      statev(idxLo : idxHi)   = reshape(this%fFElsNew, (/9/))
  
        
      arrLen = this%fNumSlipSys
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      statev(idxLo : idxHi) = this%fGammaNew
  
      arrLen = this%fNumSlipSys
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      statev(idxLo  :  idxHi)  = this%fTauCritNew
  
      arrLen = 9
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      statev(idxLo : idxHi) = reshape(this%fLPlsNew, (/9/)) 

      arrLen = 9
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      statev(idxLo : idxHi) = reshape(this%fRotMatx, (/9/)) 
  

      arrLen = 3
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      statev(idxLo: idxHi) = this%fEulerAng
  
      arrLen = 1
      idxLo  = idxHi + 1
      idxHi  = idxLo + arrLen - 1
      statev(idxLo : idxHi) = sum(this%fGammaNew)
  
    end subroutine SaveStateVar
  

  
            
    subroutine AdvanceStep(this, FTotOld, FTotNew, TempOld, TempNew, statev, nstatv, dtime)
      use algebra, only : Ten3333ToA66
      class(CrysPlasMatPt) :: this
      real(kind = RKIND),    intent(in) :: FTotOld(3, 3)
      real(kind = RKIND),    intent(in) :: FTotNew(3, 3)
      real(kind = RKIND),    intent(in) :: TempOld
      real(kind = RKIND),    intent(in) :: TempNew
      integer(kind = IKIND), intent(in) :: nstatv
      real(kind = RKIND),    intent(in) :: statev(nstatv)
      real(kind = RKIND),    intent(in) :: dtime

      real(kind = RKIND) :: stressNew(6)
      real(kind = RKIND) :: statevNew(nstatv)

      integer(kind = IKIND) :: stage


      stage     = kUpdSig
      call this%UpdateStress(FTotOld, FTotNew, TempOld, TempNew, statev, nstatv, dtime, stage)
      if (this%fNewDt < 1.0d0) return

      stressNew = this%fStress
      write(*, *) "stress = ", this%fStress 
      call this%SaveStateVar(statevNew, nstatv)

      stage = kUpdJac
      call this%UpdateJacob(FTotOld, FTotNew, TempOld, TempNew, statev, nstatv, dtime, stage)
      this%fStress = stressNew
      call this%LoadStateVar(statevNew, nstatv)

    end subroutine AdvanceStep



    function GetCauchyStress(this) result(sigma)
      class(CrysPlasMatPt) :: this
      real(kind = RKIND)   :: sigma(6)

      sigma = this%fStress

    end function GetCauchyStress



    function GetJacob(this) result(jacob)
      class(CrysPlasMatPt) :: this
      real(kind = RKIND)   :: jacob(6, 6)

      jacob = this%fDSigDEps

    end function GetJacob



    function GetNewDtScaleFactor(this) result(pNewDt)
      class(CrysPlasMatPt) :: this
      real(kind = RKIND)   :: pNewDt

      pNewDt = this%fNewDt

    end function GetNewDtScaleFactor



    subroutine UpdateStress(this, FTotOld, FTotNew, TempOld, TempNew, statev, nstatv, dtime, stage)
      class(CrysPlasMatPt) :: this
      real(kind = RKIND),    intent(in) :: FTotOld(3, 3)
      real(kind = RKIND),    intent(in) :: FTotNew(3, 3)
      real(kind = RKIND),    intent(in) :: TempOld
      real(kind = RKIND),    intent(in) :: TempNew
      integer(kind = IKIND), intent(in) :: nstatv
      real(kind = RKIND),    intent(in) :: statev(nstatv)
      real(kind = RKIND),    intent(in) :: dtime
      integer(kind = IKIND), intent(in) :: stage

      real(kind = RKIND)    :: statevBak(nstatv)
      integer :: i


      this%fFTotOld = FTotOld
      this%fFTotNew = FTotNew
      this%fTempOld = TempOld
      this%fTempNew = TempNew
      this%fDt      = dtime

      statevBak = statev
      if (stage == kUpdSig) then
      call this%LoadStateVar(statevBak, nstatv) 
      end if

      call this%PreStep(stage)
      call this%SolveStep(stage)
      call this%PostStep(stage)

    end subroutine UpdateStress




    subroutine UpdateJacob(this, FTotOld, FTotNew, TempOld, TempNew, statev, nstatv, dtime, stage)
      use utils,   only : UNITMAT
      use algebra, only : Ten3333ToA66

      class(CrysPlasMatPt) :: this
      real(kind = RKIND),    intent(in) :: FTotOld(3, 3)
      real(kind = RKIND),    intent(in) :: FTotNew(3, 3)
      real(kind = RKIND),    intent(in) :: TempOld
      real(kind = RKIND),    intent(in) :: TempNew
      integer(kind = IKIND), intent(in) :: nstatv
      real(kind = RKIND),    intent(in) :: statev(nstatv)
      real(kind = RKIND),    intent(in) :: dtime
      integer(kind = IKIND), intent(in) :: stage

      real(kind = RKIND) :: FTotInc(3, 3)
      real(kind = RKIND) :: FTotJac(3, 3)
      real(kind = RKIND) :: statevBak(nstatv)

      real(kind = RKIND), parameter :: epsInc = 5.0d-5
      real(kind = RKIND)    :: halfEpsInc
      real(kind = RKIND)    :: dltEpsJac(3, 3)
      real(kind = RKIND)    :: sigConvg(6)
      real(kind = RKIND)    :: sigJac(6)
      integer(kind = IKIND) :: idx1, idx2, idxi(6), idxj(6), iComp


      if (this%fNewDt < 1.0d0) return

      sigConvg = this%fStress

      ! General loop on the 6 perturbations to obtain Jacobian
      halfEpsInc = 0.5d0*epsInc
      dltEpsJac = reshape((/  epsInc, halfEpsInc, halfEpsInc,             &
                          halfEpsInc,     epsInc, halfEpsInc,             &
                          halfEpsInc, halfEpsInc,     epsInc/), (/3, 3/))

      idxi = (/1, 2, 3, 1, 1, 2/)
      idxj = (/1, 2, 3, 2, 3, 3/)
      do iComp = 1, 6
        idx1 = idxi(iComp)
        idx2 = idxj(iComp)
        FTotInc = UNITMAT
        FTotInc(idx1, idx2) = FTotInc(idx1, idx2) + dltEpsJac(idx1, idx2) 
        if (idx1 /= idx2) FtotInc(idx2, idx1) = FtotInc(idx1, idx2)

        FTotJac = matmul(FTotInc, FTotNew)
        call this%UpdateStress(FTotOld, FTotJac, tempOld, tempNew, statev, nstatv, dtime, stage)
        sigJac = this%fStress

        if (this%fNewDt < 1.0d0) then
          this%fDSigDEps = ten3333ToA66(this%fCijkl)
          exit
        else
          this%fDSigDEps(:, iComp) = (sigJac - sigConvg)/epsInc
        end if

      end do


    end subroutine UpdateJacob




    subroutine PreStep(this, stage)
      use algebra, only : ten4Rot, multQBQt

      class(CrysPlasMatPt) :: this
      integer(kind = IKIND), intent(in) :: stage

      real(kind = RKIND) :: Cijkl(3, 3, 3, 3)
      real(kind = RKIND) :: Schmidt(3, 3, this%fNumSlipSys)

      integer(kind = IKIND) :: i


      if (stage == kUpdJac) return

      Cijkl       = this%fCrysPlasMat%GetElasticModuli(this%fTempNew)
      this%fCijkl = ten4Rot(Cijkl, this%fRotMatx)

      Schmidt     = this%fCrysPlasMat%GetSchmidt()
      do i = 1, this%fNumSlipSys
        this%fSchmidt(:, :, i) = multQBQt(Schmidt(:, :, i), this%fRotMatx)
      end do


    end subroutine PreStep




    subroutine SolveStep(this, stage)
      use algebra, only : vecNorm, GaussJordan

      class(CrysPlasMatPt) :: this
      integer(kind = IKIND), intent(in) :: stage
  
      real(kind = RKIND), parameter :: kBigNorm = 1.0d10
  
      real(kind = RKIND) :: LHS(9+this%fNumSlipSys, 9+this%fNumSlipSys)
      real(kind = RKIND) :: RHS(9+this%fNumSlipSys)
      real(kind = RKIND) :: XX(9+this%fNumSlipSys)
      real(kind = RKIND) :: DX(9+this%fNumSlipSys)
      real(kind = RKIND) :: auxIter(18)

      real(kind = RKIND)    :: normRHS
      real(kind = RKIND)    :: normDX
      real(kind = RKIND)    :: pNewDt
      real(kind = RKIND)    :: iterTol
      integer(kind = IKIND) :: nIters
      integer(kind = IKIND) :: iter
  
      logical(kind = LKIND) :: processed
      integer(kind = IKIND) :: iflag
  

      integer(kind = IKIND), parameter :: maxItersSig = 100
      integer(kind = IKIND), parameter :: maxItersJac = 40
      real(kind    = RKIND), parameter :: tolerSig    = 1.0d-10
      real(kind    = RKIND), parameter :: tolerJac    = 1.0d-11


      call this%InitXX(stage, XX, auxIter)
  
      if (stage == kUpdSig) then
        nIters  = maxItersSig
        iterTol = tolerSig
      else if (stage == kUpdJac) then
        nIters  = maxItersJac
        iterTol = tolerJac
      end if
  
      pNewdt    = 1.0d0
      iter      = 1
      processed = .false.
      do while (iter <= nIters)
        write(*, *) "FEls   = ", this%fFElsNew
        write(*, *) "dGamma = ", this%fDGamma
        call this%FormRHS(auxIter, RHS, pNewDt)
        normRHS = vecNorm(RHS)
        write(*, *) "XX = ", XX(1:15)
        write(*, *) "iter = ", iter, " RHS = ", RHS(1:15)
        if (normRHS <= iterTol) exit
        if (normRHS > kBigNorm .and. .not. processed) then
          call this%Correction(stage, auxIter, XX, pNewDt)
          iter      = 0
          processed = .true.
         ! write(*, *) "Correction 1"
        else
          call this%FormLHS(auxIter, LHS)
          call GaussJordan(LHS, RHS, DX, iflag)
          write(*, *) "iter = ", iter, " DX = ", DX(1:15)
          if (iflag == 1)  pNewDt = min(pNewDt, 0.5d0)
          
          normDX = vecNorm(DX)
          call this%UpdateXX(XX, DX)
          if (normDX > kBigNorm)  then
            call this%Correction(stage, auxIter, XX, pNewDt)
          end if
  
        end if

        iter = iter + 1

        this%fNewDt = min(this%fNewDt, pNewDt)
        if (this%fNewDt < 1.0d0) return 

      end do
  
      if (iter >= nIters .and. stage == kUpdSig) then
        write(*, *) "Max iterations reached!"
        pNewDt      = 0.75d0
        this%fNewDt = min(this%fNewDt, pNewDt)
      end if
  
    end subroutine SolveStep



    subroutine PostStep(this, stage)
      use algebra,   only : multQBQt, matDet
      class(CrysPlasMatPt) :: this
      integer(kind = IKIND), intent(in) :: stage

      real(kind = RKIND) :: elsGrn(3, 3)
      real(kind = RKIND) :: sigPK2(3, 3)
      real(kind = RKIND) :: sigCauchy(3, 3)
      real(kind = RKIND) :: det

      integer(kind = IKIND) :: i, idxi(6), idxj(6)


      if (this%fNewDt < 1.0d0) return
      elsGrn    = this%GreenStrain(this%fFElsNew)
      sigPK2    = this%ElasticStress(elsGrn)
      sigCauchy = multQBQt(sigPK2, this%fFElsNew)
      det       = matDet(this%fFElsNew)
      sigCauchy = sigCauchy/det

      idxi = (/1, 2, 3, 1, 1, 2/)
      idxj = (/1, 2, 3, 2, 3, 3/)
      do i = 1, 6
        this%fStress(i) = sigCauchy(idxi(i), idxj(i))
      end do

    end subroutine PostStep
      



    subroutine InitXX(this, stage, XX, auxIter)
      use utils,   only : UNITMAT
      use algebra, only : matInv, polarDcmp

      class(CrysPlasMatPt) :: this
      integer(kind = IKIND), intent(in) :: stage
      real(kind = RKIND)   :: XX(9+this%fNumSlipSys)
      real(kind = RKIND)   :: auxIter(18)

      real(kind = RKIND) :: FTotInc(3, 3)
      real(kind = RKIND) :: FElsPrd(3, 3)
      real(kind = RKIND) :: RR(3, 3)
      real(kind = RKIND) :: UU(3, 3)


      FTotInc  = matmul(this%fFTotNew, matInv(this%fFTotOld))
      FElsPrd  = matmul(FTotInc, this%fFElsOld) 

      if (stage == kUpdSig) then
        this%fFElsNew = matmul(FElsPrd, UNITMAT - this%fLPlsOld*this%fDt)
        this%fDGamma  = 0.0d0
        this%f
      end if
      

      XX(1:  9) = reshape(this%fFElsNew, (/9/))
      XX(10: 9+this%fNumSlipSys) = this%fDGamma

     
      auxIter(1:9) = reshape(matInv(FElsPrd), (/9/))

      call polarDcmp(FTotInc, RR, UU)
      auxIter(10:18) = reshape(matmul(RR, this%fFTotOld), (/9/))


    end subroutine InitXX

   
      

    ! Calculate the residual (right hand section(RHS) of the equation) in each N-R iteration  
    subroutine FormRHS(this, auxIter, RHS, pNewDt)
      use utils,   only : UNITMAT
      use algebra, only : vecNorm
  
      class(CrysPlasMatPt) :: this
      real(kind = RKIND), intent(in)  :: auxIter(18)
      real(kind = RKIND), intent(out) :: RHS(9+this%fNumSlipSys)
      real(kind = RKIND)              :: pNewDt
  
      real(kind = RKIND) :: elsEGrn(3, 3)
      real(kind = RKIND) :: sigPK2(3, 3)
      real(kind = RKIND) :: LPlsCur(3, 3)
      real(kind = RKIND) :: FElsPrdInv(3, 3)
  
      elsEGrn        = this%GreenStrain(this%fFElsNew)
      sigPK2         = this%ElasticStress(elsEGrn)
      this%fTauResl  = this%ResolvedSlipStress(sigPK2)
      !write(*, *) "tauResl = ", this%fTauResl
      if (this%fCrysPlasMat%GetStressDivergenceState(this%fTauResl, this%fTauCritNew)) then
        write(*, *) "Diverged tauResl/tauCrit detected!"
        pNewDt = min(pNewDt, 0.75d0)
        return
      end if

      this%fGammaDot = this%fCrysPlasMat%GetSlipRate(this%fTauResl, this%fTauCritNew, this%fTempNew)
      LPlsCur  = this%VelocityGradientPlastic(this%fGammaDot)


      write(*, *) "FElsNew = ", this%fFElsNew
      write(*, *) "tauResl = ", this%fTauResl
      write(*, *) "tauCrit = ", this%fTauCritNew
      write(*, *) "gamaDot = ", this%fGammaDot

      FElsPrdInv    = reshape(auxIter(1:9), (/3, 3/))
      this%fLPlsNew = (UNITMAT - matmul(FElsPrdInv, this%fFElsNew))/this%fDt
  
      write(*, *) "LpPrd = ", FElsPrdInv
      write(*, *) "LpCur = ", LPlsCur
      RHS(1: 9)  = reshape(this%fLPlsNew - LPlsCur, (/9/))
      RHS(10: 9+this%fNumSlipSys) = this%fDGamma - this%fGammaDot*this%fDt
     

    end subroutine FormRHS



    subroutine FormLHS(this, auxIter, LHS)
      use utils,   only : UNITMAT
      use algebra, only : tenmul, MultAijmnBmnkl, ten3333ToA99

      class(CrysPlasMatPt) :: this
      real(kind = RKIND), intent(in) :: auxIter(18)
      real(kind = RKIND) :: LHS(9+this%fNumSlipSys, 9+this%fNumSlipSys)
  
      real(kind = RKIND) :: prtlStrain(3, 3, 3, 3)
      real(kind = RKIND) :: prtlTauResl(3, 3, this%fNumSlipSys)
  
      real(kind = RKIND) :: dGammaDotDTauResl(this%fNumSlipSys)
      real(kind = RKIND) :: dGammaDotDTauCrit(this%fNumSlipSys)
  
      real(kind = RKIND) :: dTauCritDGamma(this%fNumSlipSys, this%fNumSlipSys)
      real(kind = RKIND) :: dTauCritAux(this%fNumSlipSys, this%fNumSlipSys)
      real(kind = RKIND) :: FElsPrdInv(3, 3)
      real(kind = RKIND) :: auxTen4th1(3, 3, 3, 3), auxTen4th2(3, 3, 3, 3), auxTen4th3(3, 3, 3, 3)
      real(kind = RKIND) :: auxMatx(3, 3), aux1
      integer(kind = IKIND) :: nSlipSys
      integer(kind = IKIND) :: i, j, k, l
  
  
      nSlipSys = this%fNumSlipSys
  
      dGammaDotDTauResl = this%fCrysPlasMat%GetDSlipRateDTauResl(this%fTauResl, this%fTauCritNew, this%fTempNew)
      dGammaDotDTauCrit = this%fCrysPlasMat%GetDSlipRateDTauCrit(this%fTauResl, this%fTauCritNew, this%fTempNew)
      write(*, *) "dd1 = ", dGammaDotDTauResl
      write(*, *) "dd2 = ", dGammaDotDTauCrit
      !BOX11: partial Lp / partial Fe
      auxTen4th1 = 0.0d0
      FElsPrdInv = reshape(auxIter(1:9), (/3, 3/))
      do i = 1, 3
        do j = 1, 3
          aux1 = -FElsPrdInv(i, j)/this%fDt
          auxTen4th1(i, 1, j, 1) = aux1
          auxTen4th1(i, 2, j, 2) = aux1
          auxTen4th1(i, 3, j, 3) = aux1
        end do
      end do
  
      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              prtlStrain(i, j, k, l) = 0.5d0*(UNITMAT(i, l)*this%fFElsNew(k, j) + UNITMAT(j, l)*this%fFElsNew(k, i))
            end do
          end do
        end do
      end do
  
      auxTen4th2 = MultAijmnBmnkl(this%fCijkl, prtlStrain) 
      auxTen4th3 = 0.0d0
      do i = 1, nSlipSys
        do j = 1, 3
          do k = 1, 3
            prtlTauResl(j, k, i) = dGammaDotDTauResl(i) * sum(this%fSchmidt(:, :, i)*auxTen4th2(:, :, j, k))
          end do
        end do
        auxTen4th3 = auxTen4th3 + tenmul(this%fSchmidt(:, :, i), prtlTauResl(:, :, i))
      end do
    
      LHS(1:9, 1:9) = Ten3333ToA99(auxTen4th1 - auxTen4th3)
  
      !BOX12: Partial Lp/ partial dgamma
      dTauCritDGamma = this%fCrysPlasMat%GetDTauCritDGamma(this%fGammaNew, this%fGammaDot, this%fTauCritNew, this%fTempNew)
      do i = 1, nSlipSys
        do j = 1, nSlipSys
          dTauCritAux(i, j) = dGammaDotDTauCrit(i)*dTauCritDGamma(i, j)
        end do
      end do
  
      do i = 1, nSlipSys
        auxMatx = 0.0d0
        do j = 1, nSlipSys
          auxMatx = auxMatx + this%fSchmidt(:, :, j)*dTauCritAux(j, i)
        end do
        LHS(1:9, 9+i) = reshape(auxMatx, (/9/))
      end do
  
      !BOX21: Partial dgamma / partial Fe
      do i = 1, nSlipSys
        LHS(9+i, 1:9) = reshape(-this%fDt*prtlTauResl(:, :, i), (/9/))
      end do
  
      ! BOX22: Partial dgamma/dgamma
      LHS(10:9+nSlipSys, 10:9+nSlipSys) = dTauCritAux*this%fDt
      do i = 1, nSlipSys
        LHS(9+i, 9+i) = LHS(9+i, 9+i) + 1.0d0
      end do
  
    end subroutine FormLHS
  
  


    subroutine UpdateXX(this, XX, dX) 
      class(CrysPlasMatPt) :: this
      real(kind = RKIND)             :: XX(9+this%fNumSlipSys)
      real(kind = RKIND), intent(in) :: dX(9+this%fNumSlipSys)
  
      real(kind = RKIND) :: dTauCrit(this%fNumSlipSys)
  
  
      XX = XX - dX
  
      this%fFElsNew  = reshape(XX(1:9), (/3, 3/))
      this%fDGamma   = XX(10: 9+this%fNumSlipSys)

      this%fGammaNew = this%fGammaOld + dabs(this%fDGamma)
  
      dTauCrit = this%fCrysPlasMat%GetTauCritRate(this%fGammaNew, this%fDGamma, this%fTauCritNew, this%fTempNew)
      this%fTauCritNew = this%fTauCritOld + dTauCrit
  
    end subroutine UpdateXX
  


    ! NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
    subroutine Correction(this, stage, auxIter, XX, pNewDt)
      class(crysPlasMatPt) :: this
      integer(kind = IKIND), intent(in) :: stage
      real(kind = RKIND),    intent(in) :: auxIter(18)
      real(kind = RKIND)                :: XX(9+this%fNumSlipSys)
      real(kind = RKIND)                :: pNewDt
  
      real(kind = RKIND) :: elsEGrn(3, 3)
      real(kind = RKIND) :: sigPK2(3, 3)
      real(kind = RKIND) :: tauCritDot(this%fNumSlipSys)
  
 
      this%fFElsNew =  reshape(auxIter(10:18), (/3, 3/))
      if (stage == kUpdSig) then
        elsEGrn        = this%GreenStrain(this%fFElsNew)
        sigPK2         = this%ElasticStress(elsEGrn)
        this%fTauResl  = this%ResolvedSlipStress(sigPK2)
  
        if (this%fCrysPlasMat%GetStressDivergenceState(this%fTauResl, this%fTauCritNew)) then
          write(*, *) "Diverged tauResl/tauCrit detected!"
          pNewDt = min(pNewDt, 0.75d0)
          return
        end if
    
        this%fGammaDot = this%fCrysPlasMat%GetSlipRate(this%fTauResl, this%fTauCritNew, this%fTempNew)
        this%fDGamma   = this%fGammaDot*this%fDt
        this%fGammaNew = this%fGammaOld + dabs(this%fDGamma)
    
        tauCritDot = this%fCrysPlasMat%GetTauCritRate(this%fGammaNew, this%fGammaDot, this%fTauCritNew, this%fTempNew)
        this%fTauCritNew = this%fTauCritOld + tauCritDot*this%fDt
  
      else if (stage == kUpdJac) then
        this%fDGamma     = 0.0d0
        this%fGammaNew   = this%fGammaOld + dabs(this%fDGamma)
        this%fTauCritNew = this%fTauCritOld
      end if

      XX(1:  9) = reshape(this%fFElsNew, (/9/))
      XX(10: 9+this%fNumSlipSys) = this%fDGamma

    end subroutine Correction




    function CauchyStress(this, FElsNew) result(sigma)
      use algebra, only : multQBQt, matDet

      class(CrysPlasMatPt), intent(in) :: this
      real(kind = RKIND),   intent(in) :: FElsNew(3, 3)
      real(kind = RKIND)               :: sigma(6)

      real(kind = RKIND) :: elsGrn(3, 3)
      real(kind = RKIND) :: sigPK2(3, 3)
      real(kind = RKIND) :: sigCauchy(3, 3)
      real(kind = RKIND) :: det

      integer(kind = IKIND) :: i, idxi(6), idxj(6)


      elsGrn = this%GreenStrain(FElsNew)
      sigPK2 = this%ElasticStress(elsGrn)
      sigCauchy = multQBQt(sigPK2, FElsNew)
      det       = matDet(FElsNew)
      sigCauchy = sigCauchy/det

      idxi = (/1, 2, 3, 1, 1, 2/)
      idxj = (/1, 2, 3, 2, 3, 3/)
      do i = 1, 6
        sigma(i) = sigCauchy(idxi(i), idxj(i))
      end do

    end function CauchyStress


      

    function GreenStrain(this, FF) result(epsGrn)
      use utils, only : UNITMAT

      class(CrysPlasMatPt), intent(in) :: this
      real(kind = RKIND),   intent(in) :: FF(3, 3)
      real(kind = RKIND) :: epsGrn(3, 3)

      real(kind = RKIND) :: CC(3, 3)

      CC     = matmul(transpose(FF), FF)
      epsGrn = 0.5d0*(CC - UNITMAT) 

    end function GreenStrain



    function ElasticStress(this, epsGrn) result(sigma)
      use algebra, only : multCijklEkl

      class(CrysPlasMatPt), intent(in) :: this
      real(kind = RKIND),   intent(in) :: epsGrn(3, 3)
      real(kind = RKIND) :: sigma(3, 3)

      sigma = multCijklEkl(this%fCijkl, epsGrn)

    end function ElasticStress



    function ResolvedSlipStress(this, sigma) result(tauResl)
      class(CrysPlasMatPt), intent(in) :: this
      real(kind = RKIND),   intent(in) :: sigma(3, 3)
      real(kind = RKIND)               :: tauResl(this%fNumSlipSys)

      real(kind = RKIND)    :: auxMatx(3, 3)
      integer(kind = IKIND) :: i

      do i = 1, this%fNumSlipSys
        auxMatx = this%fSchmidt(:, :, i)
        tauResl(i) = sum(auxMatx*sigma)
      end do

    end function



        

    function VelocityGradientPlastic(this, gammaDot) result(Lp)
      class(CrysPlasMatPt), intent(in) :: this
      real(kind = RKIND),   intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND)               :: Lp(3, 3)

      integer(kind = IKIND) :: i

      Lp(1:3, 1:3) = 0.0d0
      do i = 1, this%fNumSlipSys
        Lp = Lp + gammaDot(i)*this%fSchmidt(:, :, i)
      end do

    end function VelocityGradientPlastic




    function HeatGenerationRate(this, gammaDot, tauResl) result(RPL)
      class(CrysPlasMatPt), intent(in) :: this
      real(kind = RKIND),   intent(in) :: gammaDot(this%fNumSlipSys)
      real(kind = RKIND),   intent(in) :: tauResl(this%fNumSlipSys)
      real(kind = RKIND)               :: RPL

      integer(kind = IKIND) :: i

      RPL = 0.0d0
      do i = 1, this%fNumSlipSys
        RPL = RPL + gammaDot(i)*tauResl(i)
      end do

    end function HeatGenerationRate


      
  end module typeCrysPlasMatPt

