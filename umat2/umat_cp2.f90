!> @brief  This file contains the subroutines for umat_cp2
  subroutine umat_cp2(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                      RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                      STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                      NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                      CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    use utils, only : IKIND, RKIND, UPDSIG, UPDJAC
    use algebra, only : matInv
    use init,  only : initialized, Initialization
    use expansion, only : fDfGrdTh
    implicit none

    !> Explicit Declaration of interface arguments called by ABAQUS
    include 'umatArgs.inc'

    ! local variables
    real(kind = RKIND)    :: statev0(nstatv)
    integer(kind = IKIND) :: phase

    real(kind = RKIND)    :: temp0
    real(kind = RKIND)    :: temp1
    real(kind = RKIND)    :: dfGrdTh0(3, 3)
    real(kind = RKIND)    :: dfGrdTh1(3, 3)
    real(kind = RKIND)    :: dfGrdMech0(3, 3)
    real(kind = RKIND)    :: dfGrdMech1(3, 3)

    ! Global Initialization. Executed only once.
    if (.not. initialized) then
      call Initialization(props)
      initialized = .true.
    endif
    if (time(2) == 0.0d0)  call InitStatev(statev, nstatv)


    ! Update the Cauchy stress tensor
    phase = UPDSIG
    statev0 = statev
    temp0 = temp
    temp1 = temp + dtemp
    dfGrdTh0 = fDfGrdTh(temp0)
    dfGrdTh1 = fDfGrdTh(temp1)
    dfGrdMech0 = dfGrd0*matInv(dfGrdTh0)
    dfGrdMech1 = dfGrd1*matInv(dfGrdTh1)

    write(*, *) "dfGrdTH0 = ", dfGrdTh0
    write(*, *) "dfGrdTH1 = ", dfGrdTh1
    write(*, *) "dfGrd1   = ", dfGrd1
    write(*, *) "dfGrdMe1 = ", dfGrdMech1

!    dfGrdMech0 = dfGrd0
!    dfGrdMech1 = dfGrd1
    call UpdStress(dfGrdMech0, dfGrdMech1, statev0, statev, nstatv, phase, temp1, dtime,        &
                   stress, ntens, rpl, pNewDt)
               
    if (pNewDt < 1.0d0) return


    ! Update the jacob matrix -- ddsdde
    phase = UPDJAC
    call UpdJacob (dfGrdMech0, dfGrdMech1, statev0, statev, nstatv, phase, temp1, dtime,        & 
                   stress, ntens, ddsdde, pNewDt)


  end subroutine umat_cp2



  subroutine UpdStress(dfGrd0, dfGrd1, statev0, statev, nstatv, phase, temp, dtime,    &
                       stress, ntens, rpl, pNewDt)
    use utils,     only : IKIND,  RKIND, UPDSIG, UPDJAC
    use init,      only : nSlipSys

    implicit none

    real(kind = RKIND),    intent(in)  :: dfGrd0(3, 3)
    real(kind = RKIND),    intent(in)  :: dfGrd1(3, 3)
    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(in)  :: statev0(nstatv)
    real(kind = RKIND),    intent(inout) :: statev(nstatv)
    integer(kind = IKIND), intent(in)  :: phase
    real(kind = RKIND),    intent(in)  :: temp
    real(kind = RKIND),    intent(in)  :: dtime
    integer(kind = IKIND), intent(in)  :: ntens
    real(kind = RKIND),    intent(out) :: stress(ntens)
    real(kind = RKIND),    intent(out) :: rpl
    real(kind = RKIND),    intent(out) :: pNewDt

 
    real(kind = RKIND) :: iterVec(nSlipSys+9)
    real(kind = RKIND) :: workArray(nSlipSys*8+100)


    call InitIterVec(dfGrd0, dfGrd1, statev0, statev, nstatv, phase, temp, dtime, iterVec, workArray)

    call NRSolve(iterVec, workArray, phase, temp, dtime, pNewDt)

    if (pNewDt < 1.0d0) return

    call PostConverge(iterVec, workArray, phase, temp, dtime, stress, ntens, statev, nstatv, rpl)
   

  end subroutine UpdStress



  subroutine InitIterVec(dfGrd0, dfGrd1, statev0, statev1, nstatv, phase, temp, dtime, iterVec, workArray)
    use utils,       only : IKIND,  RKIND, UNITMAT, UPDSIG, UPDJAC
    use algebra,     only : matInv, ten4Rot, polarDcmp
    use init,        only : nSlipSys, oriMatx, schmidt
    use crystal,     only : fStifLoc

    implicit none

    real(kind = RKIND),    intent(in)  :: dfGrd0(3, 3)
    real(kind = RKIND),    intent(in)  :: dfGrd1(3, 3)
    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(in)  :: statev0(nstatv)
    real(kind = RKIND),    intent(in)  :: statev1(nstatv)
    integer(kind = IKIND), intent(in)  :: phase
    real(kind = RKIND),    intent(in)  :: temp
    real(kind = RKIND),    intent(in)  :: dtime
    real(kind = RKIND),    intent(out) :: iterVec(nSlipSys+9)
    real(kind = RKIND),    intent(out) :: workArray(nSlipSys*8+100)


    real(kind = RKIND) ::  dfGrdEls0(3, 3)
    real(kind = RKIND) ::  dfGrdEls1(3, 3)
    real(kind = RKIND) ::  dGamma(nSlipSys)
    real(kind = RKIND) ::  dGamma0(nSlipSys)
    real(kind = RKIND) ::  dGamma1(nSlipSys)
    real(kind = RKIND) ::  gamma(nSlipSys)
    real(kind = RKIND) ::  gamma0(nSlipSys)
    real(kind = RKIND) ::  gamma1(nSlipSys)
    real(kind = RKIND) ::  tauCrit(nSlipSys)
    real(kind = RKIND) ::  tauCrit0(nSlipSys)
    real(kind = RKIND) ::  tauCrit1(nSlipSys)
    real(kind = RKIND) ::  Lp0(3, 3)
    real(kind = RKIND) ::  Lp1(3, 3)

    real(kind = RKIND) ::  dfGrd0Inv(3, 3)
    real(kind = RKIND) ::  dfGrdInc(3, 3)
    real(kind = RKIND) ::  dfGrdEls(3, 3)
    real(kind = RKIND) ::  dfGrdElsPrd(3, 3)
    real(kind = RKIND) ::  dfGrdElsPrdInv(3, 3)
    real(kind = RKIND) ::  dfGrdElsPrdLarDef(3, 3)
    real(kind = RKIND) ::  UU(3, 3)
    real(kind = RKIND) ::  RR(3, 3)
    real(kind = RKIND) ::  stiffLoc(3, 3, 3, 3)
    real(kind = RKIND) ::  stiffGlb(3, 3, 3, 3)


    dfGrd0Inv = matInv(dfGrd0)
    dfGrdInc  = matmul(dfGrd1, dfGrd0Inv)
    call loadStatev(statev0, nstatv, dfGrdEls0, dGamma0, gamma0, tauCrit0, Lp0)
    dfGrdElsPrd    = matmul(dfGrdInc, dfGrdEls0)
    dfGrdElsPrdInv = matInv(dfGrdElsPrd)

    call polarDcmp(dfGrdInc, UU, RR)  
    dfGrdElsPrdLarDef = matmul(RR, dfGrdEls0)
    stiffLoc   = fStifLoc(temp)
    stiffGlb   = ten4Rot(stiffLoc, transpose(oriMatx))

    if (phase == UPDSIG) then
      dfGrdEls = matmul(dfGrdElsPrd, UNITMAT - Lp0*dtime)
      dGamma   = 0.0d0
      gamma    = gamma0
      tauCrit  = tauCrit0
    else if (phase == UPDJAC) then 
      call loadStatev(statev1, nstatv, dfGrdEls1, dGamma1, gamma1, tauCrit1, Lp1)
      dfGrdEls = dfGrdEls1
      dGamma   = dGamma1
      gamma    = gamma1
      tauCrit  = tauCrit1 
    end if

    iterVec(1:  9) = reshape(dfGrdEls, (/9/))
    iterVec(10: 9+nSlipSys) = dGamma

    ! Auxliary constants and variables are stored in workArray for later use
    workArray(1:  81) = reshape(stiffGlb, (/81/))
    workArray(82: 90) = reshape(dfGrdElsPrdInv,    (/9/))
    workArray(91: 99) = reshape(dfGrdElsPrdLarDef, (/9/))

    workArray(101:            100+nSlipSys)   = gamma0
    workArray(101+nSlipSys:   100+nSlipSys*2) = tauCrit0
    workArray(101+nSlipSys*2: 100+nSlipSys*3) = gamma
    workArray(101+nSlipSys*3: 100+nSlipSys*4) = tauCrit

  end subroutine InitIterVec



  subroutine NRSolve(iterVec, workArray,  phase, temp, dtime, pNewDt)
    use utils,     only : RKIND, IKIND, LKIND, UPDSIG, UPDJAC
    use algebra,   only : vecNorm, GaussJordan
    use init,      only : nSlipSys, maxIters, toler, maxItersJac, tolerJac

    implicit none

    real(kind = RKIND),    intent(inout)  :: iterVec(nSlipSys+9)
    real(kind = RKIND),    intent(inout)  :: workArray(nSlipSys*8+100)
    integer(kind = IKIND), intent(in)     :: phase
    real(kind = RKIND),    intent(in)     :: temp
    real(kind = RKIND),    intent(in)     :: dtime
    real(kind = RKIND),    intent(out)    :: pNewDt

    ! local variables
    real(kind = RKIND), parameter :: BIGNORM = 1.0d10

    real(kind = RKIND) :: bigRes(nSlipSys+9)
    real(kind = RKIND) :: bigJac(nSlipSys+9, nSlipSys+9)
    real(kind = RKIND) :: iterVecInc(nSlipSys+9)
    real(kind = RKIND) :: normRes
    real(kind = RKIND) :: normIterVecInc
    real(kind = RKIND) :: iterTol
    integer(kind = IKIND) :: nIters
    integer(kind = IKIND) :: iter

    logical(kind = LKIND) :: processed
    integer(kind = IKIND) :: iflag


    if (phase == UPDSIG) then
      nIters  = maxIters
      iterTol = toler
    else if (phase == UPDJAC) then
      nIters  = maxItersJac
      iterTol = tolerJac
    end if

    iter = 1
    processed = .false.
    do while (iter <= nIters)
      call CalRes(iterVec, temp, dtime, bigRes, workArray, pNewDt)
      if (pNewDt < 1.0d0) return
      normRes = vecNorm(bigRes)
      if (normRes <= iterTol) exit

      if (normRes > BIGNORM .and. .not. processed) then
        call CorOnLargeDef(iterVec, workArray, phase, temp, dtime, pNewDt)
        iter = 0
        processed = .true.
      else
        call CalJAC(iterVec, workArray, temp, dtime, bigJAC)
        call GaussJordan(bigJAC, bigRes, iterVecInc, iflag)
        if (iflag == 1) then
          pNewDt = 0.5d0
          return
        end if
        
        call UpdIterVec(iterVec, iterVecInc, workArray, temp)

        normIterVecInc = vecNorm(iterVecInc)
        if (normIterVecInc > BIGNORM)  call CorOnLargeDef(iterVec, workArray, phase, temp, dtime, pNewDt)

      end if

      iter = iter + 1
    end do

    if (iter >= nIters .and. phase == UPDSIG) then
      write(*, *) "Max iterations reached!"
      pNewDt = 0.75d0
      return
    end if

  end subroutine NRSolve



  subroutine UpdIterVec(iterVec, iterVecInc, workArray, temp) 
    use utils,     only : RKIND
    use hardening, only : fTauCritDot
    use init,      only : nSlipSys

    implicit none

    real(kind = RKIND), intent(inout)  :: iterVec(nSlipSys+9)
    real(kind = RKIND), intent(in)     :: iterVecInc(nSlipSys+9)
    real(kind = RKIND), intent(inout)  :: workArray(nSlipSys*8+100)
    real(kind = RKIND), intent(in)     :: temp

    real(kind = RKIND) :: gamma0(nSlipSys)
    real(kind = RKIND) :: tauCrit0(nSlipSys)

    real(kind = RKIND) :: dGamma(nSlipSys)
    real(kind = RKIND) :: absDGamma(nSlipSys)
    real(kind = RKIND) :: gamma(nSlipSys)
    real(kind = RKIND) :: tauCrit(nSlipSys)
    real(kind = RKIND) :: tauCritInc(nSlipSys)
   


    iterVec = iterVec - iterVecInc

    dGamma    = iterVec(10: 9+nSlipSys)
    absDGamma = dabs(dGamma)

    gamma0    = workArray(101: 100+nSlipSys)
    gamma     = gamma0 + absDGamma
    workArray(101+nSlipSys*2: 100+nSlipSys*3) = gamma

    tauCrit    = workArray(101+nSlipSys*3: 100+nSlipSys*4)
    tauCritInc = fTauCritDot(gamma, absDGamma, tauCrit, temp)
    tauCrit0   = workArray(101+nSlipSys: 100+nSlipSys*2)
    workArray(101+nSlipSys*3 : 100+nSlipSys*4) = tauCrit0 + tauCritInc

  end subroutine UpdIterVec


  ! Calculate the residual (right hand section (RHS) of the equation) in each N-R iteration  
  subroutine CalRes(iterVec, temp, dtime, bigRes, workArray, pNewDt)
    use utils,       only : RKIND, IKIND, UNITMAT
    use algebra,     only : ten4Rot
    use init,        only : nSlipSys, schmidt, oriMatx
    use crystal,     only : fStressDiverged, fStifLoc, fGammaDot, fTauResl, fLp
    use deformation, only : fGreenStrain, fElasStress

    implicit none
    real(kind = RKIND), intent(in)  :: iterVec(nSlipSys+9)
    real(kind = RKIND), intent(in)  :: temp
    real(kind = RKIND), intent(in)  :: dtime
    real(kind = RKIND), intent(out) :: bigRes(nSlipSys+9)
    real(kind = RKIND), intent(out) :: workArray(nSlipSys*8+100)
    real(kind = RKIND), intent(out) :: pNewDt

    real(kind = RKIND) :: stiffGlb(3, 3, 3, 3)

    real(kind = RKIND) :: dfGrdEls(3, 3)
    real(kind = RKIND) :: dGamma(nSlipSys)
    real(kind = RKIND) :: tauCrit(nSlipSys)
    real(kind = RKIND) :: dfGrdElsPrdInv(3, 3)

    real(kind = RKIND) :: LpNew(3, 3)
    real(kind = RKIND) :: LpPrd(3, 3)
    real(kind = RKIND) :: elasStrain(3, 3)
    real(kind = RKIND) :: sKirchRot(3, 3)
    real(kind = RKIND) :: tauResl(nSlipSys)
    real(kind = RKIND) :: tauRatio(nSlipSys)
    real(kind = RKIND) :: tauSign(nSlipSys)
    real(kind = RKIND) :: gammaDot(nSlipSys)

  
    dfGrdEls        = reshape(iterVec(1:9), (/3, 3/))
    dGamma          = iterVec(10: 9+nSlipSys)
    dfGrdElsPrdInv  = reshape(workArray(82:90), (/3, 3/))
    tauCrit         = workArray(101+nSlipSys*3: 100+nSlipSys*4)

    elasStrain  = fGreenStrain(dfGrdEls)
    stiffGlb    = reshape(workArray(1:81), (/3, 3, 3, 3/))
    sKirchRot   = fElasStress(stiffGlb, elasStrain)
    tauResl     = fTauResl(sKirchRot, schmidt(:, :, 1:nSlipSys))
    tauRatio    = dabs(tauResl/tauCrit)
    if (fStressDiverged(tauRatio)) then
      write(*, *) "Diverged resolve stress detected!"
      pNewDt = min(pNewDt, 0.75d0)
      return
    end if

    tauSign   = sign(1.0d0, tauResl)
    gammaDot  = fGammaDot(tauRatio, tauSign, temp)

    LpNew = fLp(gammaDot, schmidt(:, :, 1:nSlipSys))
    LpPrd = (UNITMAT - matmul(dfGrdElsPrdInv, dfGrdEls))/dtime

    bigRes(1 : 9)          = reshape(LpPrd - LpNew, (/9/))
    bigRes(10: 9+nSlipSys) = dGamma - gammaDot*dtime
   
    workArray(101+nSlipSys*4: 100+nSlipSys*5) = tauResl 
    workArray(101+nSlipSys*5: 100+nSlipSys*6) = tauRatio
    workArray(101+nSlipSys*6: 100+nSlipSys*7) = tauSign
    workArray(101+nSlipSys*7: 100+nSlipSys*8) = gammaDot

  end subroutine CalRes



  subroutine CalJAC(iterVec, workArray, temp, dtime, bigJac)
    use utils,       only : RKIND, IKIND
    use algebra,     only : tenmul, tenCtrct44, ten3333ToA99, matNorm
    use init,        only : nSlipSys, schmidt
    use crystal,     only : fGammaDot, fGammaDotDDTauCrit
    use hardening,   only : fDDTauCritDDGamma
    use deformation, only : fPrtlGreenStrain

    implicit none
    real(kind = RKIND), intent(in)  :: iterVec(nSlipSys+9)
    real(kind = RKIND), intent(in)  :: workArray(nSlipSys*8+100)
    real(kind = RKIND), intent(in)  :: temp
    real(kind = RKIND), intent(in)  :: dtime
    real(kind = RKIND), intent(out) :: bigJac(nSlipSys+9, nSlipSys+9)

    real(kind = RKIND) :: stiffGlb(3, 3, 3, 3)
    real(kind = RKIND) :: dfGrdEls(3, 3)
    real(kind = RKIND) :: dfGrdElsPrdInv(3, 3)
    real(kind = RKIND) :: prtlStrain(3, 3, 3, 3)
    real(kind = RKIND) :: tauCrit(nSlipSys)
    real(kind = RKIND) :: gamma(nSlipSys)
    real(kind = RKIND) :: dGamma(nSlipSys)
    real(kind = RKIND) :: tauResl(nSlipSys)
    real(kind = RKIND) :: tauRatio(nSlipSys)
    real(kind = RKIND) :: tauSign(nSlipSys)
    real(kind = RKIND) :: gammaDot(nSlipSys)

    real(kind = RKIND) :: dGammaSign(nSlipSys)
    real(kind = RKIND) :: prtlTauResl(3, 3, nSlipSys)

    real(kind = RKIND) :: dfactor(nSlipSys)
    real(kind = RKIND) :: dfactor2(nSlipSys)

    real(kind = RKIND) :: ddTauCritddGamma(nSlipSys, nSlipSys)
    real(kind = RKIND) :: ddTauCritAux(nSlipSys, nSlipSys)
    real(kind = RKIND) :: aux3333(3, 3, 3, 3), aux3333_2(3, 3, 3, 3)
    real(kind = RKIND) :: aux33(3, 3), auxScalar1
    real(kind = RKIND) :: MJacob(3, 3, 3, 3)
    integer(kind = IKIND) :: i, j, k 


    dfGrdEls = reshape(iterVec(1:9), (/3, 3/))
    dGamma   = iterVec(10: 9+nSlipSys)

    stiffGlb        = reshape(workArray(1:81), (/3, 3, 3, 3/))
    dfGrdElsPrdInv  = reshape(workArray(82:90), (/3, 3/))
    gamma           = workArray(101+nSlipSys*2:  100+nSlipSys*3)
    tauCrit         = workArray(101+nSlipSys*3:  100+nSlipSys*4)
    tauResl         = workArray(101+nSlipSys*4:  100+nSlipSys*5)
    tauRatio        = workArray(101+nSlipSys*5:  100+nSlipSys*6)
    tauSign         = workArray(101+nSlipSys*6:  100+nSlipSys*7)
    gammaDot        = workArray(101+nSlipSys*7:  100+nSlipSys*8)

    dfactor  = fGammaDotDDTauCrit(tauRatio, tauCrit, temp)
    dfactor2 = dfactor*tauRatio*tauSign 
    !BOX11: partial Lp / partial Fe
    MJacob = 0.0d0
    do i = 1, 3
      do j = 1, 3
        auxScalar1 = -dfGrdElsPrdInv(i, j)/dtime
        MJacob(i, 1, j, 1) = auxScalar1
        MJacob(i, 2, j, 2) = auxScalar1
        MJacob(i, 3, j, 3) = auxScalar1
      end do
    end do

    prtlStrain  = fPrtlGreenStrain(dfGrdEls)
    aux3333 = tenCtrct44(stiffGlb, prtlStrain) 
    aux3333_2 = 0.0d0
    do i = 1, nSlipSys
      do j = 1, 3
        do k = 1, 3
          prtlTauResl(j, k, i) = dfactor(i) * sum(schmidt(:, :, i)*aux3333(:, :, j, k))
        end do
      end do
      aux3333_2 = aux3333_2 + tenmul(schmidt(:, :, i), prtlTauResl(:, :, i))
    end do
  
!    jacobNR = MJacob - aux3333_2
    bigJAC(1:9, 1:9) = Ten3333ToA99(MJacob - aux3333_2)

    !BOX12: Partial Lp/ partial dgamma
    ddTauCritDDGamma = fDDTauCritDDGamma(gamma, gammaDot, tauCrit, temp)
    
    dGammaSign = sign(1.0d0, dGamma)
    do i = 1, nSlipSys
      do j = 1, nSlipSys
        ddTauCritAux(i, j) = dfactor2(i)*dGammaSign(j)*ddTauCritDDGamma(i, j)
      end do
    end do

    do i = 1, nSlipSys
      aux33 = 0.0d0
      do j = 1, nSlipSys
        aux33 = aux33 + schmidt(:, :, j)*ddTauCritAux(j, i)
      end do
      bigJac(1:9, 9+i) = reshape(aux33, (/9/))
    end do

    !BOX21: Partial dgamma / partial Fe
    do i = 1, nSlipSys
      bigJac(9+i, 1:9) = reshape(-dtime*prtlTauResl(:, :, i), (/9/))
    end do

    ! BOX22: Partial dgamma/dgamma
    bigJac(10:9+nSlipSys, 10:9+nSlipSys) = ddTauCritAux*dtime
    do i = 1, nSlipSys
      bigJac(9+i, 9+i) = bigJac(9+i, 9+i) + 1.0d0
    end do

  end subroutine CalJAC



  ! NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
  subroutine CorOnLargeDef(iterVec, workArray,  phase, temp, dtime, pNewDt)
    use utils,       only : RKIND, IKIND, UPDSIG, UPDJAC
    use algebra,     only : polarDcmp, ten4Rot
    use init,        only : nSlipSys, oriMatx, schmidt
    use deformation, only : fGreenStrain, fElasStress
    use crystal,     only : fTauResl, fStressDiverged, fGammaDot, fStifLoc
    use hardening,   only : fTauCritDot
    implicit none

    real(kind = RKIND), intent(inout) :: iterVec(nSlipSys+9)
    real(kind = RKIND), intent(inout) :: workArray(nSlipSys*8+100)
    integer(kind = IKIND), intent(in) :: phase
    real(kind = RKIND), intent(in)    :: temp
    real(kind = RKIND), intent(in)    :: dtime
    real(kind = RKIND), intent(inout) :: pNewDt

    real(kind = RKIND) :: stiffGlb(3, 3, 3, 3)

    real(kind = RKIND) :: dfGrdEls0(3, 3)
    real(kind = RKIND) :: gamma0(nSlipSys)
    real(kind = RKIND) :: tauCrit0(nSlipSys)

    real(kind = RKIND) :: dfGrdEls(3, 3)
    real(kind = RKIND) :: elasStrain(3, 3)
    real(kind = RKIND) :: sKirchRot(3, 3)
    real(kind = RKIND) :: tauResl(nSlipSys)
    real(kind = RKIND) :: tauCrit(nSlipSys)
    real(kind = RKIND) :: tauRatio(nSlipSys)
    real(kind = RKIND) :: tauSign(nSlipSys)
    real(kind = RKIND) :: gammaDot(nSlipSys)
    real(kind = RKIND) :: dGamma(nSlipSys)
    real(kind = RKIND) :: absDGamma(nSlipSys)
    real(kind = RKIND) :: gamma(nSlipSys)
    real(kind = RKIND) :: tauCritInc(nSlipSys)


    stiffGlb  = reshape(workArray(1:81), (/3, 3, 3, 3/))
    dfGrdEls  = reshape(workArray(91:99), (/3, 3/))
    gamma0    = workArray(101: 100+nSlipSys)
    tauCrit0  = workArray(101+nSlipSys: 100+nSlipSys*2)

    if (phase == UPDSIG) then
      elasStrain = fGreenStrain(dfGrdEls)
      sKirchRot  = fElasStress(stiffGlb, elasStrain)
      tauResl    = fTauResl(sKirchRot, schmidt(:, :, 1:nSlipSys))
      tauCrit    = workArray(101+nSlipSys*3 : 100+nSlipSys*4)
      tauRatio   = dabs(tauResl/tauCrit)
  
      if (fStressDiverged(tauRatio)) then
        write(*, *) "Diverged tauResl/tauCrit detected!"
        pNewDt = min(pNewDt, 0.75d0)
        return
      end if
  
      tauSign   = sign(1.0d0, tauCrit)
      gammaDot  = fGammaDot(tauRatio, tauSign, temp)
      dGamma    = gammaDot*dtime
      absDGamma = dabs(dGamma)
      gamma     = gamma0 + absDGamma
  
      tauCritInc = fTauCritDot(gamma, absDGamma, tauCrit, temp)
      tauCrit    = tauCrit0 + tauCritInc

    else if (phase == UPDJAC) then
      dGamma  = 0.0d0
      gamma   = gamma0
      tauCrit = tauCrit0
    end if

    iterVec(1:9)           = reshape(dfGrdEls, (/9/))
    iterVec(10:9+nSlipSys) = dGamma

    workArray(101+nSlipSys*2: 100+nSlipSys*3) = gamma
    workArray(101+nSlipSys*3: 100+nSlipSys*4) = tauCrit

  end subroutine CorOnLargeDef
          

!
! ROTATE STRESSES TO DEFORMED CONFIGURATION  
  subroutine PostConverge(iterVec, workArray, phase, temp, dtime, stress, ntens, statev, nstatv, rpl)
    use utils,       only : RKIND, IKIND, UPDSIG, UPDJAC, UNITMAT
    use algebra,     only : polarDcmp, matRot, ten4Rot
    use init,        only : nSlipSys, oriMatx, schmidt
    use crystal,     only : fStifLoc, fTauResl
    use deformation, only : fGreenStrain, fElasStress

    implicit none    

    real(kind = RKIND),    intent(in)  :: iterVec(nSlipSys+9)
    real(kind = RKIND),    intent(in)  :: workArray(nSlipSys*8+100)
    integer(kind = IKIND), intent(in)  :: phase
    real(kind = RKIND),    intent(in)  :: temp
    real(kind = RKIND),    intent(in)  :: dtime
    integer(kind = IKIND), intent(in)  :: ntens
    real(kind = RKIND),    intent(out) :: stress(ntens) 
    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(out) :: statev(nstatv)
    real(kind = RKIND),    intent(out) :: rpl

    real(kind = RKIND) :: stiffGlb(3, 3, 3, 3)

    real(kind = RKIND) :: dfGrdElsNew(3, 3)
    real(kind = RKIND) :: dGammaNew(nSlipSys)
    real(kind = RKIND) :: dfGrdElsPrdInv(3, 3)
    real(kind = RKIND) :: UU(3, 3), RR(3, 3)
    real(kind = RKIND) :: elasStrain(3, 3)
    real(kind = RKIND) :: sKirchRot(3, 3)
    real(kind = RKIND) :: sKirch(3, 3)
    real(kind = RKIND) :: tauCritNew(nSlipSys)
    real(kind = RKIND) :: gammaNew(nSlipSys)
    real(kind = RKIND) :: gammaDot(nSlipSys)
    real(kind = RKIND) :: tauResl(nSlipSys)
    real(kind = RKIND) :: LpNew(3, 3)
    real(kind = RKIND) :: EulerAng(3)

    real(kind = RKIND), parameter :: disFrc = 0.90d0

    integer(kind = RKIND) :: i


    dfGrdElsNew = reshape(iterVec(1:9), (/3, 3/))

    call polarDcmp(dfGrdElsNew, UU, RR)
    elasStrain = fGreenStrain(dfGrdElsNew)
    stiffGlb   = reshape(workArray(1:81), (/3, 3, 3, 3/))
    sKirchRot  = fElasStress(stiffGlb, elasStrain)
    sKirch     = matRot(sKirchRot, transpose(RR))
    
    do i = 1, 3
      stress(i) = sKirch(i, i)
    end do

    stress(4) = sKirch(1, 2)
    stress(5) = sKirch(1, 3)
    stress(6) = sKirch(2, 3)


    if (phase == UPDSIG) then
      tauResl    = fTauResl(sKirchRot, schmidt)
      gammaDot   = workArray(101+nSlipSys*7: 100+nSlipSys*8)
      rpl        = disFrc*sum(tauResl*gammaDot)

      dGammaNew  = iterVec(10:9+nSlipSys)
      gammaNew   = workArray(101+nSlipSys*2: 100+nSlipSys*3)
      tauCritNew = workArray(101+nSlipSys*3: 100+nSlipSys*4) 

      dfGrdElsPrdInv = reshape(workArray(82:90), (/3, 3/))
      LpNew = (UNITMAT - matmul(dfGrdElsPrdInv, dfGrdElsNew))/dtime
      EulerAng   = 0.0d0

      call saveStatev(statev, nstatv, dfGrdElsNew, dGammaNew, gammaNew, tauCritNew, LpNew, EulerAng)

    end if

  end subroutine PostConverge

  

  subroutine UpdJacob(dfGrd0, dfGrd1, statev0, statev1, nstatv, phase, temp, dtime, &
                      stress, ntens, ddsdde, pNewDt)
    use utils,     only : RKIND, IKIND, UNITMAT, UPDJAC
    use init,      only : epsInc, oriMatx
    use crystal,   only : fStifLoc
    use algebra,   only : ten4Rot, ten3333ToA66

    implicit none

    real(kind = RKIND),    intent(in)  :: dfGrd0(3, 3)
    real(kind = RKIND),    intent(in)  :: dfGrd1(3, 3)
    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(in)  :: statev0(nstatv)
    real(kind = RKIND),    intent(in)  :: statev1(nstatv)
    integer(kind = IKIND), intent(in)  :: phase
    real(kind = RKIND),    intent(in)  :: temp
    real(kind = RKIND),    intent(in)  :: dtime
    integer(kind = IKIND), intent(in)  :: ntens
    real(kind = RKIND),    intent(in)  :: stress(ntens)
    real(kind = RKIND),    intent(out) :: ddsdde(ntens, ntens)
    real(kind = RKIND),    intent(out) :: pNewDt


    real(kind = RKIND) :: stressJac(ntens)
    real(kind = RKIND) :: statevJac(nstatv)
    real(kind = RKIND) :: dltEpsJac(3, 3)
    real(kind = RKIND) :: dfGrdIncAux(3, 3)
    real(kind = RKIND) :: dfGrd1Jac(3, 3)
    real(kind = RKIND) :: stiffLoc(3, 3, 3, 3)
    real(kind = RKIND) :: stiffGlb(3, 3, 3, 3)

    real(kind = RKIND)    :: halfEpsInc
    integer(kind = IKIND) :: idxi(6), idxj(6)
    integer(kind = IKIND) :: iStep, idx1, idx2
    real(kind = RKIND)    :: dummyRpl

!   General loop on the 6 perturbations to obtain Jacobian
!   On each istep a whole NR problem is solved to get the corresponding stress
    halfEpsInc = 0.5d0*epsInc
    dltEpsJAC = reshape((/    epsInc, halfEpsInc, halfEpsInc,             &
                          halfEpsInc,     epsInc, halfEpsInc,             &
                          halfEpsInc, halfEpsInc,     epsInc/), (/3, 3/))

    idxi = (/1, 2, 3, 1, 1, 2/)
    idxj = (/1, 2, 3, 2, 3, 3/)
    do iStep = 1, 6
      idx1 = idxi(iStep)
      idx2 = idxj(iStep)
      dfGrdIncAux = UNITMAT
      dfGrdIncAux(idx1, idx2)   = dfGrdIncAux(idx1, idx2) + dltEpsJac(idx1, idx2) 
      if (idx1 < idx2) then
        dfGrdIncAux(idx2, idx1) = dfGrdIncAux(idx1, idx2)
      end if
      dfGrd1Jac = matmul(dfGrdIncAux, dfGrd1)

      statevJac = statev1
      call UpdStress(dfGrd0, dfGrd1Jac, statev0, statevJac, nstatv, phase, temp, dtime,   &
                     stressJac, ntens, dummyRpl, pNewDt) 

      if (pNewDt < 1.0d0) then
        stiffLoc = fStifLoc(temp)
        stiffGlb = ten4Rot(stiffLoc, transpose(oriMatx))
        ddsdde = ten3333ToA66(stiffGlb)
        exit
      else
        ddsdde(:, iStep) = (stressJac - stress)/epsInc
      end if

    end do

           
  end subroutine UpdJacob



  subroutine InitStatev(statev, nstatv)
    use utils,     only : RKIND, IKIND, UNITMAT
    use init,      only : nSlipSys
    use hardening, only : fTauCritInit
    implicit none

    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(out) :: statev(nstatv)

    real(kind = RKIND) :: dfGrdElsInit(3, 3)
    real(kind = RKIND) :: dGammaInit(nSlipSys)
    real(kind = RKIND) :: gammaInit(nSlipSys)
    real(kind = RKIND) :: tauCritInit(nSlipSys)
    real(kind = RKIND) :: LpInit(3, 3)
    real(kind = RKIND) :: EulerAngInit(3)

    dfGrdElsInit  = UNITMAT
    dGammaInit    = 0.0d0
    gammaInit     = 0.0d0
    tauCritInit   = fTauCritInit()
    LpInit        = 0.0d0
    EulerAngInit  = 0.0d0
    
    call saveStatev(statev, nstatv, dfGrdElsInit, dGammaInit, gammaInit, tauCritInit, LpInit, EulerAngInit)

  end subroutine InitStatev



  subroutine loadStatev(statev, nstatv, dfGrdElasCur, dGamma, gammaCur, tauCritCur, LpCur)
    use utils,  only : RKIND, IKIND
    use init,   only : nSlipSys
    implicit none

    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(in)  :: statev(nstatv)
    real(kind = RKIND),    intent(out) :: dfGrdElasCur(3, 3)
    real(kind = RKIND),    intent(out) :: dGamma(nSlipSys)
    real(kind = RKIND),    intent(out) :: gammaCur(nSlipSys)
    real(kind = RKIND),    intent(out) :: tauCritCur(nSlipSys)
    real(kind = RKIND),    intent(out) :: LpCur(3, 3)

    dfGrdElasCur = reshape(statev(1:9), (/3, 3/))
    dGamma       = statev(10: 9+nSlipSys)
    gammaCur     = statev(10+nSlipSys:   9+nSlipSys*2)
    tauCritCur   = statev(10+nSlipSys*2: 9+nSlipSys*3)
    LpCur        = reshape(statev(10+nSlipSys*3: 18+nSlipSys*3), (/3, 3/))

  end subroutine loadStatev



  subroutine saveStatev(statev, nstatv, dfGrdElasNew, dGammaNew, gammaNew, tauCritNew, LpNew, EulerAng)
    use utils,     only : RKIND, IKIND
    use init,      only : nSlipSys
    implicit none

    integer(kind = IKIND), intent(in)  :: nstatv
    real(kind = RKIND),    intent(out) :: statev(nstatv)
    real(kind = RKIND),    intent(in)  :: dfGrdElasNew(3, 3)
    real(kind = RKIND),    intent(in)  :: dGammaNew(nSlipSys)
    real(kind = RKIND),    intent(in)  :: gammaNew(nSlipSys)
    real(kind = RKIND),    intent(in)  :: tauCritNew(nSlipSys)
    real(kind = RKIND),    intent(in)  :: LpNew(3, 3)
    real(kind = RKIND),    intent(in)  :: EulerAng(3)

    integer(kind = IKIND) :: idxLo, idxHi, arrLen


    arrLen = 9
    idxLo = 1
    idxHi = idxLo + arrLen - 1
    statev(idxLo : idxHi)   = reshape(dfGrdElasNew, (/9/))

    arrLen = nSlipSys
    idxLo  = idxHi + 1
    idxHi  = idxLo + arrLen - 1
    statev(idxLo : idxHi) = dGammaNew

    arrLen = nSlipSys
    idxLo  = idxHi + 1
    idxHi  = idxLo + arrLen - 1
    statev(idxLo : idxHi) = gammaNew

    arrLen = nSlipSys
    idxLo  = idxHi + 1
    idxHi  = idxLo + arrLen - 1
    statev(idxLo  :  idxHi)  = tauCritNew

    arrLen = 9
    idxLo  = idxHi + 1
    idxHi  = idxLo + arrLen - 1
    statev(idxLo : idxHi) = reshape(LpNew, (/9/)) 

    arrLen = 3
    idxLo  = idxHi + 1
    idxHi  = idxLo + arrLen - 1
    statev(idxLo: idxHi) = EulerAng

    arrLen = 1
    idxLo  = idxHi + 1
    idxHi  = idxLo + arrLen - 1
    statev(idxLo : idxHi) = sum(gammaNew)

  end subroutine saveStatev

