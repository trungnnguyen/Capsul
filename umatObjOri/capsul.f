!> @brief  This file contains the subroutines for umat_cp2
  subroutine umat(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                  RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                  NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                  CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    use utils,  only : IKIND, RKIND
    use init,   only : pCrysPlasMat, initialized, initialization
    use typeCrysPlasMatPt, only : CrysPlasMatPt 
    implicit none

    !> Explicit Declaration of interface arguments called by ABAQUS
    include 'umatArgs.inc'
    class(CrysPlasMatPt), pointer :: pCrysPlasMatPt => null() 
    real(kind = RKIND)    :: oriArray(6)
    integer(kind = IKIND) :: ansyType 

    ansyType = 0
    if (.not. initialized) then
      call initialization(ansyType, temp)
      initialized = .true.
    endif

    pCrysPlasMatPt => CrysPlasMatPt(pCrysPlasMat)
    oriArray = props(1:6)
    if (time(2) == 0.0d0) call pCrysPlasMatPt%InitStateVar(oriArray, statev, nstatv)

    call pCrysPlasMatPt%AdvanceStep(dfGrd0, dfGrd1, temp, temp+dtemp, statev, nstatv, dtime)

    stress = pCrysPlasMatPt%GetCauchyStress()
    write(*, *) "KINC = ", KINC, " NPT = ", NPT, " STRESS = ", stress
    ddsdde = pCrysPlasMatPt%GetJacob()
    pNewdt = pCrysPlasMatPt%GetNewDtScaleFactor()
    if (ansyType == 1) rpl = pCrysPlasMatPt%GetHeatGenRate()

    call pCrysPlasMatPt%SaveStateVar(statev, nstatv)

  end subroutine umat

