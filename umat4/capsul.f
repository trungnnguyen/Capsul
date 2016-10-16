!> @brief  This file contains the subroutines for umat_cp2
  subroutine umat(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                  RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                  NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                  CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    use utils,  only : IKIND, RKIND
    use init,   only : pCrysPlasMatPt, initialized, initialization
    implicit none

    !> Explicit Declaration of interface arguments called by ABAQUS
    include 'umatArgs.inc'
    real(kind = RKIND)    :: oriArray(6)
    integer(kind = IKIND) :: ansyType 

    ansyType = 0
    if (.not. initialized) then
      call initialization(ansyType, temp)
      initialized = .true.
    endif

    oriArray = props(1:6)
    if (time(2) == 0.0d0)  call pCrysPlasMatPt%InitStateVar(oriArray, statev, nstatv)

!    write(*, *) "At Inc ", KINC, " NPT = ", NPT
    call pCrysPlasMatPt%AdvanceStep(dfGrd0, dfGrd1, temp, temp+dtemp, statev, nstatv, dtime)

    stress = pCrysPlasMatPt%GetCauchyStress()
    ddsdde = pCrysPlasMatPt%GetJacob()
    pNewdt = pCrysPlasMatPt%GetNewDtScaleFactor()
    call pCrysPlasMatPt%SaveStateVar(statev, nstatv)

 !   write(*, *) "Stress  = ", stress
 !   write(*, *) "pNewDt  = ", pNewDt

  end subroutine umat

