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
    type(CrysPlasMatPt)   :: theCrysPlasMatPt 
    real(kind = RKIND)    :: oriArray(6)
    integer(kind = IKIND) :: ansyType 
    real(kind = RKIND)    :: eqvSig, eqvEps

    ansyType = 1
    if (.not. initialized) then
      call initialization(ansyType, temp)
      initialized = .true.
    endif

    call theCrysPlasMatPt%InitMatPt(pCrysPlasMat, ansyType)
    oriArray = props(1:6)
    if (time(2) == 0.0d0) call theCrysPlasMatPt%InitStateVar(oriArray, statev, nstatv)

    call theCrysPlasMatPt%AdvanceStep(dfGrd0, dfGrd1, temp, temp+dtemp, statev, nstatv, dtime)

    stress = theCrysPlasMatPt%GetCauchyStress()
    ddsdde = theCrysPlasMatPt%GetMaterialJacob()
    pNewdt = theCrysPlasMatPt%GetNewDtScaleFactor()
    if (ansyType == 1) then
      ddsddt = theCrysPlasMatPt%GetDSigDTemp()
      rpl    = theCrysPlasMatPt%GetHeatGenRate()
      dRplDT = theCrysPlasMatPt%GetDRplDTemp()
      dRplDE = theCrysPlasMatPt%GetDRplDEps()
    end if
   
!    write(*, *)  "stress = ", stress
!    write(*, *)  "ddsdde = ", ddsdde
!    write(*, *)  "ddsddt = ", ddsddt
!    write(*, *)  "rpl    = ", rpl
!    write(*, *)  "dRplDt = ", drpldt
!    write(*, *)  "dRplDe = ", drplde
    call theCrysPlasMatPt%SaveStateVar(statev, nstatv)

    eqvSig = 1.0d0/sqrt(2.0)*sqrt((stress(1)-stress(2))*(stress(1)-stress(2)) + (stress(2)-stress(3))*(stress(2)-stress(3)) &
                               +(stress(3)-stress(1))*(stress(3)-stress(1)) + 6.0d0*(stress(4)*stress(4) + stress(5)*stress(5) &
                                                                                    +stress(6)*stress(6)))
    
    eqvEps = sqrt(2.0)/3.0d0*sqrt((stran(1)-stran(2))*(stran(1)-stran(2)) + (stran(2)-stran(3))*(stran(2)-stran(3)) &
                               +(stran(3)-stran(1))*(stran(3)-stran(1)) + 6.0d0*(stran(4)*stran(4) + stran(5)*stran(5) &
                                                                                    +stran(6)*stran(6)))
    if (npt == 1) then                                                                            
      write(*, *) KINC, NPT, time(2), pnewDt, eqvSig, eqvEps, rpl, temp
    end if
  end subroutine umat

