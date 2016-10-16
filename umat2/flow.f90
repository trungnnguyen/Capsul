  module flow
    use utils, only : RKIND, IKIND, LKIND, MAX_FLOW_PARAS

    implicit none
    private

    public  :: InitFlowParas
    public  :: fGammaDot
    public  :: fDGammaDotDTauRatio
    public  :: fStressDiverged

    integer(kind = IKIND) :: flowType
    real(kind = RKIND)    :: flowProps(MAX_FLOW_PARAS)
 
  contains

    subroutine InitFlowParas(fID)
      use utils, only : ReadLine 

      character(len = 100) :: line

      call ReadLine(fID, line)
      read(line, *) flowType

      if (flowType == 1) then
        call InitFlowParas1(fID)
      else if (flowType == 2) then
        call InitFlowParas2(fID)
      else
        write(*, *) "Error! Unknown Type of Flow Rule Specified!"
        stop 
      end if
      
    end subroutine InitFlowParas


    function fStressDiverged(tauRatio)
      real(kind = RKIND), intent(in) :: tauRatio(nSlipSys)
      logical(kind = LKIND) :: fStressDiverged

      fStressDiverged = any(abs(tauRatio) > mmBig)

    end function fStressDiverged


    function fGammaDot(tauRatio, temp)
      use utils, only : SMALL

      real(kind = RKIND), intent(in) :: tauRatio(nSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: fGammaDot(nSlipSys)

      real(kind = RKIND) :: tauRatioSign(nSlipSys)
      real(kind = RKIND) :: absTauRatio(nSlipSys)
      real(kind = RKIND) :: factor
      integer :: i

    
      if (temp > SMALL) then
        factor = exp(-qActive/temp)
      else
        factor = 1.0d0
      end if

      tauRatioSign = sign(1.0d0, tauRatio)
      absTauRatio  = abs(tauRatio)
      do i = 1, nSlipSys
        fGammaDot(i) = factor*gammaDot0*tauRatioSign(i)*absTauRatio(i)**mmInv
      end do
     
    end function fGammaDot


    ! derivatives of gammaDot with respect to tauCrit in each slip systems
    function fDGammaDotDTauRatio(tauRatio, temp)
      use utils, only : SMALL

      real(kind = RKIND), intent(in) :: tauRatio(nSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: fDGammaDotDTauRatio(nSlipSys)

      real(kind = RKIND) :: absTauRatio(nSlipSys)
      real(kind = RKIND) :: factor
      integer :: i

      if (temp > SMALL) then
        factor = exp(-qActive/temp)
      else
        factor = 1.0d0
      end if

      absTauRatio  = abs(tauRatio)
      do i = 1, nSlipSys
        fDGammaDotDTauRatio(i) = factor*gammaDot0MMInv*absTauRatio(i)**(MMInv - 1)
      end do

    end function fDGammaDotDTauRatio

  end module flow


