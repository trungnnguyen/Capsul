  module hardening
    use utils, only : RKIND, IKIND, MAX_SLIP_SYSTEMS

    implicit none

    private

    public InitHardParas, fTauCritInit, fTauCritDot, fDDTauCritDDGamma

    type ANHard
      real(kind = RKIND) :: q(MAX_SLIP_SYSTEMS, MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: h0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: tau0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: taus(MAX_SLIP_SYSTEMS)
    end type ANHard

    type VTHard
      real(kind = RKIND) :: q(MAX_SLIP_SYSTEMS, MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: h0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: h1(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: tau0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: taus(MAX_SLIP_SYSTEMS)
    end type VTHard

    type VKHard
      real(kind = RKIND) :: h0
      real(kind = RKIND) :: tau0
      real(kind = RKIND) :: taus0
      real(kind = RKIND) :: gammas0
      real(kind = RKIND) :: ma
    end type VKHard


    integer(kind = IKIND) :: nSlipSys

    character(len = 16) :: hardType
    type(ANHard)  :: assaNeed
    type(VTHard)  :: voceTome
    type(VKHard)  :: voceKock

    real(kind = RKIND)  :: aux1(MAX_SLIP_SYSTEMS+5)

  contains
  

  
    subroutine InitHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      use utils, only : ReadLine, toLower
  
      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: nSlipSystems
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSystems)

      character(len = 100) :: line
      

      call ReadLine(fID, line)
      read(line, *) hardType
      call toLower(hardType)

      if (hardType(1:8) == "assaneed") then
        call InitANHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      else if (hardType(1:8) == "vocetome") then
        call InitVTHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      else if (hardType(1:8) == "vocekock") then
        call InitVKHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      end if

    end subroutine InitHardParas

    

    subroutine InitANHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: nSlipSystems
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSystems)

      character(len = 100) :: line

      real(kind = RKIND) :: qsys(6, 6)
      real(kind = RKIND) :: h0_set(6), tau0_set(6), taus_set(6)

      integer(kind = IKIND) :: ii, jj, iSet1, iSet2
      integer(kind = IKIND) :: stat1, stat2, stat3, stat4, stat5


      nSlipSys = nSlipSystems

      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) (qsys(ii, jj), jj = 1, nSlipSets)
      end do
   
      do ii = 1, nSlipSys
        do jj = 1, nSlipSys
          iSet1 = systemSet(ii)
          iSet2 = systemSet(jj)
          assaNeed%q(ii, jj) = qsys(iSet1, iSet2)
          if(ii == jj) assaNeed%q(ii,jj) = 1.0d0
        end do
      end do
      
      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) tau0_set(ii), taus_set(ii), h0_set(ii)
      end do

      do ii = 1, nSlipSys
        iSet1 = systemSet(ii)
        assaNeed%tau0(ii) = tau0_set(iSet1)
        assaNeed%taus(ii) = taus_set(iSet1)
        assaNeed%h0(ii)   = h0_set(iSet1)
      end do

      aux1(1:nSlipSys) = assaNeed%h0(1:nSlipSys) / (assaNeed%taus(1:nSlipSys) - assaNeed%tau0(1:nSlipSys))
     
    end subroutine InitANHardParas



    subroutine InitVTHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: nSlipSystems
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSystems)

      character(len = 100) :: line

      real(kind = RKIND) :: qsys(6, 6)
      real(kind = RKIND) :: h0_set(6), h1_set(6), tau0_set(6), taus_set(6)

      integer(kind = IKIND) :: ii, jj, iSet1, iSet2
      integer(kind = IKIND) :: stat1, stat2, stat3, stat4, stat5, stat6


      nSlipSys = nSlipSystems

      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) (qsys(ii, jj), jj = 1, nSlipSets)
      end do

      do ii = 1, nSlipSys
        do jj = 1, nSlipSys
          iSet1 = systemSet(ii)
          iSet2 = systemSet(jj)
          voceTome%q(ii, jj) = qsys(iSet1, iSet2)
          if(ii == jj) voceTome%q(ii,jj) = 1.0d0
        end do
      end do

      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) tau0_set(ii), taus_set(ii), h0_set(ii), h1_set(ii)
      end do

      do ii = 1, nSlipSys
        iset1 = systemSet(ii)
        voceTome%tau0(ii) = tau0_set(iSet1)
        voceTome%taus(ii) = taus_set(iSet1) - tau0_set(iSet1)
        voceTome%h0(ii)   = h0_set(iSet1)
        voceTome%h1(ii)   = h1_set(iSet1)
      end do

      if(any(voceTome%taus(1:nSlipSys) < 0.0d0)) then
        write(*, *) 'ERROR IN HARDENING', voceTome%taus(1:nSlipSys)
        stop
      end if

      aux1(1:nSlipSys) = voceTome%h0(1:nSlipSys) / voceTome%taus(1:nSlipSys)

    end subroutine InitVTHardParas



    subroutine InitVKHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: nSlipSystems
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSystems)

      character(len = 100) :: line


      nSlipSys = nSlipSystems

      call readLine(fID, line)
      read(line, *) voceKock%h0, voceKock%tau0, voceKock%taus0,  &
                    voceKock%gammas0, voceKock%ma

    end subroutine InitVKHardParas



    function fTauCritInit()
      real(kind = RKIND) :: fTauCritInit(nSlipSys)

      if (hardType(1:8) == "assaneed") then
        fTauCritInit = assaNeed%tau0(1:nSlipSys)
      else if (hardType(1:8) == "vocetome") then
        fTauCritInit = voceTome%tau0(1:nSlipSys)
      else if (hardType(1:8) == "vocekock") then
        fTauCritInit = voceKock%tau0
      end if

    end function fTauCritInit




    function fTauCritDot(gamma, gammaDot, tauCrit, temp)
      
      real(kind = RKIND), intent(in)  :: gamma(nSlipSys)
      real(kind = RKIND), intent(in)  :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in)  :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in)  :: temp

      real(kind = RKIND) :: fTauCritDot(nSlipSys)

      if (hardType(1:8) == "assaneed") then
        fTauCritDot = fTauCritDotAN(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocetome") then
        fTauCritDot = fTauCritDotVT(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocekock") then
        fTauCritDot = fTauCritDotVK(gammaDot, tauCrit, temp)
      end if

    end function fTauCritDot 


    !> \[ \frac{\partial \dot{\tau_i}}{\partial \dot{\gamma_j}} \]
    function fDDTauCritDDGamma(gamma, gammaDot, tauCrit, temp)
      
      real(kind = RKIND), intent(in)  :: gamma(nSlipSys)
      real(kind = RKIND), intent(in)  :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in)  :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in)  :: temp

      real(kind = RKIND) :: fDDTauCritDDGamma(nSlipSys, nSlipSys)

      if (hardType(1:8) == "assaneed") then
        fDDTauCritDDGamma = fDDTauCritDDGammaAN(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocetome") then
        fDDTauCritDDGamma = fDDTauCritDDGammaVT(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocekock") then
        fDDTauCritDDGamma = fDDTauCritDDGammaVK(gammaDot, tauCrit, temp)
      end if

    end function fDDTauCritDDGamma 



    function fTauCritDotAN(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fTauCritDotAN(nSlipSys)

      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: aux11(nSlipSys)
      real(kind = RKIND) :: hh(nSlipSys)
      real(kind = RKIND) :: qq(nSlipSys, nSlipSys)


      gammaTot = sum(gamma)
      hh = hAN(gammaTot)
      qq = assaNeed%q(1:nSlipSys, 1:nSlipSys)

      aux11 = reshape(matmul(qq, gammaDot), (/nSlipSys/))
      fTauCritDotAN = aux11 * hh

    end function fTauCritDotAN


      
    function fDDTauCritDDGammaAN(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND)  :: fDDTauCritDDGammaAN(nSlipSys, nSlipSys)

      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: hh(nSlipSys)
      real(kind = RKIND) :: qq(nSlipSys, nSlipSys)
      integer(kind = IKIND) :: i, j


      gammaTot = sum(gamma)
      hh = hAN(gammaTot)
      qq = assaNeed%q(1:nSlipSys, 1:nSlipSys)

      do i = 1, nSlipSys
        do j = 1, nSlipSys
          fDDTauCritDDGammaAN(i, j) = qq(i, j)*hh(i)
        end do
      end do

    end function fDDTauCritDDGammaAN


    function hAN(gammaTot)
      use algebra, only : secah

      real(kind = RKIND), intent(in) :: gammaTot

      real(kind = RKIND) :: hAN(nSlipSys)

      real(kind = RKIND) :: aux12(nSlipSys)


      aux12 = secah(aux1(1:nSlipSys) * gammaTot)
      hAN = assaNeed%h0(1:nSlipSys)*aux12*aux12

    end function hAN



    function fTauCritDotVT(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fTauCritDotVT(nSlipSys)

      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: aux11(nSlipSys)
      real(kind = RKIND) :: hh(nSlipSys)
      real(kind = RKIND) :: qq(nSlipSys, nSlipSys)


      gammaTot = sum(gamma)

      hh = hVT(gammaTot)
      qq = voceTome%q(1:nSlipSys, 1:nSlipSys)

      aux11 = reshape(matmul(qq, gammaDot), (/nSlipSys/))
      fTauCritDotVT = aux11 * hh

    end function fTauCritDotVT



    function fDDTauCritDDGammaVT(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fDDTauCritDDGammaVT(nSlipSys, nSlipSys)

      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: hh(nSlipSys)
      real(kind = RKIND) :: qq(nSlipSys, nSlipSys)
      integer(kind = IKIND) :: i, j


      gammaTot = sum(gamma)

      hh = hVT(gammaTot)
      qq = voceTome%q(1:nSlipSys, 1:nSlipSys)

      do i = 1, nSlipSys
        do j = 1, nSlipSys
          fDDTauCritDDGammaVT(i, j) = qq(i, j)*hh(i)
        end do
      end do

    end function fDDTauCritDDGammaVT


    function hVT(gammaTot)

      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: hVT(nSlipSys)
 
      real(kind = RKIND) :: expo(nSlipSys)
 
      expo  = dexp(-gammaTot * aux1(1:nSlipSys))
      hVT =  voceTome%h1(1:nSlipSys)*(1.0d0 - expo)                             &
                + (voceTome%taus(1:nSlipSys) + voceTome%h1(1:nSlipSys)*gammaTot)      &
                * aux1(1:nSlipSys)*expo 
 
    end function hVT




    function fTauCritDotVK(gammaDot, tauCrit, temp)

      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fTauCritDotVK(nSlipSys)

      real(kind = RKIND) :: mm
      real(kind = RKIND) :: taus
      real(kind = RKIND) :: gammaDotSum, factor
      integer(kind = IKIND) :: i


      mm = temp/voceKock%ma

      gammaDotSum = sum(abs(gammaDot))
      taus = voceKock%taus0 * (gammaDotSum/voceKock%gammas0) ** mm
      factor = voceKock%h0/(taus - voceKock%tau0)*gammaDotSum
      do i = 1, nSlipSys
        fTauCritDotVK(i) = factor*(taus - tauCrit(i))
      end do

    end function fTauCritDotVK



    function fDDTauCritDDGammaVK(gammaDot, tauCrit, temp)

      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fDDTauCritDDGammaVK(nSlipSys, nSlipSys) ! dTauCrit(i)/dGamma(j)

      real(kind = RKIND) :: mm
      real(kind = RKIND) :: taus
      real(kind = RKIND) :: gammaDotSum, factor, tmp1, tmp2, tmp3, tmp4
      real(kind = RKIND) :: gammaDotSign(nSlipSys)
      integer(kind = IKIND) :: i, j


      mm = temp/voceKock%ma
      
      gammaDotSign = sign(1.0d0, gammaDot)
      gammaDotSum = sum(abs(gammaDot))
      taus = voceKock%taus0 * (gammaDotSum/voceKock%gammas0) ** mm
      tmp1 = taus - voceKock%tau0
      tmp2 = -voceKock%h0*mm*taus
      factor = voceKock%h0/tmp1
      do i = 1, nSlipSys
        tmp3 = factor*(taus - tauCrit(i))
        tmp4 = tmp2*(voceKock%tau0 - tauCrit(i)) / (tmp1*tmp1)
        do j = 1, nSlipSys
          fDDTauCritDDGammaVK(i, j) = (tmp3 + tmp4)*gammaDotSign(j)
        end do
      end do

    end function fDDTauCritDDGammaVK


  end module hardening


