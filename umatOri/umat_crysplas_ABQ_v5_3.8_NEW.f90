!**************************************************************************************************
!                                   UMAT CRYSTAL PLASTICITY 
!     LARGE DEFORMATIONS AND ROTATIONS IN TOTAL DEFORMATIONS
!     FROM BOOK OF SOUZA,PERIC,OWEN
!   Author:
!       Javier Segurado
!   Version:
!       V5.3.8, Ene 2015
!       Modified on 7th of May to include in the STATEV the euler angles
!   History:
!       Subversion 5.2- 14th June 2012
!            -- EXP MAP and DEXP MAP are suppressed (linear approach)
!            -- Includes more general symmetry of stiffness matrix
!            -- write of messages supressed 
!            -- Works under parallel in all ABQ versions (6.7-6.11)
!            -- Now the former file subroutines_v4.f is included in file
!            -- orientate_tensor4 updated to use explicit rotation formulae from mathematica
!       Subversion 5.3-  18th June 2012
!                   -- Voce Hardening law included
!       Subversion 5.3.1-18th July 2012
!                   -- Improvement of NR jacobian
!       Subversion 5.3.2-June 2013
!                   -- Pseudo-line-search to improove numerical efficiency
!       Subversion 5.3.8. January 2014
!                   -- Implicit in Fe and tau by a New residual formulation on Lp & Delta_gamma_i
!                   -- Exact Jacobian and Quadratic convergency for new residual. Robust 
!
!**************************************************************************************************


  subroutine UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                  RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                  NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                  CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)


    include 'ABA_PARAM.INC'  

    implicit none

    DIMENSiON STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


!   Declaration of interface variables
    integer(kind = IKIND) ::  
    real(kind = RKIND)    :: STRESS(NTENS)
    real(kind = RKIND)    :: STATEV(NSTATV)
    real(kind = RKIND)    :: DDSDDE(NTENS, NTENS)
    real(kind = RKIND)    :: DDSDDT(NTENS)
    real(kind = RKIND)    :: DRPLDE(NTENS)
    real(kind = RKIND)    :: DRPLDT(NTENS)
    real(kind = RKIND)    :: DRPLDT(NTENS)
    real(kind = RKIND)    :: DRPLDT(NTENS)
    real(kind = RKIND)    :: DRPLDT(NTENS)
    real(kind = RKIND)    :: DRPLDT(NTENS)



!   To set a maximum of slip systems
    integer max_nsystems, nsets, nsystems
    parameter (max_nsystems = 30)

!     FILE OPENING
    character(len = 2)   ::  FLABEL
    character(len = 80)  ::  outdir, filename
    integer(len = IKIND) ::  lenoutdir,nfiles,ns,UR0

!   MATERIAL PARAMETERS

!   elastic
    real(kind = RKIND) ::  c11, c12, c44, c13, c33, c66
!   viscoplastic
    real(kind = RKIND) ::  mm, gamma_0, tautothemm, logtau, taufrac
!     hardening law
    real(kind = RKIND) ::  tau0(max_nsystems), taus(max_nsystems)
    real(kind = RKIND) ::  h0(max_nsystems), h
!     Voce hardening law
    real(kind = RKIND) ::  h1(max_nsystems)
!   latent hardening matrix
    real(kind = RKIND) ::  q(max_nsystems, max_nsystems)
!   normal planes and tangential planes
    real(kind = RKIND) ::  m(max_nsystems, 3), s(max_nsystems, 3)
    real(kind = RKIND) ::  dnorm_m, dnorm_s
!   From inputs the definition of schmidt matrices
    real(kind = RKIND) ::  schmidt(max_nsystems, 3, 3)
!   tau_max
    real(kind = RKIND) ::  tau_max

   
    real(kind = RKIND) ::  dfactor(max_nsystems), dfactor2(max_nsystems)
    real(kind = RKIND) ::  stiff(3,3,3,3),stiff_orient(3,3,3,3)
    real(kind = RKIND) ::  stiff_rotated(3,3,3,3)

 
    real(kind = RKIND) ::  dfgrd_inc(3, 3)
    real(kind = RKIND) ::  dfgrd0_inv(3, 3)

    real(kind = RKIND) ::  dfgrd_elas(3, 3)
    real(kind = RKIND) ::  dfgrd_elas_t(3, 3)
    real(kind = RKIND) ::  dfgrd_elas_pred(3, 3),     dfgrd_elas_pred0(3, 3)
    real(kind = RKIND) ::  dfgrd_elas_pred_inv(3, 3), dfgrd_elas_pred0_inv(3, 3)
    real(kind = RKIND) ::  epsilon_el(3, 3)
   
    
    real(kind = RKIND)  :: RR_tot(3,3), RR_tot_tr(3,3)
    real(kind = RKIND)  :: RR(3,3), RR_tr(3,3), UU(3,3), VV(3,3), CC(3,3)
   

    real(kind = RKIND)  :: skirchoff(3, 3),jota
    real(kind = RKIND)  :: skirchoff_rot(3, 3)
    
    real(kind = RKIND)  :: L_p(3,3), L_p_new(3,3), EXP_L_p(3,3), Lterm, term
   
    real(kind = RKIND)  :: gamma_pred(max_nsystems)
    real(kind = RKIND)  :: gamma_dot(max_nsystems),gamma_act(max_nsystems)
    real(kind = RKIND)  :: gamma_dot_old(max_nsystems),gama_inc(max_nsystems)
    real(kind = RKIND)  :: gamma_dot_res(max_nsystems),error_hard
    real(kind = RKIND)  :: gamma_TOT,gamma_tot_pred,gamma_tot_act
    real(kind = RKIND)  :: tau_act(max_nsystems),tau_pred(max_nsystems)
    real(kind = RKIND)  :: tau_pred0(max_nsystems),tau_predn(max_nsystems)
    real(kind = RKIND)  :: tau_resolved(max_nsystems)
   

    real(kind = RKIND)  :: orient1(3),orient2(3),orient3(3)
    real(kind = RKIND)  :: orient_MAT(3,3),orient_MAT_T(3,3)
    real(kind = RKIND)  :: ROTATED_ORIENT(3,3),ROTATED_ORIENT_T(3,3)
    real(kind = RKIND)  :: orient_MAT_X(3,3),rotation_axisX
    real(kind = RKIND)  :: rotation_X(3,3)
    real(kind = RKIND)  :: vfixed(3,1),vred(3,1),vred_orient(3,1)
    
    real(kind = RKIND)  :: aux_escalar
    real(kind = RKIND)  :: aux3_1(3),aux3_2(3)
    real(kind = RKIND)  :: aux33(3,3),aux33_2(3,3),aux33_3(3,3),aux33_4(3,3)
    real(kind = RKIND)  :: sqrt2,sqrt3

    real(kind = RKIND)  :: tensor4(3,3,3,3)
    real(kind = RKIND)  :: JACOB_NR(3,3,3,3),jacob_NR_i(3,3,3,3)
    real(kind = RKIND)  :: AUX66(6,6),aux99(9,9)
    real(kind = RKIND)  :: aux3333(3,3,3,3),aux3333_2(3,3,3,3)
    real(kind = RKIND)  :: aux3333j(3,3,3,3)
    real(kind = RKIND)  :: dexpx(3,3,3,3)
    real(kind = RKIND)  :: RES(3,3),dnorm,dnorm_NR,dnorm_old,toler,toler_jac
    real(kind = RKIND)  :: dnorm_hard
    integer(kind = IKIND) ::  nincmax, nincmax_jac
    
    

    integer i,j,ii,jj,kk,ll,pp,qq,nm,nn,iii,jjj
    integer ia,ib,ic,id,ii1,jj1
    integer isystem,jsystem
    logical noconv
    integer index, iter, intern, kroneker, istep, mod
    integer implicit_hard,iflag
    integer intv(2)

!   for jacobian
    real(kind = RKIND)  :: EXP_L_P2(3,3)
    real(kind = RKIND)  :: DFGRD_ELAS_tinv(3,3)
    real(kind = RKIND)  :: PARTIALFEF(3,3,3,3),PARTIALFFE(3,3,3,3)
   
    real(kind = RKIND)  :: JACOB3333(3,3,3,3)
    real(kind = RKIND)  :: sum_s_m_dG(3,3,3,3)
    real(kind = RKIND)  :: dnorm_inc1,dnorm_inc0
    real(kind = RKIND)  :: skirchoff2(3,3),delta_eps_jac(3,3)
    real(kind = RKIND)  :: dtime2,strain_increment
    real(kind = RKIND)  :: gamma_act_jacob(max_nsystems)
    real(kind = RKIND)  :: gamma_pred_jacob(max_nsystems)
    real(kind = RKIND)  :: tau_pred_jacob(max_nsystems),tau_act_jacob(max_nsystems)
    real(kind = RKIND)  :: dfgrd_elas_pred_jacob(3,3)
    real(kind = RKIND)  :: dfgrd_inc0(3,3)
    real(kind = RKIND)  :: tinterpol,detF

    real(kind = RKIND)  :: BIGJAC(max_nsystems+9,max_nsystems+9)
    real(kind = RKIND)  :: BIGRES(max_nsystems+9)
    real(kind = RKIND)  :: BIGCORR(max_nsystems+9),R1_FE(9,9),deter
    integer(kind = IKIND) ::  ntot,nloc
    real(kind = RKIND)  :: dgamma(max_nsystems), dgamma_new(max_nsystems)
    real(kind = RKIND)  :: dgamma_jacob(max_nsystems),signog,signogT
    real(kind = RKIND)  :: mjacob(3,3,3,3),partial_sigma(3,3,3,3)
    real(kind = RKIND)  :: partial_tau(max_nsystems,3,3)
    real(kind = RKIND)  :: partial_tau2(3,3)

!    real(kind = RKIND)  :: pi
    
!   Variables for damping NR
    integer(kind = IKIND) ::  toolarge

    logical init

!   Common Variables saved once at the beginning
    save init
    save mm,gamma_0,tau0,taus,h0,h1,tau_max
    save q,stiff    
    save s,m
    save nsystems
    save pi,sqrt2,sqrt3,tensor4
    save toler,toler_jac,nincmax,nincmax_jac,strain_increment
    save implicit_hard

!   CALCULATION OF CRYSTAL PROPERTIES, ONLY ONCE AT BEGINNING
!   Reads 'crystal.prop' and saves: 
!   mm ,gamma_0,tau0,taus,h0,h1
!   q,schmidt,stiff
!   orient_mat,stiff_orient 
!   nsystems
    
!   To prevent problems in multithread when reading properties
!   Note that element 1010101 have to exist!
    if(.not. init ) then
!      init=.true.
       if(noel /= 1010101) then
          call sleep(1)
       endif
    endif

    if(.not. init ) then


!   Some constants
!       pi = 4d0*datan(1d0)         
!       sqrt2 = 1/dsqrt(2D0)
!       sqrt3 = 1/dsqrt(3D0)
!   Unit 4 order tensor

       do ii = 1, 3
          do jj = 1, 3
             do kk = 1, 3
                do ll = 1, 3
                   tensor4(ii,jj,kk,ll) = 0.0d0
                   if(ii == kk .AND. jj == ll) tensor4(ii,jj,kk,ll) = 1.0d0
                enddo
             enddo
          enddo
       enddo
       
!   OPEN THE FILE WITH THE CRYSTAL DEFINITION
       CALL GETOUTDIR(OUTDIR, LENOUTDIR)
       FILENAME = OUTDIR(1:LENOUTDIR)//'/crystal.prop'
!       write(*,*) filename
       UR0 = 74
       OPEN(UR0, FILE = FILENAME, STATUS = 'OLD')

       CALL READPROPERTIES(UR0, max_nsystems, c11, c12, c44, c13, c33, c66,       &
                           gamma_0,mm,nsets,nsystems,s,m,q,tau0,taus,h0,h1,       &
                           toler,toler_jac,nincmax,nincmax_jac,strain_increment,
                           implicit_hard)
       
       CLOSE(unit = UR0)

!   Generation of second order stiffness matrix STIFF

       if(abs(c33) > 1.0d-10) then
          CALL STIFF6(c11, c12, c44, c13, c33, c66, stiff)
       else
          CALL STIFF4(c11, c12, c44, stiff)
       endif      
       
       tau_max = 0d0
       do ii = 1, nsystems
          tau_max = max(tau_max, tau0(ii))
       end do

       init = .true.

    end if

!   DEFINITION OF THE ORIENTATION MATRIX 

    orient1(1) = props(1)
    orient1(2) = props(2)
    orient1(3) = props(3)
    call norm(dnorm, orient1, 3, 1)
    do i = 1, 3
       orient1(i) = orient1(i)/dnorm
    end do

    orient2(1) = props(4)
    orient2(2) = props(5)
    orient2(3) = props(6)
    call norm(dnorm, orient2, 3, 1)
    do i = 1, 3
       orient2(i) = orient2(i)/dnorm
    end do

    call pvect(orient3, orient1, orient2)
    call norm(dnorm, orient3, 3, 1)
    do i = 1, 3
       orient3(i) = orient3(i)/dnorm
    end do

    do i = 1, 3
       orient_MAT(i, 1) = orient1(i)
       orient_MAT(i, 2) = orient2(i)
       orient_MAT(i, 3) = orient3(i)
    end do
   
    orient_MAT_T = transpose(orient_mat)
!    call transpose(orient_MAT_T, orient_mat, 3, 3)
    
!   DEFINITION OF THE SCHMDIT MATRICES 
    
    do ii = 1, nsystems
       do i = 1, 3
          aux3_1(i) = s(ii, i)
          aux3_2(i) = m(ii, i)
       end do
       CALL PTENS(aux33, aux3_1, aux3_2)
       CALL PMAT(aux33_2,ORIENT_MAT,aux33,3,3,3)
       CALL PMAT(aux33,aux33_2,ORIENT_MAT_T,3,3,3)

       do i = 1, 3
          do j = 1, 3
             schmidt(ii,i,j) = aux33(i, j)
          end do 
       end do
   
    end do         
      
!   Orientation of stiffness matrix to material orientation, STIFF_ORIENT
  
    CALL ORIENTATE_TENSOR4(STIFF_ORIENT, STIFF, ORIENT_MAT_T)

!   READ STATE VARIABLES:
  
!   If time=0 initialize variables

    if(time(2) == 0.0d0) then
      
       call dzeros(DFGRD_ELAS,3,3)
       call dzeros(DFGRD_ELAS_t,3,3)
       call dzeros(l_p,3,3)
       
       do ii = 1, 3
          DFGRD_ELAS(ii, ii)  = 1d0
          DFGRD_ELAS_t(ii,ii) = 1d0            
       end do
       
       call dzeros(gamma_act,max_nsystems,1)
       call dzeros(tau_act,max_nsystems,1)

       do ii = 1,nsystems
          tau_act(ii) = tau0(ii)
       end do

    else
       index = 0
       do ii = 1, 3
          do jj = 1, 3
             index = index + 1
             DFGRD_ELAS_t(ii, jj) = statev(index)
          end do
       end do

       do i = 1, nsystems
          index = index + 1
          gamma_act(i) = statev(index)
       end do

       do i = 1, nsystems
          index = index + 1
          tau_act(i) = statev(index)
       end do

       do ii = 1,3
          do jj = 1,3
             index = index + 1
             L_p(ii,jj) = statev(index)
          end do
       end do
    end if

!   PREDICTOR:

!   Deformation gradient increment

    CALL MATINV3(DFGRD0_inv, DFGRD0, iflag)
    if(iflag /= 0) WRITE(*,*)'ERROR INVERTIGN DFGRD0_inv'

    CALL PMAT(DFGRD_INC,DFGRD1,DFGRD0_inv,3,3,3)   

!   Elastic deformation gradient predictor, DGGRD_ELAS_pred0
    
    CALL PMAT(DFGRD_ELAS_pred0,DFGRD_INC,DFGRD_ELAS_t,3,3,3)  

!   Elastic deformation gradient prediction based on last L_p
       
    do ii = 1, 3
       do jj = 1, 3
          L_p_new(ii,jj) = L_p(ii,jj)*(-dtime)
       enddo
    enddo
       
!   CALL EXPMAP(EXP_L_P,noconv,L_P)
    do ii = 1,3
       do jj = 1, 3
          EXP_L_P(ii, jj) = kroneker(ii, jj) + L_p_new(ii, jj)
       enddo
    enddo
    
    CALL PMAT(DFGRD_ELAS_PRED,dfgrd_elas_pred0,exp_l_p,3,3,3)
          
               
!   gamma and CRSS initialization

    call dzeros(dgamma_new, max_nsystems, 1)
    call dzeros(dgamma,     max_nsystems, 1)

    gamma_TOT_pred = 0d0
    do isystem = 1, nsystems
       tau_pred(isystem)   = tau_act(isystem)
       gamma_pred(isystem) = gamma_act(isystem)
       gamma_TOT_pred      = gamma_TOT_pred + gamma_pred(isystem)
       gamma_TOT_act       = gamma_TOT_act  + gamma_act(isystem)
!   dgamma_prediction = 0

       dgamma(isystem) = 0.

    enddo
              
!   GLOBAL NEWTON-RAPHSON

    dnorm_NR = 1d20
    iter     = 0    
    toolarge = 0
    
    DO 100 WHILE(dnorm_NR > TOLER)
      
       iter = iter + 1


!   1.1: Decomposition of Fe in U and R   
!   1.3: Lagrangian deformation

!   Lagrangian deformation
    CALL GREEN_LAGRANGE(DFGRD_ELAS_pred,EPSILON_EL)
   
!   1.4: Calculation of Kirchoff stress

!   skirchoff_rot DEFINED ON crystal axes (undeformed configuration)
!   skirchoff DEFINED ON final configuration

    CALL PCONTRACT2(skirchoff_rot, STIFF_ORIENT, epsilon_el, 3)
    
!   2: PROYECTION OF STRESS IN SYSTEMS

!   2.2: Obtention of resolved tau on all systems

    do isystem = 1, nsystems
       do ii = 1, 3
          do jj = 1,3
             AUX33_2(ii,jj) = schmidt(isystem,ii,jj)
          enddo
       enddo
       CALL PCONTRACT(tau_resolved(isystem),skirchoff_rot,AUX33_2,3)

    enddo

!   3: CALCULATION OF SLIP RATES


!cUEVO ALGORITMO PARA tau_pred
    
    do isystem = 1, nsystems
       
       CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,                    &
                     tau_resolved(isystem),tau_pred(isystem), noconv)
       
       if(noconv .EQV. .TRUE.) then
          write(*,*) 'ERROR IN GAMMADOT',isystem
          pnewdt = .75
          return             
       endif
       
       dgamma_new(isystem) = gamma_dot(isystem)*dtime
                 
    enddo
       

!    write(*,*) 'after Viscolaw'

!   4: CALCULATION OF FP

!   4.1 Form L_P

    call dzeros(L_p_new,3,3)

    do isystem = 1, nsystems

       do ii = 1,3
          do jj = 1,3
             L_p_new(ii,jj) = L_p_new(ii,jj) + gamma_dot(isystem)*schmidt(isystem,ii,jj)
          enddo
       enddo

    enddo

!   5: CALCULATION OF RESIDUAL: RES=DFGRD_ELAS_PRED-DFGRD_ELAS

    call matinv3(DFGRD_ELAS_PRED0_inv,DFGRD_ELAS_PRED0,iflag)

    do ii = 1, 3
       do jj = 1, 3
          Lterm = 0d0
          do pp = 1, 3
             Lterm = Lterm + dfgrd_elas_pred0_inv(ii,pp)*dfgrd_elas_pred(pp,jj)
          enddo
          L_p(ii,jj) = (-Lterm + kroneker(ii,jj))/dtime
       enddo
    enddo


    BIGRES(1) = L_p(1,1) - L_p_new(1,1)
    BIGRES(2) = L_p(2,2) - L_p_new(2,2)
    BIGRES(3) = L_p(3,3) - L_p_new(3,3) 
    BIGRES(4) = L_p(1,2) - L_p_new(1,2)
    BIGRES(5) = L_p(1,3) - L_p_new(1,3)
    BIGRES(6) = L_p(2,3) - L_p_new(2,3)
    BIGRES(7) = L_p(2,1) - L_p_new(2,1)
    BIGRES(8) = L_p(3,1) - L_p_new(3,1)
    BIGRES(9) = L_p(3,2) - L_p_new(3,2)
    
    do isystem = 1, nsystems
       BIGRES(9+isystem) = dgamma(isystem)-dgamma_new(isystem)
    enddo

    CALL norm(dnorm_NR,BIGRES,9+nsystems,1)
    
!    write(*,*)'kinc,iter,dnorm',kinc,iter,dnorm_NR
    

    if(dnorm_NR < TOLER) THEN
       dnorm = 0d0
       goto 201
    endif
 
    if(iter >= nincmax) THEN
       write(*,*)'ERROR, too many iterations'
       PNEWDT = 0.75
       return
    endif

    if(dnorm_NR > 1d10 .and. toolarge <= 1) then
       toolarge = 2
       iter = 0
!       write(*,*)'toolarge RES',noel,npt,kinc

       dnorm = 1d50
       goto 201        
    endif


!   6: NON LINEAR SYSTEM RESOLUTION
!   6.1 Form Jacobian  analytic J=\partial(RESIDUAL \partial(DFGRD_ELAS_PRED) 
        
    CALL  dzeros_tensor4(Mjacob,3,3,3,3)
    CALL  dzeros_tensor4(partial_sigma,3,3,3,3)
          
    do ii = 1, 3
       do jj = 1, 3
          do kk = 1, 3
             do ll = 1, 3
                term = DFGRD_ELAS_PRED0_inv(ii,kk)*kroneker(jj,ll)
                Mjacob(ii,jj,kk,ll) = -(1/dtime)*term
             enddo
          enddo
       enddo
    enddo

    do isystem = 1, nsystems
       do ii = 1, 3
          do jj = 1, 3
             partial_tau(isystem,ii,jj) = 0d0                                  
          enddo
       enddo
    enddo

!   partial Epsilon/Fe (rs,ss,nm,nn)
    do ii = 1, 3
       do jj = 1, 3
          do nm = 1, 3
             do nn = 1, 3
                aux3333(ii,jj,nm,nn) =  0.5d0*(kroneker(ii,nn)*DFGRD_ELAS_PRED(nm,jj)  &
                                     +  kroneker(jj,nn)*DFGRD_ELAS_PRED(nm,ii))
             enddo
          enddo
       enddo
    enddo  

!   partial sigmaij/Fmn
    
    do ii = 1, 3
       do jj = 1, 3
          do nm = 1, 3
             do nn = 1, 3
                do kk = 1, 3
                   do ll = 1, 3
                      partial_sigma(ii,jj,nm,nn) = partial_sigma(ii,jj,nm,nn) &
                                                 + stiff_orient(ii,jj,kk,ll)  &
                                                 * aux3333(kk,ll,nm,nn)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    do isystem = 1,nsystems

       dfactor(isystem) = (1/mm)*gamma_0*                               &
           (abs(tau_resolved(isystem)/tau_pred(isystem)))**((1/mm)-1)   &
           *(1d0/tau_pred(isystem))

       if(tau_resolved(isystem) >= 0d0) then
          signogT = 1.0           
       else
          signogT = -1.0
       endif

       dfactor2(isystem) = dfactor(isystem)*                            &    
           abs(tau_resolved(isystem)/tau_pred(isystem))*signoT
              
       do nm = 1, 3
          do nn =1, 3
             do ii = 1,3
                do jj = 1,3
                   partial_tau(isystem,nm,nn) = partial_tau(isystem,nm,nn)  &
                       + dfactor(isystem)*partial_sigma(ii,jj,nm,nn)*       &
                        schmidt(isystem,ii,jj)
                enddo
             enddo
          enddo
       enddo
       
    enddo

    CALL  dzeros_tensor4(jacob_NR,3,3,3,3)
  
    do ii = 1, 3
       do jj = 1, 3
          do nm = 1, 3
             do nn = 1, 3
                term = 0d0
                do isystem = 1, nsystems ! sum on isystem
                   term = term + partial_tau(isystem,nm,nn)*schmidt(isystem,ii,jj)
                enddo         ! end isystem
                jacob_NR(ii,jj,nm,nn) = Mjacob(ii,jj,nm,nn) - term
             enddo
          enddo
       enddo
    enddo                    
  
    call dzeros(BIGJAC,9+max_nsystems,9+max_nsystems)
    call dzeros(BIGCORR,9+max_nsystems,1)

!   BOX11: partial Lp / partial Fe

    CALL TENS3333(jacob_NR,R1_FE)
    do ii = 1, 9
       do jj = 1, 9
          BIGJAC(ii,jj) = R1_FE(ii, jj)
       enddo
    enddo

!   BOX12: Partial Lp/ partial dgamma

    do jsystem = 1, nsystems
      call dzeros(aux33, 3, 3)
      do isystem = 1, nsystems            
        if(dgamma(jsystem) >= 0) then
           signog = 1.0           
        else
           signog = -1.0
        endif

        aux_escalar = dfactor2(isystem)*q(isystem,jsystem)*    &
             h(gamma_tot_pred,tau0(isystem),taus(isystem),     &
             h0(isystem),h1(isystem))*signog
        do ii = 1, 3
           do jj = 1, 3
              aux33(ii,jj) = aux33(ii,jj) + aux_escalar*schmidt(isystem,ii,jj)
           enddo
        enddo
      enddo

      BIGJAC(1, 9 + jsystem) = aux33(1,1)
      BIGJAC(2, 9 + jsystem) = aux33(2,2)
      BIGJAC(3, 9 + jsystem) = aux33(3,3)
      BIGJAC(4, 9 + jsystem) = aux33(1,2)
      BIGJAC(5, 9 + jsystem) = aux33(1,3)
      BIGJAC(6, 9 + jsystem) = aux33(2,3)
      BIGJAC(7, 9 + jsystem) = aux33(2,1)
      BIGJAC(8, 9 + jsystem) = aux33(3,1)
      BIGJAC(9, 9 + jsystem) = aux33(3,2)                      
    enddo

!   BOX21: Partial dgamma / partial Fe
    do isystem = 1, nsystems
      BIGJAC(9 + isystem, 1) = -dtime*partial_tau(isystem,1,1)
      BIGJAC(9 + isystem, 2) = -dtime*partial_tau(isystem,2,2)
      BIGJAC(9 + isystem, 3) = -dtime*partial_tau(isystem,3,3)
      BIGJAC(9 + isystem, 4) = -dtime*partial_tau(isystem,1,2)
      BIGJAC(9 + isystem, 5) = -dtime*partial_tau(isystem,1,3)
      BIGJAC(9 + isystem, 6) = -dtime*partial_tau(isystem,2,3)
      BIGJAC(9 + isystem, 7) = -dtime*partial_tau(isystem,2,1)
      BIGJAC(9 + isystem, 8) = -dtime*partial_tau(isystem,3,1)
      BIGJAC(9 + isystem, 9) = -dtime*partial_tau(isystem,3,2)
    enddo

!   BOX22: Partial dgamma/dgamma
    do isystem = 1,nsystems
      do jsystem = 1,nsystems
        if(dgamma(jsystem) >= 0d0) then
           signog = 1.0           
        else
           signog = -1.0
        endif

        BIGJAC(9+isystem, 9+jsystem) = kroneker(isystem,jsystem)                  &
                                     + dfactor2(isystem)                          & 
                                     * h(gamma_tot_pred,tau0(isystem),            &
                                         taus(isystem),h0(isystem),h1(isystem))   &
                                     * signog*q(isystem,jsystem)*dtime
      enddo
    enddo
    
      
    ntot = max_nsystems + 9
    nloc = 9 + nsystems
    call gauss_3(BIGJAC,BIGRES,BIGCORR,ntot,nloc,deter,iflag)
    if(iflag == 1) then
       write(*,*) 'error in system of eqs'
       PNEWDT = 0.5
       return
    endif
   
  
 1  FORMAT(21(F10.2,' '))
    
    DFGRD_ELAS_PRED(1,1) = DFGRD_ELAS_PRED(1,1) - BIGCORR(1)
    DFGRD_ELAS_PRED(2,2) = DFGRD_ELAS_PRED(2,2) - BIGCORR(2)
    DFGRD_ELAS_PRED(3,3) = DFGRD_ELAS_PRED(3,3) - BIGCORR(3)
    DFGRD_ELAS_PRED(1,2) = DFGRD_ELAS_PRED(1,2) - BIGCORR(4)
    DFGRD_ELAS_PRED(1,3) = DFGRD_ELAS_PRED(1,3) - BIGCORR(5)
    DFGRD_ELAS_PRED(2,3) = DFGRD_ELAS_PRED(2,3) - BIGCORR(6)
    DFGRD_ELAS_PRED(2,1) = DFGRD_ELAS_PRED(2,1) - BIGCORR(7)
    DFGRD_ELAS_PRED(3,1) = DFGRD_ELAS_PRED(3,1) - BIGCORR(8)
    DFGRD_ELAS_PRED(3,2) = DFGRD_ELAS_PRED(3,2) - BIGCORR(9)
    
    gamma_tot_pred = 0d0
    gamma_tot_act  = 0d0
    do isystem = 1, nsystems
      dgamma(isystem)     = dgamma(isystem)    - BIGCORR(9+isystem)
      gamma_pred(isystem) = gamma_act(isystem) + abs(dgamma(isystem))
      gamma_tot_pred      = gamma_tot_pred     + gamma_pred(isystem)
      gamma_tot_act       = gamma_tot_act      + gamma_act(isystem)
    enddo

    do isystem = 1, nsystems
      tau_pred(isystem) = tau_act(isystem)
      do jsystem = 1, nsystems
        tau_pred(isystem) = tau_pred(isystem) + q(isystem,jsystem)  &
              *h(gamma_tot_pred,tau0(isystem),taus(isystem),        &
               h0(isystem),h1(isystem))*abs(dgamma(jsystem))
      enddo
    enddo
!   
!   NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
!   
          
    call norm(dnorm,BIGCORR,nsystems+9,1)
          
 2  IF(dnorm > 1d10) THEN

      call polar_decomp(aux33,rr_tot,aux33_2,DFGRD_inc)  
      do ii = 1, 3
        do jj = 1, 3
          DFGRD_ELAS_PRED(ii,jj) = 0d0
          do kk = 1, 3
            DFGRD_ELAS_PRED(ii,jj) = DFGRD_ELAS_PRED(ii,jj) +    &
                  rr_tot(ii,kk)*DFGRD_ELAS_t(kk,jj)
            
          enddo
        enddo
      enddo

      CALL GREEN_LAGRANGE(DFGRD_ELAS_pred,EPSILON_EL)
      CALL PCONTRACT2(skirchoff_rot,STIFF_ORIENT,epsilon_el,3)     
      do isystem = 1, nsystems
         do ii = 1, 3
            do jj = 1, 3
               AUX33_2(ii,jj) = schmidt(isystem,ii,jj)
            enddo
         enddo
         CALL PCONTRACT(tau_resolved(isystem),skirchoff_rot,AUX33_2,3)           
      enddo

      gamma_tot_pred = 0d0
      gamma_tot_act  = 0d0
      do isystem = 1, nsystems
         CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,             &
                  tau_resolved(isystem),tau_pred(isystem),noconv)
         
         if(noconv .EQV. .TRUE.) then
            write(*,*) 'ERROR IN GAMMADOT', isystem
            pnewdt = .75
            return                            
         endif
         
         dgamma(isystem)     = gamma_dot(isystem)*dtime
         gamma_pred(isystem) = gamma_act(isystem) + abs(dgamma(isystem))   
         gamma_tot_pred      = gamma_tot_pred     + gamma_pred(isystem)
         gamma_tot_act       = gamma_tot_pred     + gamma_act(isystem)

      enddo

      do isystem = 1, nsystems
         tau_pred(isystem) = tau_act(isystem)
                       
         do jsystem = 1,nsystems
            tau_pred(isystem) = tau_pred(isystem)+ q(isystem,jsystem)    &
                 *h(gamma_tot_pred,tau0(isystem),taus(isystem),          &
                 h0(isystem),h1(isystem))*abs(dgamma(jsystem))
         enddo
      enddo
    endif               

 100  CONTINUE

!      ROTATE STRESSES TO DEFORMED CONFIGURATION  

      CALL POLAR_DECOMP(UU,RR,CC,DFGRD_ELAS_pred)
      CALL TRANSPOSE(RR_tr,RR,3,3)
      CALL ROTATE_TENS3(skirchoff,skirchoff_ROT,RR_tr)

      do ii = 1, 3
         STRESS(ii) = skirchoff(ii,ii)
      end do
      STRESS(4) = skirchoff(1, 2)
      STRESS(5) = skirchoff(1, 3)
      STRESS(6) = skirchoff(2, 3)

     
!     SAVE INTERNAL VARIABLES
!     AND REDEFINE INTIAL STATE FOR JACOBIAN CALCULATION
 
      index = 0
      do ii = 1, 3
        do jj = 1, 3
          index = index + 1      
          statev(index) = DFGRD_ELAS_pred(ii,jj) 
          dfgrd_elas_pred_jacob(ii,jj) = dfgrd_elas_pred(ii,jj)
        enddo
      enddo

      do i = 1, nsystems
         index = index + 1

!     Actualization of gamma_act
         statev(index) = gamma_pred(i)

!     Definition of gammas for Jacobian
         gamma_act_jacob(i) = gamma_act(i)
         gamma_pred_jacob(i)= gamma_pred(i)                   
         dgamma_jacob(i)    = dgamma(i)

      enddo

      do i = 1, nsystems
         index = index + 1

!     Actualization of tau_act
         statev(index) = tau_pred(i)

!     Definition of tau_act_jacob and tau_pred_jacob
         tau_act_jacob(i)  = tau_act(i)       
         tau_pred_jacob(i) = tau_pred(i)

      enddo

      do ii = 1,3
         do jj = 1, 3
            index = index + 1
            statev(index) = L_p(ii,jj)
         enddo
      enddo

!     SAVE THE EULER ANGLES OF EACH IP
!     ROTATED_ORIENT=RR*ORIENT_MAT

      call PMAT(ROTATED_ORIENT_T,RR,ORIENT_MAT,3,3,3)
      call transpose(ROTATED_ORIENT,ROTATED_ORIENT_T,3,3)
      call euler(1,statev(index+1),statev(index+2),statev(index+3),rotated_orient)

!     SAVE acumulated shear strain
      index = index + 4
      statev(index) = gamma_tot_pred
      

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!     CALCULATION OF JACOBIAN
!     numerical derivation
!     From DFGRD0 to DFGRD1+delta_epsilon
!     Prediction based on actual converged solution of NR: DFGRD_ELAS_PRED

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
!     SAME TENSORS AS MAIN NR:
!
!      DFGRD0
!      DFGRD_ELAS_t
!      DFGRD_ELAS_PRED

!     NEW TENSORS
!      dfgrd_inc
!      DFGRD_ELAS_PRED0
     
   
      dtime2 = dtime


!     Saving initial values before each step

      do ii = 1, 3
        do jj = 1, 3
          dfgrd_inc0(ii,jj) = dfgrd_inc(ii,jj)           
        enddo
      enddo
      
!     Initialize to cero the DDSDDE in 3x3x3x3 (aux3333j)
      do ii = 1, 3
        do jj = 1, 3
          do kk = 1, 3
            do ll = 1, 3
              aux3333j(ii,jj,kk,ll) = 0d0
            enddo
          enddo
        enddo
      enddo      

    
!     General loop on the 6 perturbations to obtain Jacobian
!     On each istep a whole NR problem is solved

      do istep = 1, 6
         
!         write(*,*)'istep',istep
         call dzeros(delta_eps_jac,3,3)

         if(istep <= 3) then
            iii = istep
            jjj = istep
         else if(istep == 4) then
            iii = 1
            jjj = 2
         else if(istep == 5) then
            iii = 1
            jjj = 3
         else if(istep == 6) then
            iii = 2
            jjj = 3
         endif

!     Definition of perturbation strain_increment
         
         if(istep <= 3) then
            delta_eps_jac(jjj,iii) =    strain_increment
         else
            delta_eps_jac(iii,jjj) = .5*strain_increment
            delta_eps_jac(jjj,iii) = .5*strain_increment
         endif

!     Obtention of new deformation gradient increment DFGRD_INC for this perturbation

!$$$         do,ii=1,3
!$$$            do,jj=1,3
!$$$               DFGRD_INC(ii,jj)=DFGRD_INC0(ii,jj)+delta_eps_jac(ii,jj)
!$$$c               DFGRD_ELAS_PRED(II,JJ)=DFGRD_ELAS_PRED1(II,JJ)              
!$$$            enddo
!$$$         enddo  
!$$$         write(*,*)'DFGRD_INC_forma1',DFGRD_INC
          do ii = 1, 3
            do jj = 1, 3
               DFGRD_INC(ii,jj) = 0.0              
               do pp = 1, 3
                  DFGRD_INC(ii,jj) = DFGRD_INC(ii,jj) +                      &
                   (kroneker(ii,pp)+delta_eps_jac(ii,pp))*DFGRD_INC0(pp,jj)
               enddo
            enddo
         enddo  
!$$$         write(*,*)'DFGRD_INC_forma2',DFGRD_INC

!     NEW dfgrd_elas_pred0 

         CALL PMAT(DFGRD_ELAS_pred0,DFGRD_INC,DFGRD_ELAS_t,3,3,3)   
           
!     RELOAD INTERNAL VARIABLES
         
         index=0
         
         do,ii=1,3
            do,jj=1,3
               DFGRD_ELAS_PRED(ii,jj)=DFGRD_ELAS_PRED_JACOB(ii,jj)
            enddo
         enddo

         do,i=1,nsystems
            gamma_act(i)=gamma_act_jacob(i)
            gamma_pred(i)=gamma_pred_jacob(i)
            dgamma(i)=dgamma_jacob(i)
            tau_act(i)=tau_act_jacob(i)
            tau_pred(i)=tau_pred_jacob(i)   
         enddo
         
         
!     HERE ALL THE NR LOOP
!     Enter with DFGRD_INC and DFGRD_ELAS_t(ii,jj)
         
!     GLOBAL NEWTON-RAPHSON
         
         dnorm_NR = 1d0
         iter     = 0
         toolarge = 0
!
         DO 101 WHILE(dnorm_NR > TOLER_jac .and. iter <= nincmax_jac)
            
            iter = iter + 1
                    
            
!     1. COMPUTE LAGRANGIAN ELASTIC STRAIN
!     1.1: Decomposition of Fe in U and R, 
                       
           
            CALL green_lagrange(DFGRD_ELAS_pred,EPSILON_EL)
                                 
            CALL PCONTRACT2(skirchoff_rot,STIFF_ORIENT,epsilon_el,3)
            
!     2: PROYECTION OF STRESS IN SYSTEMS
            
            do isystem = 1,nsystems
               do ii = 1, 3
                  do jj = 1, 3
                     AUX33_2(ii,jj) = schmidt(isystem,ii,jj)
                  enddo
               enddo
               CALL PCONTRACT(tau_resolved(isystem),skirchoff_rot, AUX33_2,3)
            enddo
            
     

!     3: CALCULATION OF SLIP RATES
   
            do isystem = 1, nsystems
               
               CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,               &
                   tau_resolved(isystem),tau_pred(isystem),noconv)
               if(noconv .EQV. .TRUE.) then
                 write(*,*)'ERROR IN GAMMADOT',isystem
                 pnewdt = .75
                 return   

               endif
               
               dgamma_new(isystem) = gamma_dot(isystem)*dtime

            enddo
            
            
!     4: CALCULATION OF FP
            
!     4.1 Form L_P
            
            call dzeros(L_p_new,3,3)
            
            do isystem = 1, nsystems
               
               do ii = 1, 3
                  do jj = 1, 3
                     L_p_new(ii,jj) = L_p_new(ii,jj) + gamma_dot(isystem)*schmidt(isystem,ii,jj)
                  enddo
               enddo
               
            enddo
            
!     5: CALCULATION OF RESIDUAL: RES=DFGRD_ELAS_PRED-DFGRD_ELAS

            CALL matinv3(DFGRD_ELAS_PRED0_inv,DFGRD_ELAS_PRED0,iflag) 

          
      
            do ii = 1, 3
               do jj = 1, 3
                  Lterm = 0d0
                  do pp = 1,3
                     Lterm = Lterm + dfgrd_elas_pred0_inv(ii,pp)*dfgrd_elas_pred(pp,jj)
                  enddo
                  L_p(ii,jj) = (-Lterm + kroneker(ii,jj))/dtime
               enddo
            enddo
            
            call dzeros(BIGRES, 9 + max_nsystems, 1)
            
            BIGRES(1) = L_p(1,1) - L_p_new(1,1)
            BIGRES(2) = L_p(2,2) - L_p_new(2,2)
            BIGRES(3) = L_p(3,3) - L_p_new(3,3) 
            BIGRES(4) = L_p(1,2) - L_p_new(1,2)
            BIGRES(5) = L_p(1,3) - L_p_new(1,3)
            BIGRES(6) = L_p(2,3) - L_p_new(2,3)
            BIGRES(7) = L_p(2,1) - L_p_new(2,1)
            BIGRES(8) = L_p(3,1) - L_p_new(3,1)
            BIGRES(9) = L_p(3,2) - L_p_new(3,2)
            
            do isystem = 1,nsystems
               BIGRES(9+isystem) = dgamma(isystem) - dgamma_new(isystem)
            enddo

            CALL norm(dnorm_NR,BIGRES,9+nsystems,1)

!             write(*,*)'JAC: kinc,iter,dnorm',kinc,iter,dnorm_NR
            if(dnorm_NR < TOLER_jac) then
               dnorm = 0d0
               GOTO 202
            endif
            
            if(iter >= nincmax_jac - 1) then          
               write(*,*)'ERROR_JAC, demasiadas iteraciones'
               write(*,*) 'iter,dnorm_NR',iter,dnorm_NR
            endif

            if(dnorm_NR > 1E10 .and. toolarge <= 1) then
               toolarge = 2
               iter     = 0
               dnorm    = 1d50
               goto 202
               
            endif
      


!     6: NON LINEAR SYSTEM RESOLUTION
!     6.1 Form Jacobian  analytic J=\partial(RESIDUAL \partial(DFGRD_ELAS_PRED) 
          
            CALL  dzeros_tensor4(Mjacob,3,3,3,3)
            CALL  dzeros_tensor4(partial_sigma,3,3,3,3)
            
            do ii = 1,3
               do jj = 1,3
                  do kk = 1,3
                     do ll = 1,3
                        term = DFGRD_ELAS_PRED0_inv(ii,kk)*kroneker(jj,ll)
                        Mjacob(ii,jj,kk,ll) = -(1/dtime)*term
                     enddo
                  enddo
               enddo
            enddo
            
            do isystem = 1, nsystems
               do ii = 1, 3
                  do jj = 1, 3
                     partial_tau(isystem,ii,jj) = 0d0                                  
                  enddo
               enddo
            enddo
            
!     partial Epsilon/Fe (rs,ss,nm,nn)
            do ii = 1, 3
               do jj = 1, 3
                  do nm = 1, 3
                     do nn = 1, 3
                        aux3333(ii,jj,nm,nn) = 0.5d0*(kroneker(ii,nn)*DFGRD_ELAS_PRED(nm,jj)  &
                                                    + kroneker(jj,nn)*DFGRD_ELAS_PRED(nm,ii))
                     enddo
                  enddo
               enddo
            enddo  
            
!     partial sigmaij/Fmn
            
            do ii = 1, 3
               do jj = 1, 3
                  do nm = 1, 3
                     do nn = 1, 3
                        do kk = 1, 3
                           do ll = 1, 3
                              partial_sigma(ii,jj,nm,nn) =  partial_sigma(ii,jj,nm,nn) & 
                                                         +  stiff_orient(ii,jj,kk,ll)  &
                                                         *  aux3333(kk,ll,nm,nn)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            
            do isystem = 1, nsystems
               dfactor(isystem) = (1/mm)*gamma_0*                      &
                   (abs(tau_resolved(isystem)/tau_pred(isystem)))      &
                   **((1/mm)-1)*(1d0/tau_pred(isystem))

               do nm = 1, 3
                  do  nn = 1, 3
                     do ii = 1, 3
                        do jj = 1, 3
                           partial_tau(isystem,nm,nn) = partial_tau(isystem,nm,nn)+  &
                                dfactor(isystem)*partial_sigma(ii,jj,nm,nn)*         &
                                schmidt(isystem,ii,jj)
                        enddo
                     enddo
                  enddo
               enddo
               
            enddo
            
            CALL  dzeros_tensor4(jacob_NR,3,3,3,3)
            
            
            do ii = 1, 3
               do jj = 1, 3
                  do nm = 1, 3
                     do nn = 1, 3
                        term = 0d0
                        do isystem = 1, nsystems ! sum on isystem
!     
                           term = term + partial_tau(isystem,nm,nn)*schmidt(isystem,ii,jj)
                        enddo   ! end isystem
                        
                        jacob_NR(ii,jj,nm,nn) = Mjacob(ii,jj,nm,nn) - term
                     enddo
                  enddo
               enddo
            enddo          
            
            call dzeros(BIGJAC,9+max_nsystems,9+max_nsystems)
            call dzeros(BIGCORR,9+max_nsystems,1)
            
!     BOX11: partial Lp / partial Fe
            
            CALL TENS3333(jacob_NR,R1_FE)
            do ii = 1, 9
               do jj = 1, 9
                  BIGJAC(ii,jj) = R1_FE(ii, jj)
               enddo
            enddo
            
!     BOX12: Partial Lp/ partial dgamma

            do jsystem = 1,nsystems
!     
               call dzeros(aux33,3,3)
               do  isystem = 1,nsystems
                  
                  if(tau_resolved(isystem) >= 0d0) then
                     signogT = 1.0           
                  else
                     signogT = -1.0
                  endif
                  
                  if(dgamma(jsystem) >= 0) then
                     signog = 1.0           
                  else
                     signog = -1.0
                  endif

                  aux_escalar = (1/mm)*gamma_0*                                    &
                       (abs(tau_resolved(isystem)/tau_pred(isystem)))**(1/mm)      &
                       *(1d0/tau_pred(isystem))* q(isystem,jsystem)*               &
                       h(gamma_tot_pred,tau0(isystem),taus(isystem),               &
                       h0(isystem),h1(isystem))*signog*signogT

                  do ii = 1, 3
                     do jj = 1, 3
                        aux33(ii,jj) = aux33(ii,jj) + aux_escalar*schmidt(isystem,ii,jj)
                     enddo
                  enddo
               enddo
               
               BIGJAC(1, 9+jsystem) = aux33(1,1)
               BIGJAC(2, 9+jsystem) = aux33(2,2)
               BIGJAC(3, 9+jsystem) = aux33(3,3)
               BIGJAC(4, 9+jsystem) = aux33(1,2)
               BIGJAC(5, 9+jsystem) = aux33(1,3)
               BIGJAC(6, 9+jsystem) = aux33(2,3)
               BIGJAC(7, 9+jsystem) = aux33(2,1)
               BIGJAC(8, 9+jsystem) = aux33(3,1)
               BIGJAC(9, 9+jsystem) = aux33(3,2)                      
            enddo
            
!     BOX21: Partial dgamma / partial Fe
            do isystem = 1,nsystems
               BIGJAC(9+isystem, 1) = -dtime*partial_tau(isystem,1,1)
               BIGJAC(9+isystem, 2) = -dtime*partial_tau(isystem,2,2)
               BIGJAC(9+isystem, 3) = -dtime*partial_tau(isystem,3,3)
               BIGJAC(9+isystem, 4) = -dtime*partial_tau(isystem,1,2)
               BIGJAC(9+isystem, 5) = -dtime*partial_tau(isystem,1,3)
               BIGJAC(9+isystem, 6) = -dtime*partial_tau(isystem,2,3)
               BIGJAC(9+isystem, 7) = -dtime*partial_tau(isystem,2,1)
               BIGJAC(9+isystem, 8) = -dtime*partial_tau(isystem,3,1)
               BIGJAC(9+isystem, 9) = -dtime*partial_tau(isystem,3,2)
               
            enddo

!     BOX22: Partial dgamma/dgamma
            do isystem = 1,nsystems
               do jsystem = 1,nsystems
                  
                  if(tau_resolved(isystem) >= 0d0) then
                     signogT = 1.0           
                  else
                     signogT = -1.0
                  endif
                  
                  if(dgamma(jsystem) >= 0) then
                     signog = 1.0           
                  else
                     signog = -1.0
                  endif
                  
                  BIGJAC(9+isystem,9+jsystem)=kroneker(isystem,jsystem)     &
                       +dfactor(isystem)*                                   &
                       abs(tau_resolved(isystem)/tau_pred(isystem))*        &
                       h(gamma_tot_pred,tau0(isystem),taus(isystem),        &
                       h0(isystem),h1(isystem))*                            &
                       signog*signogT*q(isystem,jsystem)*dtime
               enddo
            enddo
            
            ntot = max_nsystems + 9
            nloc = 9 + nsystems
            call gauss_3(BIGJAC,BIGRES,BIGCORR,ntot,nloc,deter,iflag)
            if(iflag == 1) then
               write(*,*) 'error in system of eqs'
               PNEWDT = 0.5
               return
            endif

            DFGRD_ELAS_PRED(1,1) = DFGRD_ELAS_PRED(1,1) - BIGCORR(1)
            DFGRD_ELAS_PRED(2,2) = DFGRD_ELAS_PRED(2,2) - BIGCORR(2)
            DFGRD_ELAS_PRED(3,3) = DFGRD_ELAS_PRED(3,3) - BIGCORR(3)
            DFGRD_ELAS_PRED(1,2) = DFGRD_ELAS_PRED(1,2) - BIGCORR(4)
            DFGRD_ELAS_PRED(1,3) = DFGRD_ELAS_PRED(1,3) - BIGCORR(5)
            DFGRD_ELAS_PRED(2,3) = DFGRD_ELAS_PRED(2,3) - BIGCORR(6)
            DFGRD_ELAS_PRED(2,1) = DFGRD_ELAS_PRED(2,1) - BIGCORR(7)
            DFGRD_ELAS_PRED(3,1) = DFGRD_ELAS_PRED(3,1) - BIGCORR(8)
            DFGRD_ELAS_PRED(3,2) = DFGRD_ELAS_PRED(3,2) - BIGCORR(9)
            
            gamma_tot_pred = 0d0
            gamma_tot_act  = 0d0
            do isystem = 1, nsystems
               dgamma(isystem) = dgamma(isystem) - BIGCORR(9+isystem)
               gamma_pred(isystem) = gamma_act(isystem)+ abs(dgamma(isystem))
               gamma_tot_pred      = gamma_tot_pred    + gamma_pred(isystem)
               gamma_tot_act       = gamma_tot_act     + gamma_act(isystem)
            enddo


            do isystem = 1, nsystems
               tau_pred(isystem) = tau_act(isystem)
            
               do jsystem = 1,nsystems
                  tau_pred(isystem) = tau_pred(isystem) +                 &
                       q(isystem,jsystem)                                 &
                       *h(gamma_tot_pred,tau0(isystem),taus(isystem),     &
                       h0(isystem),h1(isystem))*abs(dgamma(jsystem))
               enddo
               
            enddo

            
!     
!     NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
!     
               
            call norm(dnorm,BIGCORR,max_nsystems+9,1)
            
 202        IF(dnorm > 1d10) then
               write(*,*) 'In jacob in 202'
               call polar_decomp(aux33,rr_tot,aux33_2,DFGRD_inc)  
               do ii = 1, 3
                  do jj = 1, 3
                     DFGRD_ELAS_PRED(ii,jj) = 0d0
                     do kk = 1, 3
                        DFGRD_ELAS_PRED(ii,jj) = DFGRD_ELAS_PRED(ii,jj)  &
                                               + rr_tot(ii,kk)*DFGRD_ELAS_t(kk,jj)
                        
                     enddo
                  enddo
               enddo
               
               do isystem = 1, nsystems            
                  gamma_pred(isystem) = gamma_act(isystem)
                  dgamma(isystem)     = 0.0
               enddo

               do isystem = 1, nsystems            
                  tau_pred(isystem) = tau_act(isystem)            
               enddo
               
            endif
               
               
 101     CONTINUE
         
!     ROTATE STRESSES TO DEFORMED CONFIGURATION 

         CALL POLAR_DECOMP(UU,RR,CC,DFGRD_ELAS_pred)
         CALL TRANSPOSE(RR_tr,RR,3,3)  
         CALL ROTATE_TENS3(skirchoff2,skirchoff_rot,RR_tr)

!     Numerical Jacobian is J=skirchoff2-skirchoff/strain_increment

         do ii = 1,3
            do jj = 1,3
               if(ii == jj) then
                 aux3333j(ii,jj,iii,jjj) = (skirchoff2(ii,jj) - skirchoff(ii,jj))/strain_increment
                 aux3333j(ii,jj,jjj,iii) = aux3333j(ii,jj,iii,jjj)                    
               else
                 aux3333j(ii,jj,iii,jjj) = (skirchoff2(ii,jj) - skirchoff(ii,jj))/strain_increment
                 aux3333j(jj,ii,iii,jjj) = (skirchoff2(ii,jj) - skirchoff(ii,jj))/strain_increment
                 aux3333j(ii,jj,iii,jjj) = (skirchoff2(ii,jj) - skirchoff(ii,jj))/strain_increment
                 aux3333j(ii,jj,jjj,iii) = (skirchoff2(ii,jj) - skirchoff(ii,jj))/strain_increment
               endif
               
            enddo
        enddo
         
      enddo                     ! END OF 6 perturbation steps
      

      
      CALL TENS3333_sym(aux3333j,ddsdde)

!$$$      write(*,*) 'DDSDDE'
!$$$      do i=1,6
!$$$         write(*,103) (DDSDDE(i,j), j=1,6)
!$$$      enddo
 
      
 110  IF(PNEWDT < 1) then
         CALL TENS3333_sym(STIFF_orient, ddsdde)
      ENDIF

      RETURN 
      END

