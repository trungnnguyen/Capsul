C
c     CAPSUL: CrystAl Plasticity UtiLities tool set 
C     Multiscale Materials Modelling Group, IMDEA Materials
C
C     UMAT CRYSTAL PLASTICITY MODEL 
c     umat_crysplas_CAPSUL_isoh.f
C
c     Short description:
c              -Phenomenlogical hardening (Asaro-Needleman + Voce)
c              -Double residual (Lp + Dgamma)-->Quadratic convergency and robust
c              -Integration of internal variables using several subincrements (now commented)
c              -NR Residual exact, numerical Jacobian 
c              -EXP MAP and derivatives not used and delted. To reincorporate, see older versions

c     Author:  Javier Segurado 

c     Version controls and dates
c     Update 18/01/2016: Split into library with common routines and actual umat
c     Update 22/06/2015: Cleaning, rewritting, comments, from old version 5.3.8
c    
C     Subversion 5.2- 14th June 2012
c                   -- EXP MAP and DEXP MAP are suppressed (linear approach)
c                   -- Includes more general symmetry of stiffness matrix
c                   -- write of messages supressed 
c                   -- Works under parallel in all ABQ versions (6.7-6.11)
c                   -- Now the former file subroutines_v4.f is included in file
C                   -- orientate_tensor4 updated to use explicit rotation formulae from mathematica
c     Subversion 5.3-  18th June 2012
c                   -- Voce Hardening law included
C     Subversion 5.3.1-18th July 2012
c                   -- Improvement of NR jacobian
C     Subversion 5.3.2-June 2013
c                   -- Pseudo-line-search to improove numerical efficiency
C     Version 5.3.8. January 2015
c                   -- Implicit in Fe and tau by a New residual formulation on Lp & Delta_gamma_i
c                   -- Exact Jacobian and Quadratic convergency for new residual. Robust
C     
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c

       INCLUDE 'ABA_PARAM.INC' 
 


      DIMENSiON STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C     List of internal variables:

C     To set the maximum of slip systems    
      integer max_nsystems,nsets,nsystems
      PARAMETER (max_nsystems=30)

c     File management
      CHARACTER*2  FLABEL
      CHARACTER*80 outdir,filename
      INTEGER lenoutdir,nfiles,ns,UR0

c     MATERIAL PARAMETERS

c     Elastic constants
      REAL*8 c11,c12,c44,c13,c33,c66
c     Viscoplastic parameters
      REAL*8 mm,gamma_0
c     Hardening law
      real*8 tau0(max_nsystems),taus(max_nsystems)
      real*8 h0(max_nsystems),h
c     Voce hardening law
      real*8 h1(max_nsystems)
c     Latent hardening matrix
      real*8 q(max_nsystems,max_nsystems)
c     Normal planes and tangential planes
      real*8 m(max_nsystems,3),s(max_nsystems,3)
      real*8 dnorm_m,dnorm_s
c     Schmidt matrices
      real*8 schmidt(max_nsystems,3,3)

C     Stiffness  matrices
      REAL*8 STIFF(3,3,3,3),stiff_orient(3,3,3,3)
              
c     Kinematic tensors: deformation gradients and strains 
      REAL*8 DFGRD_INC(3,3)
      REAL*8 DFGRD0_inv(3,3)
      REAL*8 DFGRD_elas(3,3)
      REAL*8 DFGRD_elas_t(3,3)
      REAL*8 DFGRD_elas_pred(3,3),DFGRD_elas_pred0(3,3)
      REAL*8 DFGRD_elas_pred_inv(3,3),DFGRD_elas_pred0_inv(3,3)
      REAL*8 EPSILON_el(3,3)      
      real*8 RR_tot(3,3),RR_tot_tr(3,3)
      real*8 RR(3,3),RR_tr(3,3),UU(3,3),VV(3,3),CC(3,3)
     
c     Stresses: skirchoff_rot: 2nd Piola-Kirchoff (intermediate configuration)
c               skirchoff: Kirchoff stress (final configuration), Cauchy 
      real*8 skirchoff(3,3)
      REAL*8 skirchoff_rot(3,3)
      
c     Velocity gradients:
      REAL*8 L_p(3,3),L_p_new(3,3),EXP_L_p(3,3),Lterm,term

c     Shear on slip systems:      
      real*8 gamma_pred(max_nsystems)
      real*8 gamma_dot(max_nsystems),gamma_act(max_nsystems)
      real*8 gamma_dot_old(max_nsystems),gama_inc(max_nsystems)
      real*8 gamma_dot_res(max_nsystems),error_hard
      real*8 gamma_TOT,gamma_tot_pred,gamma_tot_act

c     Tau on slip systems:
      real*8 tau_act(max_nsystems),tau_pred(max_nsystems)
      real*8 tau_pred0(max_nsystems)
      real*8 tau_resolved(max_nsystems)
     
c     Definition of orientation 
      REAL*8 orient1(3),orient2(3),orient3(3)
      REAL*8 orient_MAT(3,3),orient_MAT_tr(3,3)
      real*8 ROTATED_ORIENT(3,3),ROTATED_ORIENT_tr(3,3)
         
c     Auxiliar (tmp) integers, scalars, vectors, tensors
      real*8 aux_escalar
      real*8 aux3_1(3),aux3_2(3)
      real*8 aux33(3,3),aux33_2(3,3)
      real*8 aux3333(3,3,3,3)


      integer i,j,ii,jj,kk,ll,pp,qq,nm,nn,iii,jjj
      integer ia,ib,ic,id,ii1,jj1
      integer isystem,jsystem
      logical noconv
      real*8 pi

c     NR resolution, Jacobian, etc
      integer index , iter, kroneker,istep
 
      integer iflag

      real*8 JACOB_NR(3,3,3,3),jacob_NR_i(3,3,3,3)
      REAL*8 dfactor(max_nsystems),dfactor2(max_nsystems)
      real*8 mjacob(3,3,3,3),partial_sigma(3,3,3,3)
      real*8 partial_tau(max_nsystems,3,3)
      real*8 partial_tau2(3,3)
      real*8 dnorm,dnorm_NR,dnorm_old,toler,toler_jac
      real*8 dnorm_hard
      integer nincmax,nincmax_jac

      REAL*8 JACOB3333(3,3,3,3)

      real*8 BIGJAC(max_nsystems+9,max_nsystems+9)
      real*8 BIGRES(max_nsystems+9)
      real*8 BIGCORR(max_nsystems+9),R1_FE(9,9)
      integer ntot,nloc
      real*8 dgamma(max_nsystems), dgamma_new(max_nsystems)
      real*8 dgamma_jacob(max_nsystems),signog,signot

      integer toolarge

c     Numerical tangent stiffness matrix (DDSDDE)
      real*8 skirchoff2(3,3),delta_eps_jac(3,3)
      real*8 dtime2,strain_increment
      real*8 gamma_act_jacob(max_nsystems)
      real*8 gamma_pred_jacob(max_nsystems)
      real*8 tau_pred_jacob(max_nsystems),tau_act_jacob(max_nsystems)
      real*8 dfgrd_elas_pred_jacob(3,3)
      real*8 dfgrd_inc0(3,3)
      real*8 tinterpol,detF
        

c     Common Variables saved once at the beginning
      save mm,gamma_0,tau0,taus,h0,h1
      save q,stiff    
      save s,m
      save nsystems
      save pi
      save toler,toler_jac,nincmax,nincmax_jac,strain_increment
      save implicit_hard

c     Variables for paralelization
      integer noelprocess(10),numprocess
      logical init(10)
      logical init1(10)
      save init
      save init1
      save noelprocess


C     STEP 0: READ CRYSTAL PROPERTIES AND COMPUTE AND SAVE SOME VARIABLES
C     Reads 'crystal.prop' and saves some static variables 

c     This STEP has to be run ONLY ONCE AT BEGINNING (t=0) for the whole model (no loop on elements or GP)
C     'init1' and 'init' variable controls if the STEP has been run, if init=.true. the block is jumped
c     To ensure compatibility for mpi and thread paralelization, paralelization subroutine is called

      call paralelization(init,init1,noel,numprocess,noelprocess)

      if(.not. init1(numprocess+1).and.noel.eq.noelprocess(numprocess+1))   
     1     THEN
         
c     Some constants
         pi=4d0*datan(1d0)
         
c     OPEN THE FILE WITH THE CRYSTAL DEFINITION
         CALL GETOUTDIR( OUTDIR, LENOUTDIR )
         FILENAME = OUTDIR(1:LENOUTDIR)//'/crystal.prop' ! Opens crystal.prop
         UR0=74
         OPEN(UR0,FILE=FILENAME,STATUS='OLD')
         
         CALL READPROPERTIES(UR0,max_nsystems,c11,c12,c44,c13,c33,c66, ! Subroutine readproperties read crystal.prop
     1        gamma_0,mm,nsets,nsystems,s,m,q,tau0,taus,h0,h1,
     1        toler,toler_jac,nincmax,nincmax_jac,strain_increment,
     1        implicit_hard)
         
         CLOSE(unit=UR0)

c     Generation of a fourth order stiffness matrix STIFF
c     The symmetry of the crystal taken from c33, 
c     if c33=0--> cubic (3 constants) if c33!=0 ortrhotropic (6 constants)

         if(abs(c33).GT.1D-10) THEN
            CALL STIFF6(c11,c12,c44,c13,c33,c66,stiff)
         else
            CALL STIFF4(c11,c12,c44,stiff)
         endif      
         
        call sleep(1)
        init1(numprocess+1)=.true.
        call sleep(1)

      endif


C     STEP 1: READ INTEGRATION POINT PROPERTIES (here only orientation)
C     DONE FOR EACH INTEGRATION POINT AT EVERY TIME (using static variables could reduce operation but will result in large data)

      CALL form_orient_MAT(props,orient_MAT,orient_MAT_tr)
      
c     DEFINITION OF THE SCHMDIT MATRICES 
      
      do,ii=1,nsystems
         do,i=1,3
            aux3_1(i)=s(ii,i)
            aux3_2(i)=m(ii,i)
         enddo

c     Orientation of schmidt tensors to material orientation, schmidt

         CALL PTENS(aux33_2,aux3_1,aux3_2)
         CALL ROTATE_TENS3(aux33,aux33_2,ORIENT_MAT_tr) ! Schmidt_cart = Orient_MAT  Schmidt_Crys Orient_MAT_T

         do,i=1,3
            do,j=1,3
               schmidt(ii,i,j)=aux33(i,j)
            enddo 
         enddo
    
      enddo         
      
c     Orientation of stiffness matrix to material orientation, STIFF_ORIENT
   
      CALL ORIENTATE_TENSOR4(STIFF_ORIENT,STIFF,orient_MAT_tr) ! STIFF_ORIENT is STIFF in cartesians

c     READ STATE VARIABLES:
    
c     If time=0 initialize variables

      if(time(2).EQ.0d0) THEN
        
         call dzeros(DFGRD_ELAS,3,3)
         call dzeros(DFGRD_ELAS_t,3,3)
         call dzeros(L_p,3,3)
         
         do,ii=1,3
            DFGRD_ELAS(ii,ii)=1d0
            DFGRD_ELAS_t(ii,ii)=1d0            
         enddo
         
         call dzeros(gamma_act,max_nsystems,1)
         call dzeros(tau_act,max_nsystems,1)

         do,ii=1,nsystems
            tau_act(ii)=tau0(ii)
         enddo

      else

         index=0
         do,ii=1,3
            do,jj=1,3
               index=index+1
               DFGRD_ELAS_t(ii,jj)=statev(index)
            enddo
         enddo
         do,i=1,nsystems
            index=index+1
            gamma_act(i)=statev(index)
         enddo
         do,i=1,nsystems
            index=index+1
            tau_act(i)=statev(index)
         enddo
         do,ii=1,3
            do,jj=1,3
               index=index+1
               L_p(ii,jj)=statev(index)
            enddo
         enddo
      endif

c     PREDICTOR:

c     Deformation gradient increment

      CALL MATINV3(DFGRD0_inv,DFGRD0,iflag)

      if(iflag.NE.0) WRITE(*,*)'ERROR INVERTIGN DFGRD0_inv'

      CALL PMAT(DFGRD_INC,DFGRD1,DFGRD0_inv,3,3,3)   

c     Elastic deformation gradient predictor, DGGRD_ELAS_pred0
      
      CALL PMAT(DFGRD_ELAS_pred0,DFGRD_INC,DFGRD_ELAS_t,3,3,3)  
      cALL matinv3(DFGRD_ELAS_PRED0_inv,DFGRD_ELAS_PRED0,iflag)

c     Jacobian initialization using DFGRD_ELAS_PRED0_inv

      CALL form_Mjacob(dtime,DFGRD_ELAS_PRED0_inv,Mjacob)

c     Elastic deformation gradient prediction based on last L_p
         
      do,ii=1,3
         do,jj=1,3
            L_p_new(ii,jj)=L_p(ii,jj)*(-dtime)
         enddo
      enddo
         
c     CALL EXPMAP(EXP_L_P,noconv,L_P)

      do,ii=1,3
         do,jj=1,3
            EXP_L_P(ii,jj)=kroneker(ii,jj)+L_p_new(ii,jj)
         enddo
      enddo
      
      CALL PMAT(DFGRD_ELAS_PRED,dfgrd_elas_pred0,exp_l_p,3,3,3)
            
                 
c     gamma and CRSS initialization

      gamma_TOT_pred=0d0
      call dzeros(dgamma_new,max_nsystems,1)
      call dzeros(dgamma,max_nsystems,1)

      do,isystem=1,nsystems

         tau_pred(isystem)=tau_act(isystem)
         gamma_pred(isystem)=gamma_act(isystem)
         gamma_TOT_pred=gamma_TOT_pred+gamma_pred(isystem)
         gamma_TOT_act=gamma_TOT_act+gamma_act(isystem)

c     dgamma_prediction = 0

         dgamma(isystem)=0.

      enddo
                
c     GLOBAL NEWTON-RAPHSON

      dnorm_NR=1d20
      iter=0    
      toolarge=0
      
      DO 100 WHILE(dnorm_NR.GT.TOLER)
        
         iter=iter+1

C     1.1: Decomposition of Fe in U and R   
c     1.3: Lagrangian deformation

c     Lagrangian deformation
      CALL GREEN_LAGRANGE(DFGRD_ELAS_pred,EPSILON_EL)
     
c     1.4: Calculation of Kirchoff stress

c     skirchoff_rot DEFINED ON crystal axes (undeformed configuration)
c     skirchoff DEFINED ON final configuration

      CALL PCONTRACT2(skirchoff_rot,STIFF_ORIENT,epsilon_el,3)
      
c     2: PROYECTION OF STRESS IN SYSTEMS

c     2.2: Obtention of resolved tau on all systems

      do,isystem=1,nsystems
         do,ii=1,3
            do,jj=1,3
               AUX33_2(ii,jj)=schmidt(isystem,ii,jj)
            enddo
         enddo
         CALL PCONTRACT(tau_resolved(isystem),skirchoff_rot,AUX33_2,3)

      enddo

c     3: CALCULATION OF SLIP RATES


cc NUEVO ALGORITMO PARA tau_pred
      
         do,isystem=1,nsystems
            
            CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,
     1           tau_resolved(isystem),tau_pred(isystem),noconv)
            
            if(noconv.EQ..TRUE.) THEN
C               write(*,*)'ERROR IN GAMMADOT',isystem
               pnewdt=.75
               return             
             
            endif
            
            dgamma_new(isystem)=gamma_dot(isystem)*dtime
                      
         enddo
         

c     4: CALCULATION OF FP

C     4.1 Form L_P

      call dzeros(L_p_new,3,3)

      do,isystem=1,nsystems

         do,ii=1,3
            do,jj=1,3
               L_p_new(ii,jj)=L_p_new(ii,jj)+gamma_dot(isystem)* 
     1              schmidt(isystem,ii,jj)
            enddo
         enddo

      enddo

c     5: CALCULATION OF RESIDUAL: RES=DFGRD_ELAS_PRED-DFGRD_ELAS

     
      do ii=1,3
         do jj=1,3
            Lterm=0d0
            do pp=1,3
               Lterm=Lterm+dfgrd_elas_pred0_inv(ii,pp)*
     1              dfgrd_elas_pred(pp,jj)
            enddo
            L_p(ii,jj)=(-Lterm+kroneker(ii,jj))/dtime
         enddo
      enddo

c     Forms RESIDUAL VECTOR. If double residual is used keep as it is written

      call dzeros(BIGRES,9+max_nsystems,1)
      CALL form_BIGRES(L_p,L_p_new,dgamma,dgamma_new,nsystems,BIGRES)
     
      CALL norm(dnorm_NR,BIGRES,9+nsystems,1)
      
C      write(*,*)'kinc,iter,dnorm',kinc,iter,dnorm_NR
      
      if(dnorm_NR.LT.TOLER) THEN
         dnorm=0d0
         goto 201
      endif
 
      if(iter.GE.nincmax) THEN
C         write(*,*)'ERROR, too many iterations'
         PNEWDT=0.75
         return
      endif

      if(dnorm_NR.GT.1d10.and.toolarge.LE.1) THEN
	 toolarge=2
	 iter=0
C         write(*,*)'toolarge RES',noel,npt,kinc

         dnorm=1d50
         goto 201        
      endif


c     6: NON LINEAR SYSTEM RESOLUTION
c     6.1 Form Jacobian  analytic J=\partial(RESIDUAL \partial(DFGRD_ELAS_PRED) 

c     partial_sigma = \partial skirchoff / \partial F

      CALL form_partial_sigma_F(DFGRD_ELAS_PRED,stiff_orient,
     1     partial_sigma)      

      do,isystem=1,nsystems
         do,ii=1,3
            do,jj=1,3
               partial_tau(isystem,ii,jj)=0d0                                  
            enddo
         enddo
      enddo

      do,isystem=1,nsystems

c     dfactor=\partial gamma_dot \partial tau_resolved

         dfactor(isystem)=(1/mm)*gamma_0*
     1        (abs(tau_resolved(isystem)/tau_pred(isystem)))**((1/mm)-1)
     1        *(1d0/tau_pred(isystem))

c     dfactor2=\partial gamma_dot \partial tau_pred  
            
         if(tau_resolved(isystem).GE.0d0) then
            signogT=1.         
         else     
            signogT=-1    
         endif

         dfactor2(isystem)=dfactor(isystem)*
     1        abs(tau_resolved(isystem)/tau_pred(isystem))*
     1        signoT
                
         do, nm=1,3
            do, nn=1,3
               do,ii=1,3
                  do,jj=1,3
                     partial_tau(isystem,nm,nn)=
     1                    partial_tau(isystem,nm,nn)+
     1                    dfactor(isystem)*partial_sigma(ii,jj,nm,nn)*
     1                    schmidt(isystem,ii,jj)
                  enddo
               enddo
            enddo
         enddo
         
      enddo

      CALL  dzeros_tensor4(jacob_NR,3,3,3,3)
    
      do,ii=1,3
         do,jj=1,3
            do,nm=1,3
               do,nn=1,3

                  term=0d0

                  do, isystem=1, nsystems ! sum on isystem
c     
                     term=term+partial_tau(isystem,nm,nn)*  
     1                    schmidt(isystem,ii,jj)
                  enddo         ! end isystem

                  jacob_NR(ii,jj,nm,nn)= 
     1                 Mjacob(ii,jj,nm,nn)-term
               enddo
            enddo
         enddo
      enddo                    
    
      call dzeros(BIGJAC,9+max_nsystems,9+max_nsystems)
      call dzeros(BIGCORR,9+max_nsystems,1)

c     BOX11: partial Lp / partial Fe

      CALL TENS3333(jacob_NR,R1_FE)
      do,ii=1,9
         do jj=1,9
            BIGJAC(ii,jj)=R1_FE(ii,jj)
         enddo
      enddo

c     BOX12: Partial Lp/ partial dgamma

      do,jsystem=1,nsystems
c         
         call dzeros(aux33,3,3)
         do  isystem=1,nsystems            
           
            if(dgamma(jsystem).GE.0) THEN
               signog=1.           
            else
               signog=-1.
            endif

            aux_escalar=dfactor2(isystem)*
     1           q(isystem,jsystem)*
     1           h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1           h0(isystem),h1(isystem))*signog
            do,ii=1,3
               do,jj=1,3
                  aux33(ii,jj)=aux33(ii,jj)+
     1                 aux_escalar*schmidt(isystem,ii,jj)
               enddo
            enddo
         enddo
         BIGJAC(1,9+jsystem)=aux33(1,1)
         BIGJAC(2,9+jsystem)=aux33(2,2)
         BIGJAC(3,9+jsystem)=aux33(3,3)
         BIGJAC(4,9+jsystem)=aux33(1,2)
         BIGJAC(5,9+jsystem)=aux33(1,3)
         BIGJAC(6,9+jsystem)=aux33(2,3)
         BIGJAC(7,9+jsystem)=aux33(2,1)
         BIGJAC(8,9+jsystem)=aux33(3,1)
         BIGJAC(9,9+jsystem)=aux33(3,2)                      
      enddo

c     BOX21: Partial dgamma / partial Fe

      do isystem=1,nsystems

         BIGJAC(9+isystem,1)=-dtime*partial_tau(isystem,1,1)
         BIGJAC(9+isystem,2)=-dtime*partial_tau(isystem,2,2)
         BIGJAC(9+isystem,3)=-dtime*partial_tau(isystem,3,3)
         BIGJAC(9+isystem,4)=-dtime*partial_tau(isystem,1,2)
         BIGJAC(9+isystem,5)=-dtime*partial_tau(isystem,1,3)
         BIGJAC(9+isystem,6)=-dtime*partial_tau(isystem,2,3)
         BIGJAC(9+isystem,7)=-dtime*partial_tau(isystem,2,1)
         BIGJAC(9+isystem,8)=-dtime*partial_tau(isystem,3,1)
         BIGJAC(9+isystem,9)=-dtime*partial_tau(isystem,3,2)
       
      enddo

c     BOX22: Partial dgamma/dgamma

      do isystem=1,nsystems
         do jsystem=1,nsystems
                       
            if(dgamma(jsystem).GE.0d0) THEN
               signog=1.           
            else
               signog=-1.
            endif

            BIGJAC(9+isystem,9+jsystem)=kroneker(isystem,jsystem)
     1           +dfactor2(isystem)*         
     1           h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1           h0(isystem),h1(isystem))*
     1           signog*q(isystem,jsystem)*dtime
         enddo
      enddo
      
      ntot=max_nsystems+9
      nloc=9+nsystems

      call gauss_3(BIGJAC,BIGRES,BIGCORR,ntot,nloc,iflag)
      if(iflag.EQ.1) THEN
C         write(*,*)'error in system of eqs'
         PNEWDT=.5
         return
      endif
     
  
 104  FORMAT(21(F10.2,' '))
      
      DFGRD_ELAS_PRED(1,1)=DFGRD_ELAS_PRED(1,1)-BIGCORR(1)
      DFGRD_ELAS_PRED(2,2)=DFGRD_ELAS_PRED(2,2)-BIGCORR(2)
      DFGRD_ELAS_PRED(3,3)=DFGRD_ELAS_PRED(3,3)-BIGCORR(3)
      DFGRD_ELAS_PRED(1,2)=DFGRD_ELAS_PRED(1,2)-BIGCORR(4)
      DFGRD_ELAS_PRED(1,3)=DFGRD_ELAS_PRED(1,3)-BIGCORR(5)
      DFGRD_ELAS_PRED(2,3)=DFGRD_ELAS_PRED(2,3)-BIGCORR(6)
      DFGRD_ELAS_PRED(2,1)=DFGRD_ELAS_PRED(2,1)-BIGCORR(7)
      DFGRD_ELAS_PRED(3,1)=DFGRD_ELAS_PRED(3,1)-BIGCORR(8)
      DFGRD_ELAS_PRED(3,2)=DFGRD_ELAS_PRED(3,2)-BIGCORR(9)
      
      gamma_tot_pred=0d0
      gamma_tot_act=0d0
      do,isystem=1,nsystems
         dgamma(isystem)=dgamma(isystem)-BIGCORR(9+isystem)
         gamma_pred(isystem)=gamma_act(isystem)+abs(dgamma(isystem))
         gamma_tot_pred=gamma_tot_pred+gamma_pred(isystem)
         gamma_tot_act=gamma_tot_act+gamma_act(isystem)
      enddo

      do,isystem=1,nsystems
         tau_pred(isystem)=tau_act(isystem)
       
         do,jsystem=1,nsystems
            tau_pred(isystem)=tau_pred(isystem)+
     1           q(isystem,jsystem)
     1           *h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1           h0(isystem),h1(isystem))*abs(dgamma(jsystem))
         enddo
      
      enddo
c     
c$$$      do,isystem=1,nsystems
c$$$         tau_pred(isystem)=tau_act(isystem)
c$$$         if(abs(h(gamma_tot_pred,tau0(isystem),taus(isystem),
c$$$     1        h0(isystem),h1(isystem))-
c$$$     1        h(gamma_tot_act,tau0(isystem),taus(isystem),
c$$$     1        h0(isystem),h1(isystem)))/h0(isystem).GT..1) THEN
c$$$            write(*,*)'subincrementation',kinc
c$$$            do,jsystem=1,nsystems
c$$$               do,ii=1,11
c$$$                  gamma_tot=gamma_tot_act+gamma_tot_pred*(ii-1.)/11.
c$$$                  tau_pred(isystem)=tau_pred(isystem)+
c$$$     1                 q(isystem,jsystem)
c$$$     1                 *h(gamma_tot,tau0(isystem),taus(isystem),
c$$$     1                 h0(isystem),h1(isystem))*abs(dgamma(jsystem))*.1
c$$$               enddo
c$$$            enddo        
c$$$         else
c$$$            do,jsystem=1,nsystems
c$$$               tau_pred(isystem)=tau_pred(isystem)+
c$$$     1              q(isystem,jsystem)
c$$$     1              *h(gamma_tot_pred,tau0(isystem),taus(isystem),
c$$$     1              h0(isystem),h1(isystem))*abs(dgamma(jsystem))
c$$$            enddo
c$$$         endif
c$$$      enddo
c     
c     NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
c     
            
      call norm(dnorm,BIGCORR,nsystems+9,1)
            
 201  IF(dnorm.GT.1d10) THEN

C         write(*,*)'201b'

         call polar_decomp(aux33,rr_tot,aux33_2,DFGRD_inc)  
         do,ii=1,3
            do,jj=1,3
               DFGRD_ELAS_PRED(ii,jj)=0d0
               do,kk=1,3
                  DFGRD_ELAS_PRED(ii,jj)=DFGRD_ELAS_PRED(ii,jj)+
     1                 rr_tot(ii,kk)*DFGRD_ELAS_t(kk,jj)
                  
               enddo
            enddo
         enddo

         CALL GREEN_LAGRANGE(DFGRD_ELAS_pred,EPSILON_EL)
         CALL PCONTRACT2(skirchoff_rot,STIFF_ORIENT,epsilon_el,3)     
         do,isystem=1,nsystems
            do,ii=1,3
               do,jj=1,3
                  AUX33_2(ii,jj)=schmidt(isystem,ii,jj)
               enddo
            enddo
            CALL PCONTRACT(tau_resolved(isystem),skirchoff_rot,AUX33_2,3)           
         enddo

         gamma_tot_pred=0d0
         gamma_tot_act=0d0
         do,isystem=1,nsystems
            
            CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,
     1           tau_resolved(isystem),tau_pred(isystem),noconv)
            
            if(noconv.EQ..TRUE.) THEN
C               write(*,*)'ERROR IN GAMMADOT',isystem
               pnewdt=.75
               return                            
            endif
            
            dgamma(isystem)=gamma_dot(isystem)*dtime
            gamma_pred(isystem)=gamma_act(isystem)+abs(dgamma(isystem))   
            gamma_tot_pred=gamma_tot_pred+gamma_pred(isystem)
            gamma_tot_act=gamma_tot_pred+gamma_act(isystem)

         enddo

         do,isystem=1,nsystems
            tau_pred(isystem)=tau_act(isystem)
       
            do,jsystem=1,nsystems
               tau_pred(isystem)=tau_pred(isystem)+
     1              q(isystem,jsystem)
     1              *h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1              h0(isystem),h1(isystem))*abs(dgamma(jsystem))
            enddo
      
         enddo
c     
c$$$         do,isystem=1,nsystems
c$$$            tau_pred(isystem)=tau_act(isystem)
c$$$               
c$$$            if(abs(h(gamma_tot_pred,tau0(isystem),taus(isystem),
c$$$     1           h0(isystem),h1(isystem))-
c$$$     1           h(gamma_tot_act,tau0(isystem),taus(isystem),
c$$$     1           h0(isystem),h1(isystem)))/h0(isystem).GT..1) THEN
c$$$               write(*,*)'subincrementation'
c$$$               do,jsystem=1,nsystems
c$$$                  do,ii=1,11
c$$$                     gamma_tot=gamma_tot_act+gamma_tot_pred*(ii-1.)/11.
c$$$                     tau_pred(isystem)=tau_pred(isystem)+
c$$$     1                    q(isystem,jsystem)
c$$$     1                    *h(gamma_tot,tau0(isystem),taus(isystem),
c$$$     1                    h0(isystem),h1(isystem))*
c$$$     1                    abs(dgamma(jsystem))*.1
c$$$                  enddo
c$$$               enddo        
c$$$            else
c$$$               do,jsystem=1,nsystems
c$$$                  tau_pred(isystem)=tau_pred(isystem)+
c$$$     1                 q(isystem,jsystem)
c$$$     1                 *h(gamma_tot_pred,tau0(isystem),taus(isystem),
c$$$     1                 h0(isystem),h1(isystem))*abs(dgamma(jsystem))
c$$$               enddo
c$$$            endif
c$$$         enddo
      endif               

 100  CONTINUE

c      ROTATE STRESSES TO DEFORMED CONFIGURATION  

      CALL POLAR_DECOMP(UU,RR,CC,DFGRD_ELAS_pred)
      CALL TRANSPOSE(RR_tr,RR,3,3)
      CALL ROTATE_TENS3(skirchoff,skirchoff_ROT,RR_tr)

      do,ii=1,3
         STRESS(ii)=skirchoff(ii,ii)
      enddo
      STRESS(4)=skirchoff(1,2)
      STRESS(5)=skirchoff(1,3)
      STRESS(6)=skirchoff(2,3)

     
C     SAVE INTERNAL VARIABLES
C     AND REDEFINE INTIAL STATE FOR JACOBIAN CALCULATION
 
      index=0
      do,ii=1,3
         do,jj=1,3
            index=index+1      
            statev(index)=DFGRD_ELAS_pred(ii,jj) 
            dfgrd_elas_pred_jacob(ii,jj)=dfgrd_elas_pred(ii,jj)
         enddo
      enddo

      do,i=1,nsystems
         index=index+1

c     Actualization of gamma_act
         statev(index)=gamma_pred(i)

c     Definition of gammas for Jacobian
         gamma_act_jacob(i)=gamma_act(i)
         gamma_pred_jacob(i)=gamma_pred(i)                   
         dgamma_jacob(i)=dgamma(i)

      enddo

      do,i=1,nsystems
         index=index+1

c     Actualization of tau_act
         statev(index)=tau_pred(i)

c     Definition of tau_act_jacob and tau_pred_jacob
         tau_act_jacob(i)=tau_act(i)       
         tau_pred_jacob(i)=tau_pred(i)

      enddo

      do,ii=1,3
         do,jj=1,3
            index=index+1
            statev(index)=L_p(ii,jj)
         enddo
      enddo

C     SAVE THE EULER ANGLES OF EACH IP
C     ROTATED_ORIENT=RR*ORIENT_MAT

      call PMAT(ROTATED_ORIENT_tr,RR,ORIENT_MAT,3,3,3)
      call transpose(ROTATED_ORIENT,ROTATED_ORIENT_tr,3,3)
      call euler(1,statev(index+1),statev(index+2),statev(index+3)
     1     ,rotated_orient)

C     SAVE acumulated shear strain
      index=index+4
      statev(index)=gamma_tot_pred
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     CALCULATION OF JACOBIAN
C     numerical derivation
c     From DFGRD0 to DFGRD1+delta_epsilon
C     Prediction based on actual converged solution of NR: DFGRD_ELAS_PRED

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
c     SAME TENSORS AS MAIN NR:
c
c      DFGRD0
c      DFGRD_ELAS_t
c      DFGRD_ELAS_PRED

c     NEW TENSORS
c      dfgrd_inc
c      DFGRD_ELAS_PRED0
     
   
      dtime2=dtime


C     Saving initial values before each step

      do,ii=1,3
         do,jj=1,3
            dfgrd_inc0(ii,jj)=dfgrd_inc(ii,jj)           
         enddo
      enddo
      
C     Initialize to cero the DDSDDE in 3x3x3x3 (aux3333j)

      call dzeros(DDSDDE,6,6)   
    
c     General loop on the 6 perturbations to obtain Jacobian
C     On each istep a whole NR problem is solved

      do,istep=1,6
         
c         write(*,*)'istep',istep
         call dzeros(delta_eps_jac,3,3)

         if(istep.LE.3) THEN
            iii=istep
            jjj=istep
         else if(istep.EQ.4) THEN
            iii=1
            jjj=2
         else if(istep.EQ.5) THEN
            iii=1
            jjj=3
         else if(istep.EQ.6) THEN
            iii=2
            jjj=3
         endif

c     Definition of perturbation strain_increment
         
         if(istep.LE.3) THEN
            delta_eps_jac(jjj,iii)=  strain_increment
         else
            delta_eps_jac(iii,jjj)= .5*strain_increment
            delta_eps_jac(jjj,iii)= .5*strain_increment
         endif

c     Obtention of new deformation gradient increment DFGRD_INC for this perturbation

c$$$         do,ii=1,3
c$$$            do,jj=1,3
c$$$               DFGRD_INC(ii,jj)=DFGRD_INC0(ii,jj)+delta_eps_jac(ii,jj)
c$$$c               DFGRD_ELAS_PRED(II,JJ)=DFGRD_ELAS_PRED1(II,JJ)              
c$$$            enddo
c$$$         enddo  
c$$$         write(*,*)'DFGRD_INC_forma1',DFGRD_INC

          do,ii=1,3
            do,jj=1,3
           
               DFGRD_INC(ii,jj) = 0.0              
               do,pp=1,3
                  DFGRD_INC(ii,jj)=DFGRD_INC(ii,jj)+
     1              (kroneker(ii,pp)+delta_eps_jac(ii,pp))*
     1                 DFGRD_INC0(pp,jj)
               enddo
            enddo
         enddo  
c$$$         write(*,*)'DFGRD_INC_forma2',DFGRD_INC

c     NEW dfgrd_elas_pred0 

         CALL PMAT(DFGRD_ELAS_pred0,DFGRD_INC,DFGRD_ELAS_t,3,3,3)  
         cALL matinv3(DFGRD_ELAS_PRED0_inv,DFGRD_ELAS_PRED0,iflag) 

c     Jacobian initialization using DFGRD_ELAS_PRED0_inv

         CALL form_Mjacob(dtime,DFGRD_ELAS_PRED0_inv,Mjacob)
           
C     RELOAD INTERNAL VARIABLES
         
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
         
         
c     HERE ALL THE NR LOOP
c     Enter with DFGRD_INC and DFGRD_ELAS_t(ii,jj)
         
c     GLOBAL NEWTON-RAPHSON
         
         dnorm_NR=1d0
         iter=0
         toolarge=0
c
         DO 101 WHILE(dnorm_NR.GT.TOLER_jac.and.iter.LE.nincmax_jac)
            
            iter=iter+1
                    
            
c     1. COMPUTE LAGRANGIAN ELASTIC STRAIN
C     1.1: Decomposition of Fe in U and R, 
                       
           
            CALL green_lagrange(DFGRD_ELAS_pred,EPSILON_EL)
                                 
            CALL PCONTRACT2(skirchoff_rot,STIFF_ORIENT,epsilon_el,3)
            
c     2: PROYECTION OF STRESS IN SYSTEMS
            
            do,isystem=1,nsystems
               do,ii=1,3
                  do,jj=1,3
                     AUX33_2(ii,jj)=schmidt(isystem,ii,jj)
                  enddo
               enddo
               CALL PCONTRACT(tau_resolved(isystem),skirchoff_rot,
     1              AUX33_2,3)
            enddo
            
     

c     3: CALCULATION OF SLIP RATES
   
            do,isystem=1,nsystems
               
               CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,
     1              tau_resolved(isystem),tau_pred(isystem),noconv)
               if(noconv.EQ..TRUE.) THEN
C                 write(*,*)'ERROR IN GAMMADOT',isystem
                 pnewdt=.75
                 return   

               endif
               
               dgamma_new(isystem)=gamma_dot(isystem)*dtime

            enddo
            
            
c     4: CALCULATION OF FP
            
C     4.1 Form L_P
            
            call dzeros(L_p_new,3,3)
            
            do,isystem=1,nsystems
               
               do,ii=1,3
                  do,jj=1,3
                     L_p_new(ii,jj)=L_p_new(ii,jj)+gamma_dot(isystem)* 
     1                    schmidt(isystem,ii,jj)
                  enddo
               enddo
               
            enddo
            
c     5: CALCULATION OF RESIDUAL: RES=DFGRD_ELAS_PRED-DFGRD_ELAS

            
                
            do ii=1,3
               do jj=1,3
                  Lterm=0d0
                  do pp=1,3
                     Lterm=Lterm+dfgrd_elas_pred0_inv(ii,pp)*
     1                    dfgrd_elas_pred(pp,jj)
                  enddo
                  L_p(ii,jj)=(-Lterm+kroneker(ii,jj))/dtime
               enddo
            enddo
            
            call dzeros(BIGRES,9+max_nsystems,1)
            CALL form_BIGRES(L_p,L_p_new,dgamma,dgamma_new,nsystems,
     1           BIGRES)

            CALL norm(dnorm_NR,BIGRES,9+nsystems,1)

c             write(*,*)'JAC: kinc,iter,dnorm',kinc,iter,dnorm_NR
            if(dnorm_NR.LT.TOLER_jac) then
               dnorm=0d0
               GOTO 202
            endif
            
            if(iter.GE.nincmax_jac-1) THEN               
C               write(*,*)'ERROR_JAC, demasiadas iteraciones'
C               write(*,*) 'iter,dnorm_NR',iter,dnorm_NR
            endif

            if(dnorm_NR.GT.1E10.and.toolarge.LE.1) THEN
               toolarge=2
               iter=0
               dnorm=1d50
               goto 202
               
            endif
      


c     6: NON LINEAR SYSTEM RESOLUTION
c     6.1 Form Jacobian  analytic J=\partial(RESIDUAL \partial(DFGRD_ELAS_PRED) 

c     partial_sigma = \partial skirchoff / \partial F    
            
            CALL form_partial_sigma_F(DFGRD_ELAS_PRED,stiff_orient,
     1           partial_sigma)    
                                  
            do,isystem=1,nsystems
               do,ii=1,3
                  do,jj=1,3
                     partial_tau(isystem,ii,jj)=0d0                                  
                  enddo
               enddo
            enddo
    
            do,isystem=1,nsystems
               
               dfactor(isystem)=(1/mm)*gamma_0*
     1              (abs(tau_resolved(isystem)/tau_pred(isystem)))
     1              **((1/mm)-1)*(1d0/tau_pred(isystem))

c     dfactor2=\partial gamma_dot \partial tau_pred  
            
               if(tau_resolved(isystem).GE.0d0) then
                  signogT=1.         
               else     
                  signogT=-1    
               endif

               dfactor2(isystem)=dfactor(isystem)*
     1              abs(tau_resolved(isystem)/tau_pred(isystem))*
     1              signoT
               
               do, nm=1,3
                  do, nn=1,3
                     do,ii=1,3
                        do,jj=1,3
                           partial_tau(isystem,nm,nn)=
     1                          partial_tau(isystem,nm,nn)+
     1                          dfactor(isystem)*
     1                          partial_sigma(ii,jj,nm,nn)*
     1                          schmidt(isystem,ii,jj)
                        enddo
                     enddo
                  enddo
               enddo
               
            enddo
            
            CALL  dzeros_tensor4(jacob_NR,3,3,3,3)
            
            
            do,ii=1,3
               do,jj=1,3
                  do,nm=1,3
                     do,nn=1,3
                        
                        term=0d0
                        
                        do, isystem=1, nsystems ! sum on isystem
c     
                           term=term+partial_tau(isystem,nm,nn)*  
     1                          schmidt(isystem,ii,jj)
                        enddo   ! end isystem
                        
                        jacob_NR(ii,jj,nm,nn)= 
     1                       Mjacob(ii,jj,nm,nn)-term
                     enddo
                  enddo
               enddo
            enddo          
            
            call dzeros(BIGJAC,9+max_nsystems,9+max_nsystems)
            call dzeros(BIGCORR,9+max_nsystems,1)
            
c     BOX11: partial Lp / partial Fe
            
            CALL TENS3333(jacob_NR,R1_FE)
            do,ii=1,9
               do jj=1,9
                  BIGJAC(ii,jj)=R1_FE(ii,jj)
               enddo
            enddo
            
c     BOX12: Partial Lp/ partial dgamma

            do,jsystem=1,nsystems
c     
               call dzeros(aux33,3,3)
               do  isystem=1,nsystems
                  
                  if(dgamma(jsystem).GE.0) THEN
                     signog=1.           
                  else
                     signog=-1.
                  endif
                  
                  aux_escalar=dfactor2(isystem)*
     1                 q(isystem,jsystem)*
     1                 h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1                 h0(isystem),h1(isystem))*signog
              
                  do,ii=1,3
                     do,jj=1,3
                        aux33(ii,jj)=aux33(ii,jj)+
     1                       aux_escalar*schmidt(isystem,ii,jj)
                     enddo
                  enddo
               enddo
               
               BIGJAC(1,9+jsystem)=aux33(1,1)
               BIGJAC(2,9+jsystem)=aux33(2,2)
               BIGJAC(3,9+jsystem)=aux33(3,3)
               BIGJAC(4,9+jsystem)=aux33(1,2)
               BIGJAC(5,9+jsystem)=aux33(1,3)
               BIGJAC(6,9+jsystem)=aux33(2,3)
               BIGJAC(7,9+jsystem)=aux33(2,1)
               BIGJAC(8,9+jsystem)=aux33(3,1)
               BIGJAC(9,9+jsystem)=aux33(3,2)                      
            enddo
            
c     BOX21: Partial dgamma / partial Fe
            do isystem=1,nsystems
               BIGJAC(9+isystem,1)=-dtime*partial_tau(isystem,1,1)
               BIGJAC(9+isystem,2)=-dtime*partial_tau(isystem,2,2)
               BIGJAC(9+isystem,3)=-dtime*partial_tau(isystem,3,3)
               BIGJAC(9+isystem,4)=-dtime*partial_tau(isystem,1,2)
               BIGJAC(9+isystem,5)=-dtime*partial_tau(isystem,1,3)
               BIGJAC(9+isystem,6)=-dtime*partial_tau(isystem,2,3)
               BIGJAC(9+isystem,7)=-dtime*partial_tau(isystem,2,1)
               BIGJAC(9+isystem,8)=-dtime*partial_tau(isystem,3,1)
               BIGJAC(9+isystem,9)=-dtime*partial_tau(isystem,3,2)
               
            enddo

c     BOX22: Partial dgamma/dgamma
            do isystem=1,nsystems
               do jsystem=1,nsystems
                                           
                  if(dgamma(jsystem).GE.0) THEN
                     signog=1.           
                  else
                     signog=-1.
                  endif
                  BIGJAC(9+isystem,9+jsystem)=kroneker(isystem,jsystem)
     1                 +dfactor2(isystem)*         
     1                 h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1                 h0(isystem),h1(isystem))*
     1                 signog*q(isystem,jsystem)*dtime
               enddo
            enddo
            
            ntot=max_nsystems+9
            nloc=9+nsystems
            call gauss_3(BIGJAC,BIGRES,BIGCORR,ntot,nloc,iflag)
            if(iflag.EQ.1) THEN
C               write(*,*)'error in system of eqs'
               PNEWDT=.5
               return
            endif

            DFGRD_ELAS_PRED(1,1)=DFGRD_ELAS_PRED(1,1)-BIGCORR(1)
            DFGRD_ELAS_PRED(2,2)=DFGRD_ELAS_PRED(2,2)-BIGCORR(2)
            DFGRD_ELAS_PRED(3,3)=DFGRD_ELAS_PRED(3,3)-BIGCORR(3)
            DFGRD_ELAS_PRED(1,2)=DFGRD_ELAS_PRED(1,2)-BIGCORR(4)
            DFGRD_ELAS_PRED(1,3)=DFGRD_ELAS_PRED(1,3)-BIGCORR(5)
            DFGRD_ELAS_PRED(2,3)=DFGRD_ELAS_PRED(2,3)-BIGCORR(6)
            DFGRD_ELAS_PRED(2,1)=DFGRD_ELAS_PRED(2,1)-BIGCORR(7)
            DFGRD_ELAS_PRED(3,1)=DFGRD_ELAS_PRED(3,1)-BIGCORR(8)
            DFGRD_ELAS_PRED(3,2)=DFGRD_ELAS_PRED(3,2)-BIGCORR(9)
            
            gamma_tot_pred=0d0
            gamma_tot_act=0d0
            do,isystem=1,nsystems
               dgamma(isystem)=dgamma(isystem)-BIGCORR(9+isystem)
               gamma_pred(isystem)=gamma_act(isystem)+
     1              abs(dgamma(isystem))
               gamma_tot_pred=gamma_tot_pred+gamma_pred(isystem)
               gamma_tot_act=gamma_tot_act+gamma_act(isystem)
            enddo

            do,isystem=1,nsystems
               tau_pred(isystem)=tau_act(isystem)
               
               do,jsystem=1,nsystems
                  tau_pred(isystem)=tau_pred(isystem)+
     1                 q(isystem,jsystem)
     1                 *h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1                 h0(isystem),h1(isystem))*abs(dgamma(jsystem))
               enddo
               
            enddo
c     
c$$$            do,isystem=1,nsystems
c$$$               tau_pred(isystem)=tau_act(isystem)
c$$$               if(abs(h(gamma_tot_pred,tau0(isystem),taus(isystem),
c$$$     1              h0(isystem),h1(isystem))-
c$$$     1              h(gamma_tot_act,tau0(isystem),taus(isystem),
c$$$     1              h0(isystem),h1(isystem)))/h0(isystem).GT..1) THEN
c$$$                  write(*,*)'subincrementation'
c$$$                  do,jsystem=1,nsystems
c$$$                     do,ii=1,11
c$$$                        gamma_tot=
c$$$     1                       gamma_tot_act+gamma_tot_pred*(ii-1.)/11.
c$$$                        tau_pred(isystem)=tau_pred(isystem)+
c$$$     1                       q(isystem,jsystem)
c$$$     1                       *h(gamma_tot,tau0(isystem),taus(isystem),
c$$$     1                       h0(isystem),h1(isystem))*
c$$$     1                       abs(dgamma(jsystem))*.1
c$$$                     enddo
c$$$                  enddo        
c$$$               else
c$$$                  do,jsystem=1,nsystems
c$$$                     tau_pred(isystem)=tau_pred(isystem)+
c$$$     1                    q(isystem,jsystem)
c$$$     1                    *h(gamma_tot_pred,tau0(isystem),taus(isystem),
c$$$     1                    h0(isystem),h1(isystem))*abs(dgamma(jsystem))
c$$$                  enddo
c$$$               endif
c$$$            enddo

            
c     
c     NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
c     
               
            call norm(dnorm,BIGCORR,max_nsystems+9,1)
            
 202        IF(dnorm.GT.1d10) THEN

C               write(*,*)'In jacob in 202'

               call polar_decomp(aux33,rr_tot,aux33_2,DFGRD_inc)  
               do,ii=1,3
                  do,jj=1,3
                     DFGRD_ELAS_PRED(ii,jj)=0d0
                     do,kk=1,3
                        DFGRD_ELAS_PRED(ii,jj)=DFGRD_ELAS_PRED(ii,jj)+
     1                       rr_tot(ii,kk)*DFGRD_ELAS_t(kk,jj)
                        
                     enddo
                  enddo
               enddo
               
               do,isystem=1,nsystems            
                  gamma_pred(isystem)=gamma_act(isystem)
                  dgamma(isystem)=0.
               enddo     
               do,isystem=1,nsystems            
                  tau_pred(isystem)=tau_act(isystem)            
               enddo
               
            endif
               
               
 101     CONTINUE
         
c     ROTATE STRESSES TO DEFORMED CONFIGURATION 

         CALL POLAR_DECOMP(UU,RR,CC,DFGRD_ELAS_pred)
         CALL TRANSPOSE(RR_tr,RR,3,3)  
         CALL ROTATE_TENS3(skirchoff2,skirchoff_rot,RR_tr)

c     Numerical Jacobian is J=skirchoff2-skirchoff/strain_increment
c     Directly computed in Nye notation

         DDSDDE(istep,1)= (skirchoff2(1,1)-skirchoff(1,1))/
     1        strain_increment
         DDSDDE(istep,2)= (skirchoff2(2,2)-skirchoff(2,2))/
     1        strain_increment
         DDSDDE(istep,3)= (skirchoff2(3,3)-skirchoff(3,3))/
     1        strain_increment
         DDSDDE(istep,4)= (skirchoff2(1,2)-skirchoff(1,2))/
     1        strain_increment
         DDSDDE(istep,5)= (skirchoff2(1,3)-skirchoff(1,3))/
     1        strain_increment
         DDSDDE(istep,6)= (skirchoff2(2,3)-skirchoff(2,3))/
     1        strain_increment
         
      enddo                     ! END OF 6 perturbation steps
           
c$$$      write(*,*) 'DDSDDE'
c$$$      do i=1,6
c$$$         write(*,103) (DDSDDE(i,j), j=1,6)
c$$$      enddo
 
      
 110  IF(PNEWDT.LT.1) THEN
         CALL TENS3333_sym(STIFF_orient,ddsdde)
      ENDIF

      RETURN 
      END 

c     The library cointatining general purpose subroutines is included

      include 'library_crysplas_CAPSUL_v0'


c     This is the set of subroutines dependent of the constitutive equations
            
      SUBROUTINE READPROPERTIES(UR0,max_nsystems,c11,c12,c44,c13,c33
     1     ,c66,gamma_0,
     1     mm,nsets,nsystems,s,m,q,tau0,taus,h0,h1,
     1     toler,toler_jac,nincmax,nincmax_jac,strain_increment,
     1     implicit_hard)

c     Reads input file with crystal properties and orientations
c     Generates data on each system
c     Normalizes the vectors
      implicit none
      integer ii,jj,i,j,k,max_nsystems
      real*8 c11,c12,c44,c13,c33,c66
      real*8 gamma_0,mm
      integer nsets,nsystems,ur0
      real*8 q(max_nsystems,max_nsystems)
      real*8 qsys(6,6)
      real*8 h0(*),h1(*),tau0(*),taus(*)
      real*8 h0_set(6),tau0_set(6),taus_set(6),h1_set(6)
      real*8 s(max_nsystems,3),m(max_nsystems,3)
      integer system_set(max_nsystems)
 
      real*8 dnorm_m,dnorm_s
      real*8 toler,toler_jac,strain_increment
      integer nincmax,nincmax_jac,implicit_hard

      character*100 line
c     Reading all the crystal data:

      READ(UR0,*)
      READ(UR0,*)
      READ(UR0,111) line
 111  FORMAT(A100)
      c13=0d0
      c33=0d0
      c66=0d0
      READ(line,*,END=112) c11,c12,c44,c13,c33,c66
 112  write(*,*) 'c11 a c66',c11,c12,c44,c13,c33,c66
      READ(UR0,*)
      READ(UR0,*) gamma_0,mm
  
      READ(UR0,*)
      READ(UR0,*) nsets,nsystems
 
      READ(UR0,*)
      do,ii=1,nsystems
         READ(UR0,*) m(ii,1), m(ii,2), m(ii,3),s(ii,1),s(ii,2),s(ii,3),
     1        system_set(ii)
         dnorm_m=dsqrt(m(ii,1)*m(ii,1)+m(ii,2)*m(ii,2)
     1        +m(ii,3)*m(ii,3))
         dnorm_s=dsqrt(s(ii,1)*s(ii,1)+s(ii,2)*s(ii,2)
     1        +s(ii,3)*s(ii,3))
         do,jj=1,3
            m(ii,jj)=m(ii,jj)/dnorm_m
            s(ii,jj)=s(ii,jj)/dnorm_s
         enddo
c         write(*,*) 'system',ii
c         write(*,*) m(ii,1),m(ii,2),m(ii,3)
c         write(*,*) s(ii,1),s(ii,2),s(ii,3)
      enddo

      READ(UR0,*) 
      do,ii=1,nsets
         READ(UR0,*) (qsys(ii,jj),jj=1,nsets)
      enddo
      
      READ(UR0,*) 
      do,ii=1,nsets
         h1_set(ii)=0d0
         READ(UR0,111) line
         READ(line,*,END=113) tau0_set(ii),taus_set(ii),h0_set(ii)
     1        ,h1_set(ii)
 113     write(*,*) 'tau0,taus,h0,h1', tau0_set(ii),taus_set(ii),
     1        h0_set(ii),h1_set(ii)        
      enddo
      
      READ(UR0,*)
      READ(UR0,*) toler,toler_jac,nincmax,nincmax_jac,strain_increment,
     1     implicit_hard
c     Creating the latent-hardening matrix for each individual system
      
      call dzeros(q,max_nsystems,max_nsystems)
      call dzeros(tau0,max_nsystems,1)
      call dzeros(taus,max_nsystems,1)
      call dzeros(h0,max_nsystems,1)
      call dzeros(h1,max_nsystems,1)

      do,ii=1,nsystems
         do,jj=1,nsystems
            q(ii,jj)=qsys(system_set(ii),system_set(jj))
            if(ii.EQ.jj) q(ii,jj)=1d0
         enddo
c         write(*,101) (q(ii,jj),jj=1,nsystems)
      enddo
 101  format(20(F2.0,' '))

c     Creating the properties for each individual system
      do,ii=1,nsystems
         tau0(ii)=tau0_set(system_set(ii))
         taus(ii)=taus_set(system_set(ii))
         h0(ii)=h0_set(system_set(ii))
         h1(ii)=h1_set(system_set(ii))
      enddo

      RETURN 
      END

c     viscolaw --> viscoplastic law defining gammadot as function of 
c     resolved shear stress (tauintern), and other parameters, 
c     internal variables or stress components depending on the model

      SUBROUTINE VISCOLAW(gammadot,gammadot0,expo,tauintern,tau_crit,
     1     noconv)

      implicit none
      logical noconv
      real*8 gammadot0,expo
      real*8 tauintern,tau_crit,gammadot
      real*8 logtau
      if(abs(tauintern).GT.0d0) then
         logtau=log(abs(tauintern/tau_crit))      
      else
         logtau=0
      endif
c      write(*,*)'mm=',expo
c      write(*,*)'logtau=',logtau
      if(logtau.GT.expo*150*log(10.)) THEN
c         write(*,*)'LOGTAU>X',logtau
         gammadot=0
         noconv=.TRUE.
      else
         noconv=.FALSE.         
         gammadot=gammadot0*sign(1.,tauintern)*
     1        abs((tauintern/tau_crit))**(1/expo)     
     
      endif
      RETURN
      END

C
      FUNCTION h(gamma,tau0I,tausI,h0I,h1I)
C     Self Hardening law: Either Asaro-Needleman or Voce 
      implicit none
      REAL*8 h,gamma,tau0I,tausI,h0I,secah,h1I,expo
      REAL*8 tausVOCE

      if(h1I.LE.1D-20) then

         h=h0I*( secah (abs(h0I*gamma/(tausI-tau0I))) )**2
         
      else
         
c     In this case, taus is recalculated as taus-tau0

c        write(*,*)'tau0,taus',tau0I,tausI
        tausVOCE=tausI-tau0I
        if(tausVOCE.LT.0d0) WRITE(*,*)'ERROR IN HARDENING',tausVOCE
        expo=exp(-gamma*h0I/tausVOCE)
        h=h1I*(1d0-expo)+(tausVOCE+h1I*gamma)*(h0I/tausVOCE)*expo  
         
      endif

      RETURN 
      END

       SUBROUTINE form_BIGRES(L_p,L_p_new,dgamma,dgamma_new,
     1     nsystems,BIGRES)
      implicit none
      real*8 L_p(3,3),L_p_new(3,3),dgamma(*),dgamma_new(*),BIGRES(*)
      integer isystem,nsystems
      
      BIGRES(1)=L_p(1,1)-L_p_new(1,1)
      BIGRES(2)=L_p(2,2)-L_p_new(2,2)
      BIGRES(3)=L_p(3,3)-L_p_new(3,3) 
      BIGRES(4)=L_p(1,2)-L_p_new(1,2)
      BIGRES(5)=L_p(1,3)-L_p_new(1,3)
      BIGRES(6)=L_p(2,3)-L_p_new(2,3)
      BIGRES(7)=L_p(2,1)-L_p_new(2,1)
      BIGRES(8)=L_p(3,1)-L_p_new(3,1)
      BIGRES(9)=L_p(3,2)-L_p_new(3,2)
      
      do,isystem=1,nsystems
         BIGRES(9+isystem)=(dgamma(isystem)-dgamma_new(isystem))
      enddo

      return
      end

