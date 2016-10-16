c     UMAT CRYSTAL PLASTICITY 
c     LARGE DEFORMATIONS AND ROTATIONS IN TOTAL DEFORMATIONS
c     FROM BOOK OF SOUZA,PERIC,OWEN
c     Javier Segurado 
c     Version 5.3.8, Ene 2015
c     Modified on 7th of May to include in the STATEV the euler angles
C     Subversion 5.2- 14th June 2012(c                   -- EXP MAP and DEXP MAP are suppressed (linear approach)
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
C     Version 5.3.8. January 2014
c                   -- Implicit in Fe and tau by a New residual formulation on Lp & Delta_gamma_i
c                   -- Exact Jacobian and Quadratic convergency for new residual. Robust 
c
c

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c
!       implicit none
c        use linear_algebra
      implicit real*8 (a-h,o-z)
!       INCLUDE 'ABA_PARAM.INC'  
     
      DIMENSiON STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C     To set a maximum of slip systems
      
      integer max_nsystems,nsets,nsystems
      PARAMETER (max_nsystems=30)

c     FILE OPENING
      CHARACTER*2  FLABEL
      CHARACTER*80 outdir,filename
      INTEGER lenoutdir,nfiles,ns,UR0

c     MATERIAL PARAMETERS

c     elastic
      REAL*8 c11,c12,c44,c13,c33,c66
c     viscoplastic
      REAL*8 mm,gamma_0,tautothemm,logtau,taufrac
c     hardening law
      real*8 tau0(max_nsystems),taus(max_nsystems)
      real*8 h0(max_nsystems),h
c     Voce hardening law
      real*8 h1(max_nsystems)
c     latent hardening matrix
      real*8 q(max_nsystems,max_nsystems)
c     normal planes and tangential planes
      real*8 m(max_nsystems,3),s(max_nsystems,3)
      real*8 dnorm_m,dnorm_s
c     From inputs the definition of schmidt matrices
      real*8 schmidt(max_nsystems,3,3)
c     tau_max
      real*8 tau_max

     
      REAL*8 dfactor(max_nsystems),dfactor2(max_nsystems)
      REAL*8 STIFF(3,3,3,3),stiff_orient(3,3,3,3)
      REAL*8 STIFF_ROTATED(3,3,3,3)

 
      REAL*8 DFGRD_INC(3,3)
      REAL*8 DFGRD0_inv(3,3)

      REAL*8 DFGRD_elas(3,3)
      REAL*8 DFGRD_elas_t(3,3)
      REAL*8 DFGRD_elas_pred(3,3),DFGRD_elas_pred0(3,3)
      REAL*8 DFGRD_elas_pred_inv(3,3),DFGRD_elas_pred0_inv(3,3)
      REAL*8 EPSILON_el(3,3)
     
      
      real*8 RR_tot(3,3),RR_tot_tr(3,3)
      real*8 RR(3,3),RR_tr(3,3),UU(3,3),VV(3,3),CC(3,3)
     

      real*8 skirchoff(3,3),jota
      REAL*8 skirchoff_rot(3,3)
      
      REAL*8 L_p(3,3),L_p_new(3,3),EXP_L_p(3,3),Lterm,term
     
      real*8 gamma_pred(max_nsystems)
      real*8 gamma_dot(max_nsystems),gamma_act(max_nsystems)
      real*8 gamma_dot_old(max_nsystems),gama_inc(max_nsystems)
      real*8 gamma_dot_res(max_nsystems),error_hard
      real*8 gamma_TOT,gamma_tot_pred,gamma_tot_act
      real*8 tau_act(max_nsystems),tau_pred(max_nsystems)
      real*8 tau_pred0(max_nsystems),tau_predn(max_nsystems)
      real*8 tau_resolved(max_nsystems)
     

      REAL*8 orient1(3),orient2(3),orient3(3)
      REAL*8 orient_MAT(3,3),orient_MAT_T(3,3)
      real*8 ROTATED_ORIENT(3,3),ROTATED_ORIENT_T(3,3)
      REAL*8 orient_MAT_X(3,3),rotation_axisX
      REAL*8 rotation_X(3,3)
      REAL*8 vfixed(3,1),vred(3,1),vred_orient(3,1)
      
      real*8 aux_escalar
      real*8 aux3_1(3),aux3_2(3)
      real*8 aux33(3,3),aux33_2(3,3),aux33_3(3,3),aux33_4(3,3)
c      real*8 sqrt2,sqrt3

      real*8 tensor4(3,3,3,3)
      real*8 JACOB_NR(3,3,3,3),jacob_NR_i(3,3,3,3)
      real*8 AUX66(6,6),aux99(9,9)
      real*8 aux3333(3,3,3,3),aux3333_2(3,3,3,3)
      real*8 aux3333j(3,3,3,3)
      real*8 dexpx(3,3,3,3)
      real*8 RES(3,3),dnorm,dnorm_NR,dnorm_old,toler,toler_jac
      real*8 dnorm_hard
      integer nincmax,nincmax_jac
      
      

      integer i,j,ii,jj,kk,ll,pp,qq,nm,nn,iii,jjj
      integer ia,ib,ic,id,ii1,jj1
      integer isystem,jsystem
      logical noconv
      integer index , iter, intern,kroneker,istep,mod
      integer implicit_hard,iflag
      integer intv(2)

c     for jacobian
      real*8 EXP_L_P2(3,3)
      REAL*8 DFGRD_ELAS_tinv(3,3)
      REAL*8 PARTIALFEF(3,3,3,3),PARTIALFFE(3,3,3,3)
     
      REAL*8 JACOB3333(3,3,3,3)
      real*8 sum_s_m_dG(3,3,3,3)
      real*8 dnorm_inc1,dnorm_inc0
      real*8 skirchoff2(3,3),delta_eps_jac(3,3)
      real*8 dtime2,strain_increment
      real*8 gamma_act_jacob(max_nsystems)
      real*8 gamma_pred_jacob(max_nsystems)
      real*8 tau_pred_jacob(max_nsystems),tau_act_jacob(max_nsystems)
      real*8 dfgrd_elas_pred_jacob(3,3)
      real*8 dfgrd_inc0(3,3)
      real*8 tinterpol,detF

      real*8 BIGJAC(max_nsystems+9,max_nsystems+9)
      real*8 BIGRES(max_nsystems+9)
      real*8 BIGCORR(max_nsystems+9),R1_FE(9,9),deter
      integer ntot,nloc
      real*8 dgamma(max_nsystems), dgamma_new(max_nsystems)
      real*8 dgamma_jacob(max_nsystems),signog,signogT
      real*8 mjacob(3,3,3,3),partial_sigma(3,3,3,3)
      real*8 partial_tau(max_nsystems,3,3)
      real*8 partial_tau2(3,3)

c      real*8 pi
      
c     Variables for damping NR
      integer toolarge

c     Common Variables saved once at the beginning
      logical init
      save init
      save mm,gamma_0,tau0,taus,h0,h1,tau_max
      save q,stiff    
      save s,m
      save nsystems
c      save pi,sqrt2,sqrt3,tensor4
      save tensor4
      save toler,toler_jac,nincmax,nincmax_jac,strain_increment
      save implicit_hard

c     CALCULATION OF CRYSTAL PROPERTIES, ONLY ONCE AT BEGINNING
C     Reads 'crystal.prop' and saves: 
c     mm ,gamma_0,tau0,taus,h0,h1
c     q,schmidt,stiff
c     orient_mat,stiff_orient 
c     nsystems
      
c     To prevent problems in multithread when reading properties
c     Note that element 1010101 have to exist!
      if(.not. init ) THEN
c        init=.true.
         if(noel.NE.1010101) then
c            call sleep(1)
         endif
      endif

      if(.not. init ) THEN


c     Some constants
c         pi=4d0*datan(1d0)         
c         sqrt2=1/dsqrt(2D0)
c         sqrt3=1/dsqrt(3D0)
C     Unit 4 order tensor

         do,ii=1,3
            do,jj=1,3
               do,kk=1,3
                  do,ll=1,3
                     tensor4(ii,jj,kk,ll)=0d0
                     if(ii.EQ.kk.AND.jj.EQ.ll) tensor4(ii,jj,kk,ll)=1d0
                  enddo
               enddo
            enddo
         enddo
         
c     OPEN THE FILE WITH THE CRYSTAL DEFINITION
         CALL GETOUTDIR( OUTDIR, LENOUTDIR )
         FILENAME = OUTDIR(1:LENOUTDIR)//'/crystal.prop'
c         write(*,*) filename
         UR0=74
         OPEN(UR0,FILE=FILENAME,STATUS='OLD')

         CALL READPROPERTIES(UR0,max_nsystems,c11,c12,c44,c13,c33,c66,
     1        gamma_0,mm,nsets,nsystems,s,m,q,tau0,taus,h0,h1,
     1        toler,toler_jac,nincmax,nincmax_jac,strain_increment,
     1        implicit_hard)
         
         CLOSE(unit=UR0)

c     Generation of second order stiffness matrix STIFF

         if(abs(c33).GT.1D-10) THEN
            CALL STIFF6(c11,c12,c44,c13,c33,c66,stiff)
         else
            CALL STIFF4(c11,c12,c44,stiff)
         endif      
         
         tau_max=0d0
         do,ii=1,nsystems
            tau_max=max(tau_max,tau0(ii))
         enddo
         init=.true.

      endif

c     DEFINITION OF THE ORIENTATION MATRIX 

c      orient1(1)=props(1)
c      orient1(2)=props(2)
c      orient1(3)=props(3)
      orient1 = props(1:3)
      dnorm = dot_product(orient1, orient1)
      orient1 = orient1/dnorm
      
c      call norm(dnorm,orient1,3,1)
c      do,i=1,3
c         orient1(i)=orient1(i)/dnorm
c      enddo

      orient2(1)=props(4)
      orient2(2)=props(5)
      orient2(3)=props(6)
      call norm(dnorm,orient2,3,1)
      do,i=1,3
         orient2(i)=orient2(i)/dnorm
      enddo
      CALL pvect(orient3,orient1,orient2)
      call norm(dnorm,orient3,3,1)
      do,i=1,3
         orient3(i)=orient3(i)/dnorm
      enddo
      do,i=1,3
         orient_MAT(i,1)=orient1(i)
         orient_MAT(i,2)=orient2(i)
         orient_MAT(i,3)=orient3(i)
      enddo
    
      call transpose(orient_MAT_T,orient_mat,3,3)
      
c     DEFINITION OF THE SCHMDIT MATRICES 
      
      do,ii=1,nsystems
         do,i=1,3
            aux3_1(i)=s(ii,i)
            aux3_2(i)=m(ii,i)
         enddo
         CALL PTENS(aux33,aux3_1,aux3_2)
         CALL PMAT(aux33_2,ORIENT_MAT,aux33,3,3,3)
         CALL PMAT(aux33,aux33_2,ORIENT_MAT_T,3,3,3)

         do,i=1,3
            do,j=1,3
               schmidt(ii,i,j)=aux33(i,j)
            enddo 
         enddo
    
      enddo         
      
c     Orientation of stiffness matrix to material orientation, STIFF_ORIENT
   
      CALL ORIENTATE_TENSOR4(STIFF_ORIENT,STIFF,ORIENT_MAT_T)
c     call rot_tensor4(stiff, orient_mat_t, stiff_orient)
c     write(*, *) "PI = ", pi, " sqrt2 =", sqrt2

c     READ STATE VARIABLES:
    
c     If time=0 initialize variables

      if(time(2).EQ.0d0) THEN
        
         call dzeros(DFGRD_ELAS,3,3)
         call dzeros(DFGRD_ELAS_t,3,3)
         call dzeros(l_p,3,3)
         
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
            
            if(noconv .EQV. .TRUE.) THEN
               write(*,*)'ERROR IN GAMMADOT',isystem
               pnewdt=.75
               return             
             
            endif
            
            dgamma_new(isystem)=gamma_dot(isystem)*dtime
                      
         enddo
         

c      write(*,*) 'after Viscolaw'

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

      cALL matinv3(DFGRD_ELAS_PRED0_inv,DFGRD_ELAS_PRED0,iflag)


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

      CALL norm(dnorm_NR,BIGRES,9+nsystems,1)
      
c      write(*,*)'kinc,iter,dnorm',kinc,iter,dnorm_NR
      

      if(dnorm_NR.LT.TOLER) THEN
         dnorm=0d0
         goto 201
      endif
 
      if(iter.GE.nincmax) THEN
         write(*,*)'ERROR, too many iterations'
         PNEWDT=0.75
         return
      endif

      if(dnorm_NR.GT.1d10.and.toolarge.LE.1) THEN
        toolarge=2
        iter=0
c         write(*,*)'toolarge RES',noel,npt,kinc

         dnorm=1d50
         goto 201        
      endif


c     6: NON LINEAR SYSTEM RESOLUTION
c     6.1 Form Jacobian  analytic J=\partial(RESIDUAL \partial(DFGRD_ELAS_PRED) 
          
      CALL  dzeros_tensor4(Mjacob,3,3,3,3)
      CALL  dzeros_tensor4(partial_sigma,3,3,3,3)
            
      do ii=1,3
         do jj=1,3
            do kk=1,3
               do ll=1,3
                  term=DFGRD_ELAS_PRED0_inv(ii,kk)*
     1                 kroneker(jj,ll)
                  Mjacob(ii,jj,kk,ll)=-(1/dtime)*term
               enddo
            enddo
         enddo
      enddo

      do,isystem=1,nsystems
         do,ii=1,3
            do,jj=1,3
               partial_tau(isystem,ii,jj)=0d0                                  
            enddo
         enddo
      enddo

c     partial Epsilon/Fe (rs,ss,nm,nn)
      do,ii=1,3
         do,jj=1,3
            do,nm=1,3
               do,nn=1,3
                  aux3333(ii,jj,nm,nn)=
     1                 .5d0*(kroneker(ii,nn)* 
     1                 DFGRD_ELAS_PRED(nm,jj)+
     1                 kroneker(jj,nn)* 
     1                 DFGRD_ELAS_PRED(nm,ii))
               enddo
            enddo
         enddo
      enddo  

c     partial sigmaij/Fmn
      
      do,ii=1,3
         do,jj=1,3
            do,nm=1,3
               do,nn=1,3
                  do,kk=1,3
                     do,ll=1,3
                        partial_sigma(ii,jj,nm,nn)=
     1                       partial_sigma(ii,jj,nm,nn)+
     1                       stiff_orient(ii,jj,kk,ll)*
     1                       aux3333(kk,ll,nm,nn)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do,isystem=1,nsystems

         dfactor(isystem)=(1/mm)*gamma_0*
     1        (abs(tau_resolved(isystem)/tau_pred(isystem)))**((1/mm)-1)
     1        *(1d0/tau_pred(isystem))

         if(tau_resolved(isystem).GE.0d0) then
            signogT=1.           
         else
            signogT=-1.
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
      call gauss_3(BIGJAC,BIGRES,BIGCORR,ntot,nloc,deter,iflag)
      if(iflag.EQ.1) THEN
         write(*,*)'error in system of eqs'
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
c     NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
c     
            
      call norm(dnorm,BIGCORR,nsystems+9,1)
            
 201  IF(dnorm.GT.1d10) THEN

c         write(*,*)'201b'

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
            call PCONTRACT(tau_resolved(isystem), skirchoff_rot, 
     1                                AUX33_2, 3)
         enddo

         gamma_tot_pred=0d0
         gamma_tot_act=0d0
         do,isystem=1,nsystems
            
            CALL VISCOLAW(gamma_dot(isystem),gamma_0,mm,
     1           tau_resolved(isystem),tau_pred(isystem),noconv)
            
            if(noconv .EQV. .TRUE.) THEN
               write(*,*)'ERROR IN GAMMADOT',isystem
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

      write(*, *) "gammaDot = ", gammaDot
      write(*, *) "tauCrit = ",  tau_pred
      write(*, *) "tauResl = ",  tau_resolved
     
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

      call PMAT(ROTATED_ORIENT_T,RR,ORIENT_MAT,3,3,3)
      call transpose(ROTATED_ORIENT,ROTATED_ORIENT_T,3,3)
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
      do,ii=1,3
         do,jj=1,3
            do,kk=1,3
               do,ll=1,3
                  aux3333j(ii,jj,kk,ll)=0d0
               enddo
            enddo
         enddo
      enddo      

    
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
               if(noconv .EQV. .TRUE.) THEN
                 write(*,*)'ERROR IN GAMMADOT',isystem
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

            CALL matinv3(DFGRD_ELAS_PRED0_inv,DFGRD_ELAS_PRED0,iflag) 

          
      
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
               BIGRES(9+isystem)=dgamma(isystem)-dgamma_new(isystem)
            enddo

            CALL norm(dnorm_NR,BIGRES,9+nsystems,1)

c             write(*,*)'JAC: kinc,iter,dnorm',kinc,iter,dnorm_NR
            if(dnorm_NR.LT.TOLER_jac) then
               dnorm=0d0
               GOTO 202
            endif
            
            if(iter.GE.nincmax_jac-1) THEN               
               write(*,*)'ERROR_JAC, demasiadas iteraciones'
               write(*,*) 'iter,dnorm_NR',iter,dnorm_NR
            endif

            if(dnorm_NR.GT.1E10.and.toolarge.LE.1) THEN
               toolarge=2
               iter=0
               dnorm=1d50
               goto 202
               
            endif
      


c     6: NON LINEAR SYSTEM RESOLUTION
c     6.1 Form Jacobian  analytic J=\partial(RESIDUAL \partial(DFGRD_ELAS_PRED) 
          
            CALL  dzeros_tensor4(Mjacob,3,3,3,3)
            CALL  dzeros_tensor4(partial_sigma,3,3,3,3)
            
            do ii=1,3
               do jj=1,3
                  do kk=1,3
                     do ll=1,3
                        term=DFGRD_ELAS_PRED0_inv(ii,kk)*
     1                       kroneker(jj,ll)
                        Mjacob(ii,jj,kk,ll)=-(1/dtime)*term
                     enddo
                  enddo
               enddo
            enddo
            
            do,isystem=1,nsystems
               do,ii=1,3
                  do,jj=1,3
                     partial_tau(isystem,ii,jj)=0d0                                  
                  enddo
               enddo
            enddo
            
c     partial Epsilon/Fe (rs,ss,nm,nn)
            do,ii=1,3
               do,jj=1,3
                  do,nm=1,3
                     do,nn=1,3
                        aux3333(ii,jj,nm,nn)=
     1                       .5d0*(kroneker(ii,nn)* 
     1                       DFGRD_ELAS_PRED(nm,jj)+
     1                       kroneker(jj,nn)* 
     1                       DFGRD_ELAS_PRED(nm,ii))
                     enddo
                  enddo
               enddo
            enddo  
            
c     partial sigmaij/Fmn
            
            do,ii=1,3
               do,jj=1,3
                  do,nm=1,3
                     do,nn=1,3
                        do,kk=1,3
                           do,ll=1,3
                              partial_sigma(ii,jj,nm,nn)=
     1                             partial_sigma(ii,jj,nm,nn)+
     1                             stiff_orient(ii,jj,kk,ll)*
     1                             aux3333(kk,ll,nm,nn)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            
            do,isystem=1,nsystems
               
               dfactor(isystem)=(1/mm)*gamma_0*
     1              (abs(tau_resolved(isystem)/tau_pred(isystem)))
     1              **((1/mm)-1)*(1d0/tau_pred(isystem))

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
                  
                  if(tau_resolved(isystem).GE.0d0) then
                     signogT=1.           
                  else
                     signogT=-1.
                  endif
                  
                  if(dgamma(jsystem).GE.0) THEN
                     signog=1.           
                  else
                     signog=-1.
                  endif
                  aux_escalar=(1/mm)*gamma_0*
     1                 (abs(tau_resolved(isystem)/tau_pred(isystem)))
     1                 **(1/mm)
     1                 *(1d0/tau_pred(isystem))*
     1                 q(isystem,jsystem)*
     1                 h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1                 h0(isystem),h1(isystem))*signog*signogT
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
                  
                  if(tau_resolved(isystem).GE.0d0) then
                     signogT=1.           
                  else
                     signogT=-1.
                  endif
                  
                  if(dgamma(jsystem).GE.0) THEN
                     signog=1.           
                  else
                     signog=-1.
                  endif
                  
                  BIGJAC(9+isystem,9+jsystem)=kroneker(isystem,jsystem)
     1                 +dfactor(isystem)*
     1                 abs(tau_resolved(isystem)/tau_pred(isystem))*
     1                 h(gamma_tot_pred,tau0(isystem),taus(isystem),
     1                 h0(isystem),h1(isystem))*
     1                 signog*signogT*q(isystem,jsystem)*dtime
               enddo
            enddo
            
            ntot=max_nsystems+9
            nloc=9+nsystems
            call gauss_3(BIGJAC,BIGRES,BIGCORR,ntot,nloc,deter,iflag)
            if(iflag.EQ.1) THEN
               write(*,*)'error in system of eqs'
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
c     NUMERICAL PROBLEMS WITH LARGE CORRECTIONS-->PLASTIC PREDICTOR
c     
               
            call norm(dnorm,BIGCORR,max_nsystems+9,1)
            
 202        IF(dnorm.GT.1d10) THEN
               write(*,*)'In jacob in 202'
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

         do,ii=1,3
            do,jj=1,3

               if(ii.EQ.jj) THEN
                     aux3333j(ii,jj,iii,jjj)=
     1                    (skirchoff2(ii,jj)-skirchoff(ii,jj))/
     1                    strain_increment
                     aux3333j(ii,jj,jjj,iii)=aux3333j(ii,jj,iii,jjj)                    

               else

                 aux3333j(ii,jj,iii,jjj)=
     1                 (skirchoff2(ii,jj)-skirchoff(ii,jj))/
     1                 strain_increment
                  aux3333j(jj,ii,iii,jjj)=
     1                 (skirchoff2(ii,jj)-skirchoff(ii,jj))/
     1                 strain_increment
                  aux3333j(ii,jj,iii,jjj)=
     1                 (skirchoff2(ii,jj)-skirchoff(ii,jj))/
     1                 strain_increment
                  aux3333j(ii,jj,jjj,iii)=
     1                 (skirchoff2(ii,jj)-skirchoff(ii,jj))/
     1                 strain_increment
                  
               endif
               
            enddo
        enddo
         
      enddo                     ! END OF 6 perturbation steps
      

      
      CALL TENS3333_sym(aux3333j,ddsdde)

c$$$      write(*,*) 'DDSDDE'
c$$$      do i=1,6
c$$$         write(*,103) (DDSDDE(i,j), j=1,6)
c$$$      enddo
 
      
 110  IF(PNEWDT.LT.1) THEN
         CALL TENS3333_sym(STIFF_orient,ddsdde)
      ENDIF

      RETURN 
      END

 
      
      SUBROUTINE STIFF6(c11,c12,c44,c13,c33,c66,stiff)
      implicit none
      REAL*8 c11,c12,c44,c13,c33,c66
      REAL*8 stiff(3,3,3,3)
      integer ii,jj,kk,ll
      do,ii=1,3
         do,jj=1,3
            do,kk=1,3
               do,ll=1,3
                  STIFF(ii,jj,kk,ll)=0D0
               enddo
            enddo
         enddo
      enddo
              
      STIFF(1,1,1,1)=c11
      STIFF(2,2,2,2)=c11
      STIFF(3,3,3,3)=c33
c     
      STIFF(1,1,2,2)=c12
      STIFF(2,2,1,1)=c12
      STIFF(1,1,3,3)=c13
      STIFF(2,2,3,3)=c13
      STIFF(3,3,1,1)=c13
      STIFF(3,3,2,2)=c13
c     
      STIFF(1,2,1,2)=c66
      STIFF(2,1,2,1)=c66
      STIFF(1,2,2,1)=c66
      STIFF(2,1,1,2)=c66
c     
      STIFF(1,3,1,3)=c44
      STIFF(3,1,3,1)=c44
      STIFF(1,3,3,1)=c44
      STIFF(3,1,1,3)=c44
c     
      STIFF(2,3,2,3)=c44
      STIFF(3,2,3,2)=c44
      STIFF(2,3,3,2)=c44
      STIFF(3,2,2,3)=c44

      RETURN
      END

      SUBROUTINE STIFF4(c11,c12,c44,stiff)
      implicit none
      REAL*8 c11,c12,c44,stiff(3,3,3,3)
      integer ii,jj,kk,ll
      do,ii=1,3
         do,jj=1,3
            do,kk=1,3
               do,ll=1,3
                  STIFF(ii,jj,kk,ll)=0D0
               enddo
            enddo
         enddo
      enddo
                  
      STIFF(1,1,1,1)=c11
      STIFF(2,2,2,2)=c11
      STIFF(3,3,3,3)=c11
c     
      STIFF(1,1,2,2)=c12
      STIFF(1,1,3,3)=c12
      STIFF(2,2,1,1)=c12
      STIFF(2,2,3,3)=c12
      STIFF(3,3,1,1)=c12
      STIFF(3,3,2,2)=c12
c     
      STIFF(1,2,1,2)=c44
      STIFF(2,1,2,1)=c44
      STIFF(1,2,2,1)=c44
      STIFF(2,1,1,2)=c44
c     
      STIFF(1,3,1,3)=c44
      STIFF(3,1,3,1)=c44
      STIFF(1,3,3,1)=c44
      STIFF(3,1,1,3)=c44
c     
      STIFF(2,3,2,3)=c44
      STIFF(3,2,3,2)=c44
      STIFF(2,3,3,2)=c44
      STIFF(3,2,2,3)=c44
      RETURN
      END

      SUBROUTINE green_lagrange(FE,epsilon_el)
      implicit none
      real*8 FE(3,3),epsilon_el(3,3)
      real*8 FEt(3,3)
      integer ii,jj
      CALL TRANSPOSE(FEt,FE,3,3)
      CALL PMAT(epsilon_el,FEt,FE,3,3,3)
      
      do,ii=1,3
         do,jj=1,3
            epsilon_el(ii,jj)=epsilon_el(ii,jj)*.5d0
         enddo
         epsilon_el(ii,ii)=epsilon_el(ii,ii)-.5D0
      enddo
      return 
      end


      SUBROUTINE lagrangian_def(CC,epsilon_el)
      implicit none
      real*8 cc(3,3),epsilon_el(3,3)
      integer ii,jj
c     Lagrangian deformation

      do,ii=1,3
         do,jj=1,3
            epsilon_el(ii,jj)=cc(ii,jj)*.5d0
         enddo
         epsilon_el(ii,ii)=epsilon_el(ii,ii)-.5D0
      enddo
      RETURN 
      END


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
         gammadot=gammadot0*sign(1.0d0,tauintern)*
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

      FUNCTION secah(x)
      implicit none
      REAL*8 x,secah
       if(abs(x).GT.50)THEN
         secah=0d0
       else
      secah=2d0/(dexp(x)+dexp(-x))
       endif
      RETURN
      END


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


    
c *****************************************************************************
	subroutine euler (iopt,ph,th,tm,a)
c       
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES (in degrees) OF ca REFERRED TO sa.
c     Copied from VPSC code (Lebensohn et at 1993)
c     *****************************************************************************
	real*8 ph,th,tm,a,pi,sph,cph,sth,cth,stm,ctm
	integer iopt
	dimension a(3,3)
	pi=4.*datan(1.d0)
	
	if(iopt.eq.1) then
	   if(abs(a(3,3)).ge.0.999999d0) then
              if(a(3,3).GT.0) th=0d0
              if(a(3,3).LT.0) th=180d0
	      tm=0d0
	      ph=datan2(a(1,2),a(1,1))
	   else
              th=dacos(a(3,3))
	      sth=dsin(th)
	      tm=datan2(a(1,3)/sth,a(2,3)/sth)
	      ph=datan2(a(3,1)/sth,-a(3,2)/sth)
	   endif
	   th=th*180d0/pi
	   ph=ph*180d0/pi
	   tm=tm*180d0/pi
	else if(iopt.eq.2) then
	   sph=dsin(ph*pi/180.)
	   cph=dcos(ph*pi/180.)
	   sth=dsin(th*pi/180.)
	   cth=dcos(th*pi/180.)
	   stm=dsin(tm*pi/180.)
	   ctm=dcos(tm*pi/180.)
	   a(1,1)=ctm*cph-sph*stm*cth
	   a(2,1)=-stm*cph-sph*ctm*cth
	   a(3,1)=sph*sth
	   a(1,2)=ctm*sph+cph*stm*cth
	   a(2,2)=-sph*stm+cph*ctm*cth
	   a(3,2)=-sth*cph
	   a(1,3)=sth*stm
	   a(2,3)=ctm*sth
	   a(3,3)=cth
	endif
	
	return
	end
	




C$$$     SUBROUTINE to orientate a 4 order rank tensor using a rotation matrix
c$$$
c$$$      SUBROUTINE ORIENTATE_TENSOR4(TENSOR_R,TENSOR,ROT)
c$$$      implicit none
c$$$      REAL*8 TENSOR_R(3,3,3,3),TENSOR(3,3,3,3),ROT(3,3)
c$$$      integer ii,jj,kk,ll,mm,nn,pp,qq
c$$$ 
c$$$      
c$$$      do,ii=1,3
c$$$         do,jj=1,3
c$$$            do,kk=1,3
c$$$               do,ll=1,3
c$$$                  TENSOR_R(ii,jj,kk,ll)=0d0
c$$$
c$$$                  do,mm=1,3
c$$$                     do,nn=1,3
c$$$                        do,pp=1,3
c$$$                           do,qq=1,3
c$$$                              TENSOR_R(ii,jj,kk,ll)=
c$$$     1                             TENSOR_R(ii,jj,kk,ll)+
c$$$     2                             ROT(mm,ii)*ROT(nn,jj)*
c$$$     3                             ROT(pp,kk)*ROT(qq,ll)*
c$$$     4                             TENSOR(mm,nn,pp,qq)
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c$$$
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$      RETURN 
c$$$      END
c

c
C       SUBROUTINE to orientate a 4 order rank tensor 
c       14-Jun-2012, Vicente Herrera
c       Formulas using explicit results from Mathematica
c			
         SUBROUTINE ORIENTATE_TENSOR4(TENSOR_R,TENSOR,ROT)
      
C	Componentes Tensor de rigidez de 4º orden e inicialización dl tensor de 4º orden girado

	implicit none
        REAL*8 TENSOR_R(3,3,3,3),TENSOR(3,3,3,3),ROT(3,3)
	REAL*8 a11,a12,a13,a21,a22,a23,a31,a32,a33
        integer ii,jj,kk,ll
	REAL*8     t11
	REAL*8     t14
	REAL*8     t15
	REAL*8     t12
	REAL*8     t16
	REAL*8     t13
	REAL*8     t41
	REAL*8     t44
	REAL*8     t45
	REAL*8     t42
	REAL*8     t46
	REAL*8     t43
	REAL*8     t51
	REAL*8     t54
	REAL*8     t55
	REAL*8     t52
	REAL*8     t56
	REAL*8     t53
	REAL*8     t21
	REAL*8     t24
	REAL*8     t25
	REAL*8     t22
	REAL*8     t26
	REAL*8     t23
	REAL*8     t61
	REAL*8     t64
	REAL*8     t65
	REAL*8     t62
	REAL*8     t66
	REAL*8     t63
	REAL*8     t31
	REAL*8     t34
	REAL*8     t35
	REAL*8     t32
	REAL*8     t36
	REAL*8     t33
	real*8 UN,DOS,CUA
	
	DOS=2d0
	CUA=4d0
	
	do,ii=1,3
	   do,jj=1,3	
	      do,kk=1,3
		 do,ll=1,3
		    TENSOR_R(ii,jj,kk,ll)=0D0
		 enddo
	      enddo
	   enddo
	enddo

	t11=TENSOR(1,1,1,1)
	t14=TENSOR(1,1,1,2)
	t15=TENSOR(1,1,1,3)
	t12=TENSOR(1,1,2,2)
	t16=TENSOR(1,1,2,3)
	t13=TENSOR(1,1,3,3)
	t41=TENSOR(1,2,1,1)
	t44=TENSOR(1,2,1,2)
	t45=TENSOR(1,2,1,3)
	t42=TENSOR(1,2,2,2)
	t46=TENSOR(1,2,2,3)
	t43=TENSOR(1,2,3,3)
	t51=TENSOR(1,3,1,1)
	t54=TENSOR(1,3,1,2)
	t55=TENSOR(1,3,1,3)
	t52=TENSOR(1,3,2,2)
	t56=TENSOR(1,3,2,3)
	t53=TENSOR(1,3,3,3)
	t21=TENSOR(2,2,1,1)
	t24=TENSOR(2,2,1,2)
	t25=TENSOR(2,2,1,3)
	t22=TENSOR(2,2,2,2)
	t26=TENSOR(2,2,2,3)
	t23=TENSOR(2,2,3,3)
	t61=TENSOR(2,3,1,1)
	t64=TENSOR(2,3,1,2)
	t65=TENSOR(2,3,1,3)
	t62=TENSOR(2,3,2,2)
	t66=TENSOR(2,3,2,3)
	t63=TENSOR(2,3,3,3)
	t31=TENSOR(3,3,1,1)
	t34=TENSOR(3,3,1,2)
	t35=TENSOR(3,3,1,3)
	t32=TENSOR(3,3,2,2)
	t36=TENSOR(3,3,2,3)
	t33=TENSOR(3,3,3,3)
      
      

c 	Montaje matriz de giro ROT



	a11=ROT(1,1)
	a21=ROT(1,2)
	a31=ROT(1,3)
	a12=ROT(2,1)
	a22=ROT(2,2)
	a32=ROT(2,3)
	a13=ROT(3,1)
	a23=ROT(3,2)
	a33=ROT(3,3)
	

c      Asignación explicita de las componentes del resultado de Mathematica


             TENSOR_R(1,1,1,1)=a11**4*t11 + a11**2*a12**2*t12 + 
     1  a11**2*a13**2*t13 + DOS*a11**3*a12*t14 + 
     1  DOS*a11**3*a13*t15 + DOS*a11**2*a12*a13*t16 + 
     1  a11**2*a12**2*t21 + a12**4*t22 + 
     1  a12**2*a13**2*t23 + DOS*a11*a12**3*t24 + 
     1  DOS*a11*a12**2*a13*t25 + DOS*a12**3*a13*t26 + 
     1  a11**2*a13**2*t31 + a12**2*a13**2*t32 + 
     1  a13**4*t33 + DOS*a11*a12*a13**2*t34 + 
     1  DOS*a11*a13**3*t35 + DOS*a12*a13**3*t36 + 
     1  DOS*a11**3*a12*t41 + DOS*a11*a12**3*t42 + 
     1  DOS*a11*a12*a13**2*t43 + CUA*a11**2*a12**2*t44 + 
     1  CUA*a11**2*a12*a13*t45 + CUA*a11*a12**2*a13*t46 + 
     1  DOS*a11**3*a13*t51 + DOS*a11*a12**2*a13*t52 + 
     1  DOS*a11*a13**3*t53 + CUA*a11**2*a12*a13*t54 + 
     1  CUA*a11**2*a13**2*t55 + CUA*a11*a12*a13**2*t56 + 
     1  DOS*a11**2*a12*a13*t61 + DOS*a12**3*a13*t62 + 
     1  DOS*a12*a13**3*t63 + CUA*a11*a12**2*a13*t64 + 
     1  CUA*a11*a12*a13**2*t65 + CUA*a12**2*a13**2*t66
       TENSOR_R(1,1,1,2)=a11**3*a21*t11 + a11**2*a12*a22*t12 + 
     1  a11**2*a13*a23*t13 + a11**2*a12*a21*t14 + 
     1  a11**3*a22*t14 + a11**2*a13*a21*t15 + 
     1  a11**3*a23*t15 + a11**2*a13*a22*t16 + 
     1  a11**2*a12*a23*t16 + a11*a12**2*a21*t21 + 
     1  a12**3*a22*t22 + a12**2*a13*a23*t23 + 
     1  a12**3*a21*t24 + a11*a12**2*a22*t24 + 
     1  a12**2*a13*a21*t25 + a11*a12**2*a23*t25 + 
     1  a12**2*a13*a22*t26 + a12**3*a23*t26 + 
     1  a11*a13**2*a21*t31 + a12*a13**2*a22*t32 + 
     1  a13**3*a23*t33 + a12*a13**2*a21*t34 + 
     1  a11*a13**2*a22*t34 + a13**3*a21*t35 + 
     1  a11*a13**2*a23*t35 + a13**3*a22*t36 + 
     1  a12*a13**2*a23*t36 + DOS*a11**2*a12*a21*t41 + 
     1  DOS*a11*a12**2*a22*t42 + DOS*a11*a12*a13*a23*t43 + 
     1  DOS*a11*a12**2*a21*t44 + DOS*a11**2*a12*a22*t44 + 
     1  DOS*a11*a12*a13*a21*t45 + DOS*a11**2*a12*a23*t45 + 
     1  DOS*a11*a12*a13*a22*t46 + DOS*a11*a12**2*a23*t46 + 
     1  DOS*a11**2*a13*a21*t51 + DOS*a11*a12*a13*a22*t52 + 
     1  DOS*a11*a13**2*a23*t53 + DOS*a11*a12*a13*a21*t54 + 
     1  DOS*a11**2*a13*a22*t54 + DOS*a11*a13**2*a21*t55 + 
     1  DOS*a11**2*a13*a23*t55 + DOS*a11*a13**2*a22*t56 + 
     1  DOS*a11*a12*a13*a23*t56 + DOS*a11*a12*a13*a21*t61 + 
     1  DOS*a12**2*a13*a22*t62 + DOS*a12*a13**2*a23*t63 + 
     1  DOS*a12**2*a13*a21*t64 + DOS*a11*a12*a13*a22*t64 + 
     1  DOS*a12*a13**2*a21*t65 + DOS*a11*a12*a13*a23*t65 + 
     1  DOS*a12*a13**2*a22*t66 + DOS*a12**2*a13*a23*t66
       TENSOR_R(1,1,1,3)=a11**3*a31*t11 + a11**2*a12*a32*t12 + 
     1  a11**2*a13*a33*t13 + a11**2*a12*a31*t14 + 
     1  a11**3*a32*t14 + a11**2*a13*a31*t15 + 
     1  a11**3*a33*t15 + a11**2*a13*a32*t16 + 
     1  a11**2*a12*a33*t16 + a11*a12**2*a31*t21 + 
     1  a12**3*a32*t22 + a12**2*a13*a33*t23 + 
     1  a12**3*a31*t24 + a11*a12**2*a32*t24 + 
     1  a12**2*a13*a31*t25 + a11*a12**2*a33*t25 + 
     1  a12**2*a13*a32*t26 + a12**3*a33*t26 + 
     1  a11*a13**2*a31*t31 + a12*a13**2*a32*t32 + 
     1  a13**3*a33*t33 + a12*a13**2*a31*t34 + 
     1  a11*a13**2*a32*t34 + a13**3*a31*t35 + 
     1  a11*a13**2*a33*t35 + a13**3*a32*t36 + 
     1  a12*a13**2*a33*t36 + DOS*a11**2*a12*a31*t41 + 
     1  DOS*a11*a12**2*a32*t42 + DOS*a11*a12*a13*a33*t43 + 
     1  DOS*a11*a12**2*a31*t44 + DOS*a11**2*a12*a32*t44 + 
     1  DOS*a11*a12*a13*a31*t45 + DOS*a11**2*a12*a33*t45 + 
     1  DOS*a11*a12*a13*a32*t46 + DOS*a11*a12**2*a33*t46 + 
     1  DOS*a11**2*a13*a31*t51 + DOS*a11*a12*a13*a32*t52 + 
     1  DOS*a11*a13**2*a33*t53 + DOS*a11*a12*a13*a31*t54 + 
     1  DOS*a11**2*a13*a32*t54 + DOS*a11*a13**2*a31*t55 + 
     1  DOS*a11**2*a13*a33*t55 + DOS*a11*a13**2*a32*t56 + 
     1  DOS*a11*a12*a13*a33*t56 + DOS*a11*a12*a13*a31*t61 + 
     1  DOS*a12**2*a13*a32*t62 + DOS*a12*a13**2*a33*t63 + 
     1  DOS*a12**2*a13*a31*t64 + DOS*a11*a12*a13*a32*t64 + 
     1  DOS*a12*a13**2*a31*t65 + DOS*a11*a12*a13*a33*t65 + 
     1  DOS*a12*a13**2*a32*t66 + DOS*a12**2*a13*a33*t66
       TENSOR_R(1,1,2,2)=a11**2*a21**2*t11 + a11**2*a22**2*t12 + 
     1  a11**2*a23**2*t13 + DOS*a11**2*a21*a22*t14 + 
     1  DOS*a11**2*a21*a23*t15 + DOS*a11**2*a22*a23*t16 + 
     1  a12**2*a21**2*t21 + a12**2*a22**2*t22 + 
     1  a12**2*a23**2*t23 + DOS*a12**2*a21*a22*t24 + 
     1  DOS*a12**2*a21*a23*t25 + DOS*a12**2*a22*a23*t26 + 
     1  a13**2*a21**2*t31 + a13**2*a22**2*t32 + 
     1  a13**2*a23**2*t33 + DOS*a13**2*a21*a22*t34 + 
     1  DOS*a13**2*a21*a23*t35 + DOS*a13**2*a22*a23*t36 + 
     1  DOS*a11*a12*a21**2*t41 + DOS*a11*a12*a22**2*t42 + 
     1  DOS*a11*a12*a23**2*t43 + CUA*a11*a12*a21*a22*t44 + 
     1  CUA*a11*a12*a21*a23*t45 + CUA*a11*a12*a22*a23*t46 + 
     1  DOS*a11*a13*a21**2*t51 + DOS*a11*a13*a22**2*t52 + 
     1  DOS*a11*a13*a23**2*t53 + CUA*a11*a13*a21*a22*t54 + 
     1  CUA*a11*a13*a21*a23*t55 + CUA*a11*a13*a22*a23*t56 + 
     1  DOS*a12*a13*a21**2*t61 + DOS*a12*a13*a22**2*t62 + 
     1  DOS*a12*a13*a23**2*t63 + CUA*a12*a13*a21*a22*t64 + 
     1  CUA*a12*a13*a21*a23*t65 + CUA*a12*a13*a22*a23*t66
       TENSOR_R(1,1,2,3)=a11**2*a21*a31*t11 + a11**2*a22*a32*t12 + 
     1  a11**2*a23*a33*t13 + a11**2*a22*a31*t14 + 
     1  a11**2*a21*a32*t14 + a11**2*a23*a31*t15 + 
     1  a11**2*a21*a33*t15 + a11**2*a23*a32*t16 + 
     1  a11**2*a22*a33*t16 + a12**2*a21*a31*t21 + 
     1  a12**2*a22*a32*t22 + a12**2*a23*a33*t23 + 
     1  a12**2*a22*a31*t24 + a12**2*a21*a32*t24 + 
     1  a12**2*a23*a31*t25 + a12**2*a21*a33*t25 + 
     1  a12**2*a23*a32*t26 + a12**2*a22*a33*t26 + 
     1  a13**2*a21*a31*t31 + a13**2*a22*a32*t32 + 
     1  a13**2*a23*a33*t33 + a13**2*a22*a31*t34 + 
     1  a13**2*a21*a32*t34 + a13**2*a23*a31*t35 + 
     1  a13**2*a21*a33*t35 + a13**2*a23*a32*t36 + 
     1  a13**2*a22*a33*t36 + DOS*a11*a12*a21*a31*t41 + 
     1  DOS*a11*a12*a22*a32*t42 + DOS*a11*a12*a23*a33*t43 + 
     1  DOS*a11*a12*a22*a31*t44 + DOS*a11*a12*a21*a32*t44 + 
     1  DOS*a11*a12*a23*a31*t45 + DOS*a11*a12*a21*a33*t45 + 
     1  DOS*a11*a12*a23*a32*t46 + DOS*a11*a12*a22*a33*t46 + 
     1  DOS*a11*a13*a21*a31*t51 + DOS*a11*a13*a22*a32*t52 + 
     1  DOS*a11*a13*a23*a33*t53 + DOS*a11*a13*a22*a31*t54 + 
     1  DOS*a11*a13*a21*a32*t54 + DOS*a11*a13*a23*a31*t55 + 
     1  DOS*a11*a13*a21*a33*t55 + DOS*a11*a13*a23*a32*t56 + 
     1  DOS*a11*a13*a22*a33*t56 + DOS*a12*a13*a21*a31*t61 + 
     1  DOS*a12*a13*a22*a32*t62 + DOS*a12*a13*a23*a33*t63 + 
     1  DOS*a12*a13*a22*a31*t64 + DOS*a12*a13*a21*a32*t64 + 
     1  DOS*a12*a13*a23*a31*t65 + DOS*a12*a13*a21*a33*t65 + 
     1  DOS*a12*a13*a23*a32*t66 + DOS*a12*a13*a22*a33*t66
       TENSOR_R(1,1,3,3)=a11**2*a31**2*t11 + a11**2*a32**2*t12 + 
     1  a11**2*a33**2*t13 + DOS*a11**2*a31*a32*t14 + 
     1  DOS*a11**2*a31*a33*t15 + DOS*a11**2*a32*a33*t16 + 
     1  a12**2*a31**2*t21 + a12**2*a32**2*t22 + 
     1  a12**2*a33**2*t23 + DOS*a12**2*a31*a32*t24 + 
     1  DOS*a12**2*a31*a33*t25 + DOS*a12**2*a32*a33*t26 + 
     1  a13**2*a31**2*t31 + a13**2*a32**2*t32 + 
     1  a13**2*a33**2*t33 + DOS*a13**2*a31*a32*t34 + 
     1  DOS*a13**2*a31*a33*t35 + DOS*a13**2*a32*a33*t36 + 
     1  DOS*a11*a12*a31**2*t41 + DOS*a11*a12*a32**2*t42 + 
     1  DOS*a11*a12*a33**2*t43 + CUA*a11*a12*a31*a32*t44 + 
     1  CUA*a11*a12*a31*a33*t45 + CUA*a11*a12*a32*a33*t46 + 
     1  DOS*a11*a13*a31**2*t51 + DOS*a11*a13*a32**2*t52 + 
     1  DOS*a11*a13*a33**2*t53 + CUA*a11*a13*a31*a32*t54 + 
     1  CUA*a11*a13*a31*a33*t55 + CUA*a11*a13*a32*a33*t56 + 
     1  DOS*a12*a13*a31**2*t61 + DOS*a12*a13*a32**2*t62 + 
     1  DOS*a12*a13*a33**2*t63 + CUA*a12*a13*a31*a32*t64 + 
     1  CUA*a12*a13*a31*a33*t65 + CUA*a12*a13*a32*a33*t66
       TENSOR_R(1,2,1,1)=a11**3*a21*t11 + a11*a12**2*a21*t12 + 
     1  a11*a13**2*a21*t13 + DOS*a11**2*a12*a21*t14 + 
     1  DOS*a11**2*a13*a21*t15 + DOS*a11*a12*a13*a21*t16 + 
     1  a11**2*a12*a22*t21 + a12**3*a22*t22 + 
     1  a12*a13**2*a22*t23 + DOS*a11*a12**2*a22*t24 + 
     1  DOS*a11*a12*a13*a22*t25 + DOS*a12**2*a13*a22*t26 + 
     1  a11**2*a13*a23*t31 + a12**2*a13*a23*t32 + 
     1  a13**3*a23*t33 + DOS*a11*a12*a13*a23*t34 + 
     1  DOS*a11*a13**2*a23*t35 + DOS*a12*a13**2*a23*t36 + 
     1  a11**2*a12*a21*t41 + a11**3*a22*t41 + 
     1  a12**3*a21*t42 + a11*a12**2*a22*t42 + 
     1  a12*a13**2*a21*t43 + a11*a13**2*a22*t43 + 
     1  DOS*a11*a12**2*a21*t44 + DOS*a11**2*a12*a22*t44 + 
     1  DOS*a11*a12*a13*a21*t45 + DOS*a11**2*a13*a22*t45 + 
     1  DOS*a12**2*a13*a21*t46 + DOS*a11*a12*a13*a22*t46 + 
     1  a11**2*a13*a21*t51 + a11**3*a23*t51 + 
     1  a12**2*a13*a21*t52 + a11*a12**2*a23*t52 + 
     1  a13**3*a21*t53 + a11*a13**2*a23*t53 + 
     1  DOS*a11*a12*a13*a21*t54 + DOS*a11**2*a12*a23*t54 + 
     1  DOS*a11*a13**2*a21*t55 + DOS*a11**2*a13*a23*t55 + 
     1  DOS*a12*a13**2*a21*t56 + DOS*a11*a12*a13*a23*t56 + 
     1  a11**2*a13*a22*t61 + a11**2*a12*a23*t61 + 
     1  a12**2*a13*a22*t62 + a12**3*a23*t62 + 
     1  a13**3*a22*t63 + a12*a13**2*a23*t63 + 
     1  DOS*a11*a12*a13*a22*t64 + DOS*a11*a12**2*a23*t64 + 
     1  DOS*a11*a13**2*a22*t65 + DOS*a11*a12*a13*a23*t65 + 
     1  DOS*a12*a13**2*a22*t66 + DOS*a12**2*a13*a23*t66
       TENSOR_R(1,2,1,2)=a11**2*a21**2*t11 + a11*a12*a21*a22*t12 + 
     1  a11*a13*a21*a23*t13 + a11*a12*a21**2*t14 + 
     1  a11**2*a21*a22*t14 + a11*a13*a21**2*t15 + 
     1  a11**2*a21*a23*t15 + a11*a13*a21*a22*t16 + 
     1  a11*a12*a21*a23*t16 + a11*a12*a21*a22*t21 + 
     1  a12**2*a22**2*t22 + a12*a13*a22*a23*t23 + 
     1  a12**2*a21*a22*t24 + a11*a12*a22**2*t24 + 
     1  a12*a13*a21*a22*t25 + a11*a12*a22*a23*t25 + 
     1  a12*a13*a22**2*t26 + a12**2*a22*a23*t26 + 
     1  a11*a13*a21*a23*t31 + a12*a13*a22*a23*t32 + 
     1  a13**2*a23**2*t33 + a12*a13*a21*a23*t34 + 
     1  a11*a13*a22*a23*t34 + a13**2*a21*a23*t35 + 
     1  a11*a13*a23**2*t35 + a13**2*a22*a23*t36 + 
     1  a12*a13*a23**2*t36 + a11*a12*a21**2*t41 + 
     1  a11**2*a21*a22*t41 + a12**2*a21*a22*t42 + 
     1  a11*a12*a22**2*t42 + a12*a13*a21*a23*t43 + 
     1  a11*a13*a22*a23*t43 + a12**2*a21**2*t44 + 
     1  DOS*a11*a12*a21*a22*t44 + a11**2*a22**2*t44 + 
     1  a12*a13*a21**2*t45 + a11*a13*a21*a22*t45 + 
     1  a11*a12*a21*a23*t45 + a11**2*a22*a23*t45 + 
     1  a12*a13*a21*a22*t46 + a11*a13*a22**2*t46 + 
     1  a12**2*a21*a23*t46 + a11*a12*a22*a23*t46 + 
     1  a11*a13*a21**2*t51 + a11**2*a21*a23*t51 + 
     1  a12*a13*a21*a22*t52 + a11*a12*a22*a23*t52 + 
     1  a13**2*a21*a23*t53 + a11*a13*a23**2*t53 + 
     1  a12*a13*a21**2*t54 + a11*a13*a21*a22*t54 + 
     1  a11*a12*a21*a23*t54 + a11**2*a22*a23*t54 + 
     1  a13**2*a21**2*t55 + DOS*a11*a13*a21*a23*t55 + 
     1  a11**2*a23**2*t55 + a13**2*a21*a22*t56 + 
     1  a12*a13*a21*a23*t56 + a11*a13*a22*a23*t56 + 
     1  a11*a12*a23**2*t56 + a11*a13*a21*a22*t61 + 
     1  a11*a12*a21*a23*t61 + a12*a13*a22**2*t62 + 
     1  a12**2*a22*a23*t62 + a13**2*a22*a23*t63 + 
     1  a12*a13*a23**2*t63 + a12*a13*a21*a22*t64 + 
     1  a11*a13*a22**2*t64 + a12**2*a21*a23*t64 + 
     1  a11*a12*a22*a23*t64 + a13**2*a21*a22*t65 + 
     1  a12*a13*a21*a23*t65 + a11*a13*a22*a23*t65 + 
     1  a11*a12*a23**2*t65 + a13**2*a22**2*t66 + 
     1  DOS*a12*a13*a22*a23*t66 + a12**2*a23**2*t66
       TENSOR_R(1,2,1,3)=a11**2*a21*a31*t11 + a11*a12*a21*a32*t12 + 
     1  a11*a13*a21*a33*t13 + a11*a12*a21*a31*t14 + 
     1  a11**2*a21*a32*t14 + a11*a13*a21*a31*t15 + 
     1  a11**2*a21*a33*t15 + a11*a13*a21*a32*t16 + 
     1  a11*a12*a21*a33*t16 + a11*a12*a22*a31*t21 + 
     1  a12**2*a22*a32*t22 + a12*a13*a22*a33*t23 + 
     1  a12**2*a22*a31*t24 + a11*a12*a22*a32*t24 + 
     1  a12*a13*a22*a31*t25 + a11*a12*a22*a33*t25 + 
     1  a12*a13*a22*a32*t26 + a12**2*a22*a33*t26 + 
     1  a11*a13*a23*a31*t31 + a12*a13*a23*a32*t32 + 
     1  a13**2*a23*a33*t33 + a12*a13*a23*a31*t34 + 
     1  a11*a13*a23*a32*t34 + a13**2*a23*a31*t35 + 
     1  a11*a13*a23*a33*t35 + a13**2*a23*a32*t36 + 
     1  a12*a13*a23*a33*t36 + a11*a12*a21*a31*t41 + 
     1  a11**2*a22*a31*t41 + a12**2*a21*a32*t42 + 
     1  a11*a12*a22*a32*t42 + a12*a13*a21*a33*t43 + 
     1  a11*a13*a22*a33*t43 + a12**2*a21*a31*t44 + 
     1  a11*a12*a22*a31*t44 + a11*a12*a21*a32*t44 + 
     1  a11**2*a22*a32*t44 + a12*a13*a21*a31*t45 + 
     1  a11*a13*a22*a31*t45 + a11*a12*a21*a33*t45 + 
     1  a11**2*a22*a33*t45 + a12*a13*a21*a32*t46 + 
     1  a11*a13*a22*a32*t46 + a12**2*a21*a33*t46 + 
     1  a11*a12*a22*a33*t46 + a11*a13*a21*a31*t51 + 
     1  a11**2*a23*a31*t51 + a12*a13*a21*a32*t52 + 
     1  a11*a12*a23*a32*t52 + a13**2*a21*a33*t53 + 
     1  a11*a13*a23*a33*t53 + a12*a13*a21*a31*t54 + 
     1  a11*a12*a23*a31*t54 + a11*a13*a21*a32*t54 + 
     1  a11**2*a23*a32*t54 + a13**2*a21*a31*t55 + 
     1  a11*a13*a23*a31*t55 + a11*a13*a21*a33*t55 + 
     1  a11**2*a23*a33*t55 + a13**2*a21*a32*t56 + 
     1  a11*a13*a23*a32*t56 + a12*a13*a21*a33*t56 + 
     1  a11*a12*a23*a33*t56 + a11*a13*a22*a31*t61 + 
     1  a11*a12*a23*a31*t61 + a12*a13*a22*a32*t62 + 
     1  a12**2*a23*a32*t62 + a13**2*a22*a33*t63 + 
     1  a12*a13*a23*a33*t63 + a12*a13*a22*a31*t64 + 
     1  a12**2*a23*a31*t64 + a11*a13*a22*a32*t64 + 
     1  a11*a12*a23*a32*t64 + a13**2*a22*a31*t65 + 
     1  a12*a13*a23*a31*t65 + a11*a13*a22*a33*t65 + 
     1  a11*a12*a23*a33*t65 + a13**2*a22*a32*t66 + 
     1  a12*a13*a23*a32*t66 + a12*a13*a22*a33*t66 + 
     1  a12**2*a23*a33*t66
       TENSOR_R(1,2,2,2)=a11*a21**3*t11 + a11*a21*a22**2*t12 + 
     1  a11*a21*a23**2*t13 + DOS*a11*a21**2*a22*t14 + 
     1  DOS*a11*a21**2*a23*t15 + DOS*a11*a21*a22*a23*t16 + 
     1  a12*a21**2*a22*t21 + a12*a22**3*t22 + 
     1  a12*a22*a23**2*t23 + DOS*a12*a21*a22**2*t24 + 
     1  DOS*a12*a21*a22*a23*t25 + DOS*a12*a22**2*a23*t26 + 
     1  a13*a21**2*a23*t31 + a13*a22**2*a23*t32 + 
     1  a13*a23**3*t33 + DOS*a13*a21*a22*a23*t34 + 
     1  DOS*a13*a21*a23**2*t35 + DOS*a13*a22*a23**2*t36 + 
     1  a12*a21**3*t41 + a11*a21**2*a22*t41 + 
     1  a12*a21*a22**2*t42 + a11*a22**3*t42 + 
     1  a12*a21*a23**2*t43 + a11*a22*a23**2*t43 + 
     1  DOS*a12*a21**2*a22*t44 + DOS*a11*a21*a22**2*t44 + 
     1  DOS*a12*a21**2*a23*t45 + DOS*a11*a21*a22*a23*t45 + 
     1  DOS*a12*a21*a22*a23*t46 + DOS*a11*a22**2*a23*t46 + 
     1  a13*a21**3*t51 + a11*a21**2*a23*t51 + 
     1  a13*a21*a22**2*t52 + a11*a22**2*a23*t52 + 
     1  a13*a21*a23**2*t53 + a11*a23**3*t53 + 
     1  DOS*a13*a21**2*a22*t54 + DOS*a11*a21*a22*a23*t54 + 
     1  DOS*a13*a21**2*a23*t55 + DOS*a11*a21*a23**2*t55 + 
     1  DOS*a13*a21*a22*a23*t56 + DOS*a11*a22*a23**2*t56 + 
     1  a13*a21**2*a22*t61 + a12*a21**2*a23*t61 + 
     1  a13*a22**3*t62 + a12*a22**2*a23*t62 + 
     1  a13*a22*a23**2*t63 + a12*a23**3*t63 + 
     1  DOS*a13*a21*a22**2*t64 + DOS*a12*a21*a22*a23*t64 + 
     1  DOS*a13*a21*a22*a23*t65 + DOS*a12*a21*a23**2*t65 + 
     1  DOS*a13*a22**2*a23*t66 + DOS*a12*a22*a23**2*t66
       TENSOR_R(1,2,2,3)=a11*a21**2*a31*t11 + a11*a21*a22*a32*t12 + 
     1  a11*a21*a23*a33*t13 + a11*a21*a22*a31*t14 + 
     1  a11*a21**2*a32*t14 + a11*a21*a23*a31*t15 + 
     1  a11*a21**2*a33*t15 + a11*a21*a23*a32*t16 + 
     1  a11*a21*a22*a33*t16 + a12*a21*a22*a31*t21 + 
     1  a12*a22**2*a32*t22 + a12*a22*a23*a33*t23 + 
     1  a12*a22**2*a31*t24 + a12*a21*a22*a32*t24 + 
     1  a12*a22*a23*a31*t25 + a12*a21*a22*a33*t25 + 
     1  a12*a22*a23*a32*t26 + a12*a22**2*a33*t26 + 
     1  a13*a21*a23*a31*t31 + a13*a22*a23*a32*t32 + 
     1  a13*a23**2*a33*t33 + a13*a22*a23*a31*t34 + 
     1  a13*a21*a23*a32*t34 + a13*a23**2*a31*t35 + 
     1  a13*a21*a23*a33*t35 + a13*a23**2*a32*t36 + 
     1  a13*a22*a23*a33*t36 + a12*a21**2*a31*t41 + 
     1  a11*a21*a22*a31*t41 + a12*a21*a22*a32*t42 + 
     1  a11*a22**2*a32*t42 + a12*a21*a23*a33*t43 + 
     1  a11*a22*a23*a33*t43 + a12*a21*a22*a31*t44 + 
     1  a11*a22**2*a31*t44 + a12*a21**2*a32*t44 + 
     1  a11*a21*a22*a32*t44 + a12*a21*a23*a31*t45 + 
     1  a11*a22*a23*a31*t45 + a12*a21**2*a33*t45 + 
     1  a11*a21*a22*a33*t45 + a12*a21*a23*a32*t46 + 
     1  a11*a22*a23*a32*t46 + a12*a21*a22*a33*t46 + 
     1  a11*a22**2*a33*t46 + a13*a21**2*a31*t51 + 
     1  a11*a21*a23*a31*t51 + a13*a21*a22*a32*t52 + 
     1  a11*a22*a23*a32*t52 + a13*a21*a23*a33*t53 + 
     1  a11*a23**2*a33*t53 + a13*a21*a22*a31*t54 + 
     1  a11*a22*a23*a31*t54 + a13*a21**2*a32*t54 + 
     1  a11*a21*a23*a32*t54 + a13*a21*a23*a31*t55 + 
     1  a11*a23**2*a31*t55 + a13*a21**2*a33*t55 + 
     1  a11*a21*a23*a33*t55 + a13*a21*a23*a32*t56 + 
     1  a11*a23**2*a32*t56 + a13*a21*a22*a33*t56 + 
     1  a11*a22*a23*a33*t56 + a13*a21*a22*a31*t61 + 
     1  a12*a21*a23*a31*t61 + a13*a22**2*a32*t62 + 
     1  a12*a22*a23*a32*t62 + a13*a22*a23*a33*t63 + 
     1  a12*a23**2*a33*t63 + a13*a22**2*a31*t64 + 
     1  a12*a22*a23*a31*t64 + a13*a21*a22*a32*t64 + 
     1  a12*a21*a23*a32*t64 + a13*a22*a23*a31*t65 + 
     1  a12*a23**2*a31*t65 + a13*a21*a22*a33*t65 + 
     1  a12*a21*a23*a33*t65 + a13*a22*a23*a32*t66 + 
     1  a12*a23**2*a32*t66 + a13*a22**2*a33*t66 + 
     1  a12*a22*a23*a33*t66
       TENSOR_R(1,2,3,3)=a11*a21*a31**2*t11 + a11*a21*a32**2*t12 + 
     1  a11*a21*a33**2*t13 + DOS*a11*a21*a31*a32*t14 + 
     1  DOS*a11*a21*a31*a33*t15 + DOS*a11*a21*a32*a33*t16 + 
     1  a12*a22*a31**2*t21 + a12*a22*a32**2*t22 + 
     1  a12*a22*a33**2*t23 + DOS*a12*a22*a31*a32*t24 + 
     1  DOS*a12*a22*a31*a33*t25 + DOS*a12*a22*a32*a33*t26 + 
     1  a13*a23*a31**2*t31 + a13*a23*a32**2*t32 + 
     1  a13*a23*a33**2*t33 + DOS*a13*a23*a31*a32*t34 + 
     1  DOS*a13*a23*a31*a33*t35 + DOS*a13*a23*a32*a33*t36 + 
     1  a12*a21*a31**2*t41 + a11*a22*a31**2*t41 + 
     1  a12*a21*a32**2*t42 + a11*a22*a32**2*t42 + 
     1  a12*a21*a33**2*t43 + a11*a22*a33**2*t43 + 
     1  DOS*a12*a21*a31*a32*t44 + DOS*a11*a22*a31*a32*t44 + 
     1  DOS*a12*a21*a31*a33*t45 + DOS*a11*a22*a31*a33*t45 + 
     1  DOS*a12*a21*a32*a33*t46 + DOS*a11*a22*a32*a33*t46 + 
     1  a13*a21*a31**2*t51 + a11*a23*a31**2*t51 + 
     1  a13*a21*a32**2*t52 + a11*a23*a32**2*t52 + 
     1  a13*a21*a33**2*t53 + a11*a23*a33**2*t53 + 
     1  DOS*a13*a21*a31*a32*t54 + DOS*a11*a23*a31*a32*t54 + 
     1  DOS*a13*a21*a31*a33*t55 + DOS*a11*a23*a31*a33*t55 + 
     1  DOS*a13*a21*a32*a33*t56 + DOS*a11*a23*a32*a33*t56 + 
     1  a13*a22*a31**2*t61 + a12*a23*a31**2*t61 + 
     1  a13*a22*a32**2*t62 + a12*a23*a32**2*t62 + 
     1  a13*a22*a33**2*t63 + a12*a23*a33**2*t63 + 
     1  DOS*a13*a22*a31*a32*t64 + DOS*a12*a23*a31*a32*t64 + 
     1  DOS*a13*a22*a31*a33*t65 + DOS*a12*a23*a31*a33*t65 + 
     1  DOS*a13*a22*a32*a33*t66 + DOS*a12*a23*a32*a33*t66
       TENSOR_R(1,3,1,1)=a11**3*a31*t11 + a11*a12**2*a31*t12 + 
     1  a11*a13**2*a31*t13 + DOS*a11**2*a12*a31*t14 + 
     1  DOS*a11**2*a13*a31*t15 + DOS*a11*a12*a13*a31*t16 + 
     1  a11**2*a12*a32*t21 + a12**3*a32*t22 + 
     1  a12*a13**2*a32*t23 + DOS*a11*a12**2*a32*t24 + 
     1  DOS*a11*a12*a13*a32*t25 + DOS*a12**2*a13*a32*t26 + 
     1  a11**2*a13*a33*t31 + a12**2*a13*a33*t32 + 
     1  a13**3*a33*t33 + DOS*a11*a12*a13*a33*t34 + 
     1  DOS*a11*a13**2*a33*t35 + DOS*a12*a13**2*a33*t36 + 
     1  a11**2*a12*a31*t41 + a11**3*a32*t41 + 
     1  a12**3*a31*t42 + a11*a12**2*a32*t42 + 
     1  a12*a13**2*a31*t43 + a11*a13**2*a32*t43 + 
     1  DOS*a11*a12**2*a31*t44 + DOS*a11**2*a12*a32*t44 + 
     1  DOS*a11*a12*a13*a31*t45 + DOS*a11**2*a13*a32*t45 + 
     1  DOS*a12**2*a13*a31*t46 + DOS*a11*a12*a13*a32*t46 + 
     1  a11**2*a13*a31*t51 + a11**3*a33*t51 + 
     1  a12**2*a13*a31*t52 + a11*a12**2*a33*t52 + 
     1  a13**3*a31*t53 + a11*a13**2*a33*t53 + 
     1  DOS*a11*a12*a13*a31*t54 + DOS*a11**2*a12*a33*t54 + 
     1  DOS*a11*a13**2*a31*t55 + DOS*a11**2*a13*a33*t55 + 
     1  DOS*a12*a13**2*a31*t56 + DOS*a11*a12*a13*a33*t56 + 
     1  a11**2*a13*a32*t61 + a11**2*a12*a33*t61 + 
     1  a12**2*a13*a32*t62 + a12**3*a33*t62 + 
     1  a13**3*a32*t63 + a12*a13**2*a33*t63 + 
     1  DOS*a11*a12*a13*a32*t64 + DOS*a11*a12**2*a33*t64 + 
     1  DOS*a11*a13**2*a32*t65 + DOS*a11*a12*a13*a33*t65 + 
     1  DOS*a12*a13**2*a32*t66 + DOS*a12**2*a13*a33*t66
       TENSOR_R(1,3,1,2)=a11**2*a21*a31*t11 + a11*a12*a22*a31*t12 + 
     1  a11*a13*a23*a31*t13 + a11*a12*a21*a31*t14 + 
     1  a11**2*a22*a31*t14 + a11*a13*a21*a31*t15 + 
     1  a11**2*a23*a31*t15 + a11*a13*a22*a31*t16 + 
     1  a11*a12*a23*a31*t16 + a11*a12*a21*a32*t21 + 
     1  a12**2*a22*a32*t22 + a12*a13*a23*a32*t23 + 
     1  a12**2*a21*a32*t24 + a11*a12*a22*a32*t24 + 
     1  a12*a13*a21*a32*t25 + a11*a12*a23*a32*t25 + 
     1  a12*a13*a22*a32*t26 + a12**2*a23*a32*t26 + 
     1  a11*a13*a21*a33*t31 + a12*a13*a22*a33*t32 + 
     1  a13**2*a23*a33*t33 + a12*a13*a21*a33*t34 + 
     1  a11*a13*a22*a33*t34 + a13**2*a21*a33*t35 + 
     1  a11*a13*a23*a33*t35 + a13**2*a22*a33*t36 + 
     1  a12*a13*a23*a33*t36 + a11*a12*a21*a31*t41 + 
     1  a11**2*a21*a32*t41 + a12**2*a22*a31*t42 + 
     1  a11*a12*a22*a32*t42 + a12*a13*a23*a31*t43 + 
     1  a11*a13*a23*a32*t43 + a12**2*a21*a31*t44 + 
     1  a11*a12*a22*a31*t44 + a11*a12*a21*a32*t44 + 
     1  a11**2*a22*a32*t44 + a12*a13*a21*a31*t45 + 
     1  a11*a12*a23*a31*t45 + a11*a13*a21*a32*t45 + 
     1  a11**2*a23*a32*t45 + a12*a13*a22*a31*t46 + 
     1  a12**2*a23*a31*t46 + a11*a13*a22*a32*t46 + 
     1  a11*a12*a23*a32*t46 + a11*a13*a21*a31*t51 + 
     1  a11**2*a21*a33*t51 + a12*a13*a22*a31*t52 + 
     1  a11*a12*a22*a33*t52 + a13**2*a23*a31*t53 + 
     1  a11*a13*a23*a33*t53 + a12*a13*a21*a31*t54 + 
     1  a11*a13*a22*a31*t54 + a11*a12*a21*a33*t54 + 
     1  a11**2*a22*a33*t54 + a13**2*a21*a31*t55 + 
     1  a11*a13*a23*a31*t55 + a11*a13*a21*a33*t55 + 
     1  a11**2*a23*a33*t55 + a13**2*a22*a31*t56 + 
     1  a12*a13*a23*a31*t56 + a11*a13*a22*a33*t56 + 
     1  a11*a12*a23*a33*t56 + a11*a13*a21*a32*t61 + 
     1  a11*a12*a21*a33*t61 + a12*a13*a22*a32*t62 + 
     1  a12**2*a22*a33*t62 + a13**2*a23*a32*t63 + 
     1  a12*a13*a23*a33*t63 + a12*a13*a21*a32*t64 + 
     1  a11*a13*a22*a32*t64 + a12**2*a21*a33*t64 + 
     1  a11*a12*a22*a33*t64 + a13**2*a21*a32*t65 + 
     1  a11*a13*a23*a32*t65 + a12*a13*a21*a33*t65 + 
     1  a11*a12*a23*a33*t65 + a13**2*a22*a32*t66 + 
     1  a12*a13*a23*a32*t66 + a12*a13*a22*a33*t66 + 
     1  a12**2*a23*a33*t66
       TENSOR_R(1,3,1,3)=a11**2*a31**2*t11 + a11*a12*a31*a32*t12 + 
     1  a11*a13*a31*a33*t13 + a11*a12*a31**2*t14 + 
     1  a11**2*a31*a32*t14 + a11*a13*a31**2*t15 + 
     1  a11**2*a31*a33*t15 + a11*a13*a31*a32*t16 + 
     1  a11*a12*a31*a33*t16 + a11*a12*a31*a32*t21 + 
     1  a12**2*a32**2*t22 + a12*a13*a32*a33*t23 + 
     1  a12**2*a31*a32*t24 + a11*a12*a32**2*t24 + 
     1  a12*a13*a31*a32*t25 + a11*a12*a32*a33*t25 + 
     1  a12*a13*a32**2*t26 + a12**2*a32*a33*t26 + 
     1  a11*a13*a31*a33*t31 + a12*a13*a32*a33*t32 + 
     1  a13**2*a33**2*t33 + a12*a13*a31*a33*t34 + 
     1  a11*a13*a32*a33*t34 + a13**2*a31*a33*t35 + 
     1  a11*a13*a33**2*t35 + a13**2*a32*a33*t36 + 
     1  a12*a13*a33**2*t36 + a11*a12*a31**2*t41 + 
     1  a11**2*a31*a32*t41 + a12**2*a31*a32*t42 + 
     1  a11*a12*a32**2*t42 + a12*a13*a31*a33*t43 + 
     1  a11*a13*a32*a33*t43 + a12**2*a31**2*t44 + 
     1  DOS*a11*a12*a31*a32*t44 + a11**2*a32**2*t44 + 
     1  a12*a13*a31**2*t45 + a11*a13*a31*a32*t45 + 
     1  a11*a12*a31*a33*t45 + a11**2*a32*a33*t45 + 
     1  a12*a13*a31*a32*t46 + a11*a13*a32**2*t46 + 
     1  a12**2*a31*a33*t46 + a11*a12*a32*a33*t46 + 
     1  a11*a13*a31**2*t51 + a11**2*a31*a33*t51 + 
     1  a12*a13*a31*a32*t52 + a11*a12*a32*a33*t52 + 
     1  a13**2*a31*a33*t53 + a11*a13*a33**2*t53 + 
     1  a12*a13*a31**2*t54 + a11*a13*a31*a32*t54 + 
     1  a11*a12*a31*a33*t54 + a11**2*a32*a33*t54 + 
     1  a13**2*a31**2*t55 + DOS*a11*a13*a31*a33*t55 + 
     1  a11**2*a33**2*t55 + a13**2*a31*a32*t56 + 
     1  a12*a13*a31*a33*t56 + a11*a13*a32*a33*t56 + 
     1  a11*a12*a33**2*t56 + a11*a13*a31*a32*t61 + 
     1  a11*a12*a31*a33*t61 + a12*a13*a32**2*t62 + 
     1  a12**2*a32*a33*t62 + a13**2*a32*a33*t63 + 
     1  a12*a13*a33**2*t63 + a12*a13*a31*a32*t64 + 
     1  a11*a13*a32**2*t64 + a12**2*a31*a33*t64 + 
     1  a11*a12*a32*a33*t64 + a13**2*a31*a32*t65 + 
     1  a12*a13*a31*a33*t65 + a11*a13*a32*a33*t65 + 
     1  a11*a12*a33**2*t65 + a13**2*a32**2*t66 + 
     1  DOS*a12*a13*a32*a33*t66 + a12**2*a33**2*t66
       TENSOR_R(1,3,2,2)=a11*a21**2*a31*t11 + a11*a22**2*a31*t12 + 
     1  a11*a23**2*a31*t13 + DOS*a11*a21*a22*a31*t14 + 
     1  DOS*a11*a21*a23*a31*t15 + DOS*a11*a22*a23*a31*t16 + 
     1  a12*a21**2*a32*t21 + a12*a22**2*a32*t22 + 
     1  a12*a23**2*a32*t23 + DOS*a12*a21*a22*a32*t24 + 
     1  DOS*a12*a21*a23*a32*t25 + DOS*a12*a22*a23*a32*t26 + 
     1  a13*a21**2*a33*t31 + a13*a22**2*a33*t32 + 
     1  a13*a23**2*a33*t33 + DOS*a13*a21*a22*a33*t34 + 
     1  DOS*a13*a21*a23*a33*t35 + DOS*a13*a22*a23*a33*t36 + 
     1  a12*a21**2*a31*t41 + a11*a21**2*a32*t41 + 
     1  a12*a22**2*a31*t42 + a11*a22**2*a32*t42 + 
     1  a12*a23**2*a31*t43 + a11*a23**2*a32*t43 + 
     1  DOS*a12*a21*a22*a31*t44 + DOS*a11*a21*a22*a32*t44 + 
     1  DOS*a12*a21*a23*a31*t45 + DOS*a11*a21*a23*a32*t45 + 
     1  DOS*a12*a22*a23*a31*t46 + DOS*a11*a22*a23*a32*t46 + 
     1  a13*a21**2*a31*t51 + a11*a21**2*a33*t51 + 
     1  a13*a22**2*a31*t52 + a11*a22**2*a33*t52 + 
     1  a13*a23**2*a31*t53 + a11*a23**2*a33*t53 + 
     1  DOS*a13*a21*a22*a31*t54 + DOS*a11*a21*a22*a33*t54 + 
     1  DOS*a13*a21*a23*a31*t55 + DOS*a11*a21*a23*a33*t55 + 
     1  DOS*a13*a22*a23*a31*t56 + DOS*a11*a22*a23*a33*t56 + 
     1  a13*a21**2*a32*t61 + a12*a21**2*a33*t61 + 
     1  a13*a22**2*a32*t62 + a12*a22**2*a33*t62 + 
     1  a13*a23**2*a32*t63 + a12*a23**2*a33*t63 + 
     1  DOS*a13*a21*a22*a32*t64 + DOS*a12*a21*a22*a33*t64 + 
     1  DOS*a13*a21*a23*a32*t65 + DOS*a12*a21*a23*a33*t65 + 
     1  DOS*a13*a22*a23*a32*t66 + DOS*a12*a22*a23*a33*t66
       TENSOR_R(1,3,2,3)=a11*a21*a31**2*t11 + a11*a22*a31*a32*t12 + 
     1  a11*a23*a31*a33*t13 + a11*a22*a31**2*t14 + 
     1  a11*a21*a31*a32*t14 + a11*a23*a31**2*t15 + 
     1  a11*a21*a31*a33*t15 + a11*a23*a31*a32*t16 + 
     1  a11*a22*a31*a33*t16 + a12*a21*a31*a32*t21 + 
     1  a12*a22*a32**2*t22 + a12*a23*a32*a33*t23 + 
     1  a12*a22*a31*a32*t24 + a12*a21*a32**2*t24 + 
     1  a12*a23*a31*a32*t25 + a12*a21*a32*a33*t25 + 
     1  a12*a23*a32**2*t26 + a12*a22*a32*a33*t26 + 
     1  a13*a21*a31*a33*t31 + a13*a22*a32*a33*t32 + 
     1  a13*a23*a33**2*t33 + a13*a22*a31*a33*t34 + 
     1  a13*a21*a32*a33*t34 + a13*a23*a31*a33*t35 + 
     1  a13*a21*a33**2*t35 + a13*a23*a32*a33*t36 + 
     1  a13*a22*a33**2*t36 + a12*a21*a31**2*t41 + 
     1  a11*a21*a31*a32*t41 + a12*a22*a31*a32*t42 + 
     1  a11*a22*a32**2*t42 + a12*a23*a31*a33*t43 + 
     1  a11*a23*a32*a33*t43 + a12*a22*a31**2*t44 + 
     1  a12*a21*a31*a32*t44 + a11*a22*a31*a32*t44 + 
     1  a11*a21*a32**2*t44 + a12*a23*a31**2*t45 + 
     1  a11*a23*a31*a32*t45 + a12*a21*a31*a33*t45 + 
     1  a11*a21*a32*a33*t45 + a12*a23*a31*a32*t46 + 
     1  a11*a23*a32**2*t46 + a12*a22*a31*a33*t46 + 
     1  a11*a22*a32*a33*t46 + a13*a21*a31**2*t51 + 
     1  a11*a21*a31*a33*t51 + a13*a22*a31*a32*t52 + 
     1  a11*a22*a32*a33*t52 + a13*a23*a31*a33*t53 + 
     1  a11*a23*a33**2*t53 + a13*a22*a31**2*t54 + 
     1  a13*a21*a31*a32*t54 + a11*a22*a31*a33*t54 + 
     1  a11*a21*a32*a33*t54 + a13*a23*a31**2*t55 + 
     1  a13*a21*a31*a33*t55 + a11*a23*a31*a33*t55 + 
     1  a11*a21*a33**2*t55 + a13*a23*a31*a32*t56 + 
     1  a13*a22*a31*a33*t56 + a11*a23*a32*a33*t56 + 
     1  a11*a22*a33**2*t56 + a13*a21*a31*a32*t61 + 
     1  a12*a21*a31*a33*t61 + a13*a22*a32**2*t62 + 
     1  a12*a22*a32*a33*t62 + a13*a23*a32*a33*t63 + 
     1  a12*a23*a33**2*t63 + a13*a22*a31*a32*t64 + 
     1  a13*a21*a32**2*t64 + a12*a22*a31*a33*t64 + 
     1  a12*a21*a32*a33*t64 + a13*a23*a31*a32*t65 + 
     1  a12*a23*a31*a33*t65 + a13*a21*a32*a33*t65 + 
     1  a12*a21*a33**2*t65 + a13*a23*a32**2*t66 + 
     1  a13*a22*a32*a33*t66 + a12*a23*a32*a33*t66 + 
     1  a12*a22*a33**2*t66
       TENSOR_R(1,3,3,3)=a11*a31**3*t11 + a11*a31*a32**2*t12 + 
     1  a11*a31*a33**2*t13 + DOS*a11*a31**2*a32*t14 + 
     1  DOS*a11*a31**2*a33*t15 + DOS*a11*a31*a32*a33*t16 + 
     1  a12*a31**2*a32*t21 + a12*a32**3*t22 + 
     1  a12*a32*a33**2*t23 + DOS*a12*a31*a32**2*t24 + 
     1  DOS*a12*a31*a32*a33*t25 + DOS*a12*a32**2*a33*t26 + 
     1  a13*a31**2*a33*t31 + a13*a32**2*a33*t32 + 
     1  a13*a33**3*t33 + DOS*a13*a31*a32*a33*t34 + 
     1  DOS*a13*a31*a33**2*t35 + DOS*a13*a32*a33**2*t36 + 
     1  a12*a31**3*t41 + a11*a31**2*a32*t41 + 
     1  a12*a31*a32**2*t42 + a11*a32**3*t42 + 
     1  a12*a31*a33**2*t43 + a11*a32*a33**2*t43 + 
     1  DOS*a12*a31**2*a32*t44 + DOS*a11*a31*a32**2*t44 + 
     1  DOS*a12*a31**2*a33*t45 + DOS*a11*a31*a32*a33*t45 + 
     1  DOS*a12*a31*a32*a33*t46 + DOS*a11*a32**2*a33*t46 + 
     1  a13*a31**3*t51 + a11*a31**2*a33*t51 + 
     1  a13*a31*a32**2*t52 + a11*a32**2*a33*t52 + 
     1  a13*a31*a33**2*t53 + a11*a33**3*t53 + 
     1  DOS*a13*a31**2*a32*t54 + DOS*a11*a31*a32*a33*t54 + 
     1  DOS*a13*a31**2*a33*t55 + DOS*a11*a31*a33**2*t55 + 
     1  DOS*a13*a31*a32*a33*t56 + DOS*a11*a32*a33**2*t56 + 
     1  a13*a31**2*a32*t61 + a12*a31**2*a33*t61 + 
     1  a13*a32**3*t62 + a12*a32**2*a33*t62 + 
     1  a13*a32*a33**2*t63 + a12*a33**3*t63 + 
     1  DOS*a13*a31*a32**2*t64 + DOS*a12*a31*a32*a33*t64 + 
     1  DOS*a13*a31*a32*a33*t65 + DOS*a12*a31*a33**2*t65 + 
     1  DOS*a13*a32**2*a33*t66 + DOS*a12*a32*a33**2*t66
       TENSOR_R(2,2,1,1)=a11**2*a21**2*t11 + a12**2*a21**2*t12 + 
     1  a13**2*a21**2*t13 + DOS*a11*a12*a21**2*t14 + 
     1  DOS*a11*a13*a21**2*t15 + DOS*a12*a13*a21**2*t16 + 
     1  a11**2*a22**2*t21 + a12**2*a22**2*t22 + 
     1  a13**2*a22**2*t23 + DOS*a11*a12*a22**2*t24 + 
     1  DOS*a11*a13*a22**2*t25 + DOS*a12*a13*a22**2*t26 + 
     1  a11**2*a23**2*t31 + a12**2*a23**2*t32 + 
     1  a13**2*a23**2*t33 + DOS*a11*a12*a23**2*t34 + 
     1  DOS*a11*a13*a23**2*t35 + DOS*a12*a13*a23**2*t36 + 
     1  DOS*a11**2*a21*a22*t41 + DOS*a12**2*a21*a22*t42 + 
     1  DOS*a13**2*a21*a22*t43 + CUA*a11*a12*a21*a22*t44 + 
     1  CUA*a11*a13*a21*a22*t45 + CUA*a12*a13*a21*a22*t46 + 
     1  DOS*a11**2*a21*a23*t51 + DOS*a12**2*a21*a23*t52 + 
     1  DOS*a13**2*a21*a23*t53 + CUA*a11*a12*a21*a23*t54 + 
     1  CUA*a11*a13*a21*a23*t55 + CUA*a12*a13*a21*a23*t56 + 
     1  DOS*a11**2*a22*a23*t61 + DOS*a12**2*a22*a23*t62 + 
     1  DOS*a13**2*a22*a23*t63 + CUA*a11*a12*a22*a23*t64 + 
     1  CUA*a11*a13*a22*a23*t65 + CUA*a12*a13*a22*a23*t66
       TENSOR_R(2,2,1,2)=a11*a21**3*t11 + a12*a21**2*a22*t12 + 
     1  a13*a21**2*a23*t13 + a12*a21**3*t14 + 
     1  a11*a21**2*a22*t14 + a13*a21**3*t15 + 
     1  a11*a21**2*a23*t15 + a13*a21**2*a22*t16 + 
     1  a12*a21**2*a23*t16 + a11*a21*a22**2*t21 + 
     1  a12*a22**3*t22 + a13*a22**2*a23*t23 + 
     1  a12*a21*a22**2*t24 + a11*a22**3*t24 + 
     1  a13*a21*a22**2*t25 + a11*a22**2*a23*t25 + 
     1  a13*a22**3*t26 + a12*a22**2*a23*t26 + 
     1  a11*a21*a23**2*t31 + a12*a22*a23**2*t32 + 
     1  a13*a23**3*t33 + a12*a21*a23**2*t34 + 
     1  a11*a22*a23**2*t34 + a13*a21*a23**2*t35 + 
     1  a11*a23**3*t35 + a13*a22*a23**2*t36 + 
     1  a12*a23**3*t36 + DOS*a11*a21**2*a22*t41 + 
     1  DOS*a12*a21*a22**2*t42 + DOS*a13*a21*a22*a23*t43 + 
     1  DOS*a12*a21**2*a22*t44 + DOS*a11*a21*a22**2*t44 + 
     1  DOS*a13*a21**2*a22*t45 + DOS*a11*a21*a22*a23*t45 + 
     1  DOS*a13*a21*a22**2*t46 + DOS*a12*a21*a22*a23*t46 + 
     1  DOS*a11*a21**2*a23*t51 + DOS*a12*a21*a22*a23*t52 + 
     1  DOS*a13*a21*a23**2*t53 + DOS*a12*a21**2*a23*t54 + 
     1  DOS*a11*a21*a22*a23*t54 + DOS*a13*a21**2*a23*t55 + 
     1  DOS*a11*a21*a23**2*t55 + DOS*a13*a21*a22*a23*t56 + 
     1  DOS*a12*a21*a23**2*t56 + DOS*a11*a21*a22*a23*t61 + 
     1  DOS*a12*a22**2*a23*t62 + DOS*a13*a22*a23**2*t63 + 
     1  DOS*a12*a21*a22*a23*t64 + DOS*a11*a22**2*a23*t64 + 
     1  DOS*a13*a21*a22*a23*t65 + DOS*a11*a22*a23**2*t65 + 
     1  DOS*a13*a22**2*a23*t66 + DOS*a12*a22*a23**2*t66
       TENSOR_R(2,2,1,3)=a11*a21**2*a31*t11 + a12*a21**2*a32*t12 + 
     1  a13*a21**2*a33*t13 + a12*a21**2*a31*t14 + 
     1  a11*a21**2*a32*t14 + a13*a21**2*a31*t15 + 
     1  a11*a21**2*a33*t15 + a13*a21**2*a32*t16 + 
     1  a12*a21**2*a33*t16 + a11*a22**2*a31*t21 + 
     1  a12*a22**2*a32*t22 + a13*a22**2*a33*t23 + 
     1  a12*a22**2*a31*t24 + a11*a22**2*a32*t24 + 
     1  a13*a22**2*a31*t25 + a11*a22**2*a33*t25 + 
     1  a13*a22**2*a32*t26 + a12*a22**2*a33*t26 + 
     1  a11*a23**2*a31*t31 + a12*a23**2*a32*t32 + 
     1  a13*a23**2*a33*t33 + a12*a23**2*a31*t34 + 
     1  a11*a23**2*a32*t34 + a13*a23**2*a31*t35 + 
     1  a11*a23**2*a33*t35 + a13*a23**2*a32*t36 + 
     1  a12*a23**2*a33*t36 + DOS*a11*a21*a22*a31*t41 + 
     1  DOS*a12*a21*a22*a32*t42 + DOS*a13*a21*a22*a33*t43 + 
     1  DOS*a12*a21*a22*a31*t44 + DOS*a11*a21*a22*a32*t44 + 
     1  DOS*a13*a21*a22*a31*t45 + DOS*a11*a21*a22*a33*t45 + 
     1  DOS*a13*a21*a22*a32*t46 + DOS*a12*a21*a22*a33*t46 + 
     1  DOS*a11*a21*a23*a31*t51 + DOS*a12*a21*a23*a32*t52 + 
     1  DOS*a13*a21*a23*a33*t53 + DOS*a12*a21*a23*a31*t54 + 
     1  DOS*a11*a21*a23*a32*t54 + DOS*a13*a21*a23*a31*t55 + 
     1  DOS*a11*a21*a23*a33*t55 + DOS*a13*a21*a23*a32*t56 + 
     1  DOS*a12*a21*a23*a33*t56 + DOS*a11*a22*a23*a31*t61 + 
     1  DOS*a12*a22*a23*a32*t62 + DOS*a13*a22*a23*a33*t63 + 
     1  DOS*a12*a22*a23*a31*t64 + DOS*a11*a22*a23*a32*t64 + 
     1  DOS*a13*a22*a23*a31*t65 + DOS*a11*a22*a23*a33*t65 + 
     1  DOS*a13*a22*a23*a32*t66 + DOS*a12*a22*a23*a33*t66
       TENSOR_R(2,2,2,2)=a21**4*t11 + a21**2*a22**2*t12 + 
     1  a21**2*a23**2*t13 + DOS*a21**3*a22*t14 + 
     1  DOS*a21**3*a23*t15 + DOS*a21**2*a22*a23*t16 + 
     1  a21**2*a22**2*t21 + a22**4*t22 + 
     1  a22**2*a23**2*t23 + DOS*a21*a22**3*t24 + 
     1  DOS*a21*a22**2*a23*t25 + DOS*a22**3*a23*t26 + 
     1  a21**2*a23**2*t31 + a22**2*a23**2*t32 + 
     1  a23**4*t33 + DOS*a21*a22*a23**2*t34 + 
     1  DOS*a21*a23**3*t35 + DOS*a22*a23**3*t36 + 
     1  DOS*a21**3*a22*t41 + DOS*a21*a22**3*t42 + 
     1  DOS*a21*a22*a23**2*t43 + CUA*a21**2*a22**2*t44 + 
     1  CUA*a21**2*a22*a23*t45 + CUA*a21*a22**2*a23*t46 + 
     1  DOS*a21**3*a23*t51 + DOS*a21*a22**2*a23*t52 + 
     1  DOS*a21*a23**3*t53 + CUA*a21**2*a22*a23*t54 + 
     1  CUA*a21**2*a23**2*t55 + CUA*a21*a22*a23**2*t56 + 
     1  DOS*a21**2*a22*a23*t61 + DOS*a22**3*a23*t62 + 
     1  DOS*a22*a23**3*t63 + CUA*a21*a22**2*a23*t64 + 
     1  CUA*a21*a22*a23**2*t65 + CUA*a22**2*a23**2*t66
       TENSOR_R(2,2,2,3)=a21**3*a31*t11 + a21**2*a22*a32*t12 + 
     1  a21**2*a23*a33*t13 + a21**2*a22*a31*t14 + 
     1  a21**3*a32*t14 + a21**2*a23*a31*t15 + 
     1  a21**3*a33*t15 + a21**2*a23*a32*t16 + 
     1  a21**2*a22*a33*t16 + a21*a22**2*a31*t21 + 
     1  a22**3*a32*t22 + a22**2*a23*a33*t23 + 
     1  a22**3*a31*t24 + a21*a22**2*a32*t24 + 
     1  a22**2*a23*a31*t25 + a21*a22**2*a33*t25 + 
     1  a22**2*a23*a32*t26 + a22**3*a33*t26 + 
     1  a21*a23**2*a31*t31 + a22*a23**2*a32*t32 + 
     1  a23**3*a33*t33 + a22*a23**2*a31*t34 + 
     1  a21*a23**2*a32*t34 + a23**3*a31*t35 + 
     1  a21*a23**2*a33*t35 + a23**3*a32*t36 + 
     1  a22*a23**2*a33*t36 + DOS*a21**2*a22*a31*t41 + 
     1  DOS*a21*a22**2*a32*t42 + DOS*a21*a22*a23*a33*t43 + 
     1  DOS*a21*a22**2*a31*t44 + DOS*a21**2*a22*a32*t44 + 
     1  DOS*a21*a22*a23*a31*t45 + DOS*a21**2*a22*a33*t45 + 
     1  DOS*a21*a22*a23*a32*t46 + DOS*a21*a22**2*a33*t46 + 
     1  DOS*a21**2*a23*a31*t51 + DOS*a21*a22*a23*a32*t52 + 
     1  DOS*a21*a23**2*a33*t53 + DOS*a21*a22*a23*a31*t54 + 
     1  DOS*a21**2*a23*a32*t54 + DOS*a21*a23**2*a31*t55 + 
     1  DOS*a21**2*a23*a33*t55 + DOS*a21*a23**2*a32*t56 + 
     1  DOS*a21*a22*a23*a33*t56 + DOS*a21*a22*a23*a31*t61 + 
     1  DOS*a22**2*a23*a32*t62 + DOS*a22*a23**2*a33*t63 + 
     1  DOS*a22**2*a23*a31*t64 + DOS*a21*a22*a23*a32*t64 + 
     1  DOS*a22*a23**2*a31*t65 + DOS*a21*a22*a23*a33*t65 + 
     1  DOS*a22*a23**2*a32*t66 + DOS*a22**2*a23*a33*t66
       TENSOR_R(2,2,3,3)=a21**2*a31**2*t11 + a21**2*a32**2*t12 + 
     1  a21**2*a33**2*t13 + DOS*a21**2*a31*a32*t14 + 
     1  DOS*a21**2*a31*a33*t15 + DOS*a21**2*a32*a33*t16 + 
     1  a22**2*a31**2*t21 + a22**2*a32**2*t22 + 
     1  a22**2*a33**2*t23 + DOS*a22**2*a31*a32*t24 + 
     1  DOS*a22**2*a31*a33*t25 + DOS*a22**2*a32*a33*t26 + 
     1  a23**2*a31**2*t31 + a23**2*a32**2*t32 + 
     1  a23**2*a33**2*t33 + DOS*a23**2*a31*a32*t34 + 
     1  DOS*a23**2*a31*a33*t35 + DOS*a23**2*a32*a33*t36 + 
     1  DOS*a21*a22*a31**2*t41 + DOS*a21*a22*a32**2*t42 + 
     1  DOS*a21*a22*a33**2*t43 + CUA*a21*a22*a31*a32*t44 + 
     1  CUA*a21*a22*a31*a33*t45 + CUA*a21*a22*a32*a33*t46 + 
     1  DOS*a21*a23*a31**2*t51 + DOS*a21*a23*a32**2*t52 + 
     1  DOS*a21*a23*a33**2*t53 + CUA*a21*a23*a31*a32*t54 + 
     1  CUA*a21*a23*a31*a33*t55 + CUA*a21*a23*a32*a33*t56 + 
     1  DOS*a22*a23*a31**2*t61 + DOS*a22*a23*a32**2*t62 + 
     1  DOS*a22*a23*a33**2*t63 + CUA*a22*a23*a31*a32*t64 + 
     1  CUA*a22*a23*a31*a33*t65 + CUA*a22*a23*a32*a33*t66
       TENSOR_R(2,3,1,1)=a11**2*a21*a31*t11 + a12**2*a21*a31*t12 + 
     1  a13**2*a21*a31*t13 + DOS*a11*a12*a21*a31*t14 + 
     1  DOS*a11*a13*a21*a31*t15 + DOS*a12*a13*a21*a31*t16 + 
     1  a11**2*a22*a32*t21 + a12**2*a22*a32*t22 + 
     1  a13**2*a22*a32*t23 + DOS*a11*a12*a22*a32*t24 + 
     1  DOS*a11*a13*a22*a32*t25 + DOS*a12*a13*a22*a32*t26 + 
     1  a11**2*a23*a33*t31 + a12**2*a23*a33*t32 + 
     1  a13**2*a23*a33*t33 + DOS*a11*a12*a23*a33*t34 + 
     1  DOS*a11*a13*a23*a33*t35 + DOS*a12*a13*a23*a33*t36 + 
     1  a11**2*a22*a31*t41 + a11**2*a21*a32*t41 + 
     1  a12**2*a22*a31*t42 + a12**2*a21*a32*t42 + 
     1  a13**2*a22*a31*t43 + a13**2*a21*a32*t43 + 
     1  DOS*a11*a12*a22*a31*t44 + DOS*a11*a12*a21*a32*t44 + 
     1  DOS*a11*a13*a22*a31*t45 + DOS*a11*a13*a21*a32*t45 + 
     1  DOS*a12*a13*a22*a31*t46 + DOS*a12*a13*a21*a32*t46 + 
     1  a11**2*a23*a31*t51 + a11**2*a21*a33*t51 + 
     1  a12**2*a23*a31*t52 + a12**2*a21*a33*t52 + 
     1  a13**2*a23*a31*t53 + a13**2*a21*a33*t53 + 
     1  DOS*a11*a12*a23*a31*t54 + DOS*a11*a12*a21*a33*t54 + 
     1  DOS*a11*a13*a23*a31*t55 + DOS*a11*a13*a21*a33*t55 + 
     1  DOS*a12*a13*a23*a31*t56 + DOS*a12*a13*a21*a33*t56 + 
     1  a11**2*a23*a32*t61 + a11**2*a22*a33*t61 + 
     1  a12**2*a23*a32*t62 + a12**2*a22*a33*t62 + 
     1  a13**2*a23*a32*t63 + a13**2*a22*a33*t63 + 
     1  DOS*a11*a12*a23*a32*t64 + DOS*a11*a12*a22*a33*t64 + 
     1  DOS*a11*a13*a23*a32*t65 + DOS*a11*a13*a22*a33*t65 + 
     1  DOS*a12*a13*a23*a32*t66 + DOS*a12*a13*a22*a33*t66
       TENSOR_R(2,3,1,2)=a11*a21**2*a31*t11 + a12*a21*a22*a31*t12 + 
     1  a13*a21*a23*a31*t13 + a12*a21**2*a31*t14 + 
     1  a11*a21*a22*a31*t14 + a13*a21**2*a31*t15 + 
     1  a11*a21*a23*a31*t15 + a13*a21*a22*a31*t16 + 
     1  a12*a21*a23*a31*t16 + a11*a21*a22*a32*t21 + 
     1  a12*a22**2*a32*t22 + a13*a22*a23*a32*t23 + 
     1  a12*a21*a22*a32*t24 + a11*a22**2*a32*t24 + 
     1  a13*a21*a22*a32*t25 + a11*a22*a23*a32*t25 + 
     1  a13*a22**2*a32*t26 + a12*a22*a23*a32*t26 + 
     1  a11*a21*a23*a33*t31 + a12*a22*a23*a33*t32 + 
     1  a13*a23**2*a33*t33 + a12*a21*a23*a33*t34 + 
     1  a11*a22*a23*a33*t34 + a13*a21*a23*a33*t35 + 
     1  a11*a23**2*a33*t35 + a13*a22*a23*a33*t36 + 
     1  a12*a23**2*a33*t36 + a11*a21*a22*a31*t41 + 
     1  a11*a21**2*a32*t41 + a12*a22**2*a31*t42 + 
     1  a12*a21*a22*a32*t42 + a13*a22*a23*a31*t43 + 
     1  a13*a21*a23*a32*t43 + a12*a21*a22*a31*t44 + 
     1  a11*a22**2*a31*t44 + a12*a21**2*a32*t44 + 
     1  a11*a21*a22*a32*t44 + a13*a21*a22*a31*t45 + 
     1  a11*a22*a23*a31*t45 + a13*a21**2*a32*t45 + 
     1  a11*a21*a23*a32*t45 + a13*a22**2*a31*t46 + 
     1  a12*a22*a23*a31*t46 + a13*a21*a22*a32*t46 + 
     1  a12*a21*a23*a32*t46 + a11*a21*a23*a31*t51 + 
     1  a11*a21**2*a33*t51 + a12*a22*a23*a31*t52 + 
     1  a12*a21*a22*a33*t52 + a13*a23**2*a31*t53 + 
     1  a13*a21*a23*a33*t53 + a12*a21*a23*a31*t54 + 
     1  a11*a22*a23*a31*t54 + a12*a21**2*a33*t54 + 
     1  a11*a21*a22*a33*t54 + a13*a21*a23*a31*t55 + 
     1  a11*a23**2*a31*t55 + a13*a21**2*a33*t55 + 
     1  a11*a21*a23*a33*t55 + a13*a22*a23*a31*t56 + 
     1  a12*a23**2*a31*t56 + a13*a21*a22*a33*t56 + 
     1  a12*a21*a23*a33*t56 + a11*a21*a23*a32*t61 + 
     1  a11*a21*a22*a33*t61 + a12*a22*a23*a32*t62 + 
     1  a12*a22**2*a33*t62 + a13*a23**2*a32*t63 + 
     1  a13*a22*a23*a33*t63 + a12*a21*a23*a32*t64 + 
     1  a11*a22*a23*a32*t64 + a12*a21*a22*a33*t64 + 
     1  a11*a22**2*a33*t64 + a13*a21*a23*a32*t65 + 
     1  a11*a23**2*a32*t65 + a13*a21*a22*a33*t65 + 
     1  a11*a22*a23*a33*t65 + a13*a22*a23*a32*t66 + 
     1  a12*a23**2*a32*t66 + a13*a22**2*a33*t66 + 
     1  a12*a22*a23*a33*t66
       TENSOR_R(2,3,1,3)=a11*a21*a31**2*t11 + a12*a21*a31*a32*t12 + 
     1  a13*a21*a31*a33*t13 + a12*a21*a31**2*t14 + 
     1  a11*a21*a31*a32*t14 + a13*a21*a31**2*t15 + 
     1  a11*a21*a31*a33*t15 + a13*a21*a31*a32*t16 + 
     1  a12*a21*a31*a33*t16 + a11*a22*a31*a32*t21 + 
     1  a12*a22*a32**2*t22 + a13*a22*a32*a33*t23 + 
     1  a12*a22*a31*a32*t24 + a11*a22*a32**2*t24 + 
     1  a13*a22*a31*a32*t25 + a11*a22*a32*a33*t25 + 
     1  a13*a22*a32**2*t26 + a12*a22*a32*a33*t26 + 
     1  a11*a23*a31*a33*t31 + a12*a23*a32*a33*t32 + 
     1  a13*a23*a33**2*t33 + a12*a23*a31*a33*t34 + 
     1  a11*a23*a32*a33*t34 + a13*a23*a31*a33*t35 + 
     1  a11*a23*a33**2*t35 + a13*a23*a32*a33*t36 + 
     1  a12*a23*a33**2*t36 + a11*a22*a31**2*t41 + 
     1  a11*a21*a31*a32*t41 + a12*a22*a31*a32*t42 + 
     1  a12*a21*a32**2*t42 + a13*a22*a31*a33*t43 + 
     1  a13*a21*a32*a33*t43 + a12*a22*a31**2*t44 + 
     1  a12*a21*a31*a32*t44 + a11*a22*a31*a32*t44 + 
     1  a11*a21*a32**2*t44 + a13*a22*a31**2*t45 + 
     1  a13*a21*a31*a32*t45 + a11*a22*a31*a33*t45 + 
     1  a11*a21*a32*a33*t45 + a13*a22*a31*a32*t46 + 
     1  a13*a21*a32**2*t46 + a12*a22*a31*a33*t46 + 
     1  a12*a21*a32*a33*t46 + a11*a23*a31**2*t51 + 
     1  a11*a21*a31*a33*t51 + a12*a23*a31*a32*t52 + 
     1  a12*a21*a32*a33*t52 + a13*a23*a31*a33*t53 + 
     1  a13*a21*a33**2*t53 + a12*a23*a31**2*t54 + 
     1  a11*a23*a31*a32*t54 + a12*a21*a31*a33*t54 + 
     1  a11*a21*a32*a33*t54 + a13*a23*a31**2*t55 + 
     1  a13*a21*a31*a33*t55 + a11*a23*a31*a33*t55 + 
     1  a11*a21*a33**2*t55 + a13*a23*a31*a32*t56 + 
     1  a12*a23*a31*a33*t56 + a13*a21*a32*a33*t56 + 
     1  a12*a21*a33**2*t56 + a11*a23*a31*a32*t61 + 
     1  a11*a22*a31*a33*t61 + a12*a23*a32**2*t62 + 
     1  a12*a22*a32*a33*t62 + a13*a23*a32*a33*t63 + 
     1  a13*a22*a33**2*t63 + a12*a23*a31*a32*t64 + 
     1  a11*a23*a32**2*t64 + a12*a22*a31*a33*t64 + 
     1  a11*a22*a32*a33*t64 + a13*a23*a31*a32*t65 + 
     1  a13*a22*a31*a33*t65 + a11*a23*a32*a33*t65 + 
     1  a11*a22*a33**2*t65 + a13*a23*a32**2*t66 + 
     1  a13*a22*a32*a33*t66 + a12*a23*a32*a33*t66 + 
     1  a12*a22*a33**2*t66
       TENSOR_R(2,3,2,2)=a21**3*a31*t11 + a21*a22**2*a31*t12 + 
     1  a21*a23**2*a31*t13 + DOS*a21**2*a22*a31*t14 + 
     1  DOS*a21**2*a23*a31*t15 + DOS*a21*a22*a23*a31*t16 + 
     1  a21**2*a22*a32*t21 + a22**3*a32*t22 + 
     1  a22*a23**2*a32*t23 + DOS*a21*a22**2*a32*t24 + 
     1  DOS*a21*a22*a23*a32*t25 + DOS*a22**2*a23*a32*t26 + 
     1  a21**2*a23*a33*t31 + a22**2*a23*a33*t32 + 
     1  a23**3*a33*t33 + DOS*a21*a22*a23*a33*t34 + 
     1  DOS*a21*a23**2*a33*t35 + DOS*a22*a23**2*a33*t36 + 
     1  a21**2*a22*a31*t41 + a21**3*a32*t41 + 
     1  a22**3*a31*t42 + a21*a22**2*a32*t42 + 
     1  a22*a23**2*a31*t43 + a21*a23**2*a32*t43 + 
     1  DOS*a21*a22**2*a31*t44 + DOS*a21**2*a22*a32*t44 + 
     1  DOS*a21*a22*a23*a31*t45 + DOS*a21**2*a23*a32*t45 + 
     1  DOS*a22**2*a23*a31*t46 + DOS*a21*a22*a23*a32*t46 + 
     1  a21**2*a23*a31*t51 + a21**3*a33*t51 + 
     1  a22**2*a23*a31*t52 + a21*a22**2*a33*t52 + 
     1  a23**3*a31*t53 + a21*a23**2*a33*t53 + 
     1  DOS*a21*a22*a23*a31*t54 + DOS*a21**2*a22*a33*t54 + 
     1  DOS*a21*a23**2*a31*t55 + DOS*a21**2*a23*a33*t55 + 
     1  DOS*a22*a23**2*a31*t56 + DOS*a21*a22*a23*a33*t56 + 
     1  a21**2*a23*a32*t61 + a21**2*a22*a33*t61 + 
     1  a22**2*a23*a32*t62 + a22**3*a33*t62 + 
     1  a23**3*a32*t63 + a22*a23**2*a33*t63 + 
     1  DOS*a21*a22*a23*a32*t64 + DOS*a21*a22**2*a33*t64 + 
     1  DOS*a21*a23**2*a32*t65 + DOS*a21*a22*a23*a33*t65 + 
     1  DOS*a22*a23**2*a32*t66 + DOS*a22**2*a23*a33*t66
       TENSOR_R(2,3,2,3)=a21**2*a31**2*t11 + a21*a22*a31*a32*t12 + 
     1  a21*a23*a31*a33*t13 + a21*a22*a31**2*t14 + 
     1  a21**2*a31*a32*t14 + a21*a23*a31**2*t15 + 
     1  a21**2*a31*a33*t15 + a21*a23*a31*a32*t16 + 
     1  a21*a22*a31*a33*t16 + a21*a22*a31*a32*t21 + 
     1  a22**2*a32**2*t22 + a22*a23*a32*a33*t23 + 
     1  a22**2*a31*a32*t24 + a21*a22*a32**2*t24 + 
     1  a22*a23*a31*a32*t25 + a21*a22*a32*a33*t25 + 
     1  a22*a23*a32**2*t26 + a22**2*a32*a33*t26 + 
     1  a21*a23*a31*a33*t31 + a22*a23*a32*a33*t32 + 
     1  a23**2*a33**2*t33 + a22*a23*a31*a33*t34 + 
     1  a21*a23*a32*a33*t34 + a23**2*a31*a33*t35 + 
     1  a21*a23*a33**2*t35 + a23**2*a32*a33*t36 + 
     1  a22*a23*a33**2*t36 + a21*a22*a31**2*t41 + 
     1  a21**2*a31*a32*t41 + a22**2*a31*a32*t42 + 
     1  a21*a22*a32**2*t42 + a22*a23*a31*a33*t43 + 
     1  a21*a23*a32*a33*t43 + a22**2*a31**2*t44 + 
     1  DOS*a21*a22*a31*a32*t44 + a21**2*a32**2*t44 + 
     1  a22*a23*a31**2*t45 + a21*a23*a31*a32*t45 + 
     1  a21*a22*a31*a33*t45 + a21**2*a32*a33*t45 + 
     1  a22*a23*a31*a32*t46 + a21*a23*a32**2*t46 + 
     1  a22**2*a31*a33*t46 + a21*a22*a32*a33*t46 + 
     1  a21*a23*a31**2*t51 + a21**2*a31*a33*t51 + 
     1  a22*a23*a31*a32*t52 + a21*a22*a32*a33*t52 + 
     1  a23**2*a31*a33*t53 + a21*a23*a33**2*t53 + 
     1  a22*a23*a31**2*t54 + a21*a23*a31*a32*t54 + 
     1  a21*a22*a31*a33*t54 + a21**2*a32*a33*t54 + 
     1  a23**2*a31**2*t55 + DOS*a21*a23*a31*a33*t55 + 
     1  a21**2*a33**2*t55 + a23**2*a31*a32*t56 + 
     1  a22*a23*a31*a33*t56 + a21*a23*a32*a33*t56 + 
     1  a21*a22*a33**2*t56 + a21*a23*a31*a32*t61 + 
     1  a21*a22*a31*a33*t61 + a22*a23*a32**2*t62 + 
     1  a22**2*a32*a33*t62 + a23**2*a32*a33*t63 + 
     1  a22*a23*a33**2*t63 + a22*a23*a31*a32*t64 + 
     1  a21*a23*a32**2*t64 + a22**2*a31*a33*t64 + 
     1  a21*a22*a32*a33*t64 + a23**2*a31*a32*t65 + 
     1  a22*a23*a31*a33*t65 + a21*a23*a32*a33*t65 + 
     1  a21*a22*a33**2*t65 + a23**2*a32**2*t66 + 
     1  DOS*a22*a23*a32*a33*t66 + a22**2*a33**2*t66
       TENSOR_R(2,3,3,3)=a21*a31**3*t11 + a21*a31*a32**2*t12 + 
     1  a21*a31*a33**2*t13 + DOS*a21*a31**2*a32*t14 + 
     1  DOS*a21*a31**2*a33*t15 + DOS*a21*a31*a32*a33*t16 + 
     1  a22*a31**2*a32*t21 + a22*a32**3*t22 + 
     1  a22*a32*a33**2*t23 + DOS*a22*a31*a32**2*t24 + 
     1  DOS*a22*a31*a32*a33*t25 + DOS*a22*a32**2*a33*t26 + 
     1  a23*a31**2*a33*t31 + a23*a32**2*a33*t32 + 
     1  a23*a33**3*t33 + DOS*a23*a31*a32*a33*t34 + 
     1  DOS*a23*a31*a33**2*t35 + DOS*a23*a32*a33**2*t36 + 
     1  a22*a31**3*t41 + a21*a31**2*a32*t41 + 
     1  a22*a31*a32**2*t42 + a21*a32**3*t42 + 
     1  a22*a31*a33**2*t43 + a21*a32*a33**2*t43 + 
     1  DOS*a22*a31**2*a32*t44 + DOS*a21*a31*a32**2*t44 + 
     1  DOS*a22*a31**2*a33*t45 + DOS*a21*a31*a32*a33*t45 + 
     1  DOS*a22*a31*a32*a33*t46 + DOS*a21*a32**2*a33*t46 + 
     1  a23*a31**3*t51 + a21*a31**2*a33*t51 + 
     1  a23*a31*a32**2*t52 + a21*a32**2*a33*t52 + 
     1  a23*a31*a33**2*t53 + a21*a33**3*t53 + 
     1  DOS*a23*a31**2*a32*t54 + DOS*a21*a31*a32*a33*t54 + 
     1  DOS*a23*a31**2*a33*t55 + DOS*a21*a31*a33**2*t55 + 
     1  DOS*a23*a31*a32*a33*t56 + DOS*a21*a32*a33**2*t56 + 
     1  a23*a31**2*a32*t61 + a22*a31**2*a33*t61 + 
     1  a23*a32**3*t62 + a22*a32**2*a33*t62 + 
     1  a23*a32*a33**2*t63 + a22*a33**3*t63 + 
     1  DOS*a23*a31*a32**2*t64 + DOS*a22*a31*a32*a33*t64 + 
     1  DOS*a23*a31*a32*a33*t65 + DOS*a22*a31*a33**2*t65 + 
     1  DOS*a23*a32**2*a33*t66 + DOS*a22*a32*a33**2*t66
       TENSOR_R(3,3,1,1)=a11**2*a31**2*t11 + a12**2*a31**2*t12 + 
     1  a13**2*a31**2*t13 + DOS*a11*a12*a31**2*t14 + 
     1  DOS*a11*a13*a31**2*t15 + DOS*a12*a13*a31**2*t16 + 
     1  a11**2*a32**2*t21 + a12**2*a32**2*t22 + 
     1  a13**2*a32**2*t23 + DOS*a11*a12*a32**2*t24 + 
     1  DOS*a11*a13*a32**2*t25 + DOS*a12*a13*a32**2*t26 + 
     1  a11**2*a33**2*t31 + a12**2*a33**2*t32 + 
     1  a13**2*a33**2*t33 + DOS*a11*a12*a33**2*t34 + 
     1  DOS*a11*a13*a33**2*t35 + DOS*a12*a13*a33**2*t36 + 
     1  DOS*a11**2*a31*a32*t41 + DOS*a12**2*a31*a32*t42 + 
     1  DOS*a13**2*a31*a32*t43 + CUA*a11*a12*a31*a32*t44 + 
     1  CUA*a11*a13*a31*a32*t45 + CUA*a12*a13*a31*a32*t46 + 
     1  DOS*a11**2*a31*a33*t51 + DOS*a12**2*a31*a33*t52 + 
     1  DOS*a13**2*a31*a33*t53 + CUA*a11*a12*a31*a33*t54 + 
     1  CUA*a11*a13*a31*a33*t55 + CUA*a12*a13*a31*a33*t56 + 
     1  DOS*a11**2*a32*a33*t61 + DOS*a12**2*a32*a33*t62 + 
     1  DOS*a13**2*a32*a33*t63 + CUA*a11*a12*a32*a33*t64 + 
     1  CUA*a11*a13*a32*a33*t65 + CUA*a12*a13*a32*a33*t66
       TENSOR_R(3,3,1,2)=a11*a21*a31**2*t11 + a12*a22*a31**2*t12 + 
     1  a13*a23*a31**2*t13 + a12*a21*a31**2*t14 + 
     1  a11*a22*a31**2*t14 + a13*a21*a31**2*t15 + 
     1  a11*a23*a31**2*t15 + a13*a22*a31**2*t16 + 
     1  a12*a23*a31**2*t16 + a11*a21*a32**2*t21 + 
     1  a12*a22*a32**2*t22 + a13*a23*a32**2*t23 + 
     1  a12*a21*a32**2*t24 + a11*a22*a32**2*t24 + 
     1  a13*a21*a32**2*t25 + a11*a23*a32**2*t25 + 
     1  a13*a22*a32**2*t26 + a12*a23*a32**2*t26 + 
     1  a11*a21*a33**2*t31 + a12*a22*a33**2*t32 + 
     1  a13*a23*a33**2*t33 + a12*a21*a33**2*t34 + 
     1  a11*a22*a33**2*t34 + a13*a21*a33**2*t35 + 
     1  a11*a23*a33**2*t35 + a13*a22*a33**2*t36 + 
     1  a12*a23*a33**2*t36 + DOS*a11*a21*a31*a32*t41 + 
     1  DOS*a12*a22*a31*a32*t42 + DOS*a13*a23*a31*a32*t43 + 
     1  DOS*a12*a21*a31*a32*t44 + DOS*a11*a22*a31*a32*t44 + 
     1  DOS*a13*a21*a31*a32*t45 + DOS*a11*a23*a31*a32*t45 + 
     1  DOS*a13*a22*a31*a32*t46 + DOS*a12*a23*a31*a32*t46 + 
     1  DOS*a11*a21*a31*a33*t51 + DOS*a12*a22*a31*a33*t52 + 
     1  DOS*a13*a23*a31*a33*t53 + DOS*a12*a21*a31*a33*t54 + 
     1  DOS*a11*a22*a31*a33*t54 + DOS*a13*a21*a31*a33*t55 + 
     1  DOS*a11*a23*a31*a33*t55 + DOS*a13*a22*a31*a33*t56 + 
     1  DOS*a12*a23*a31*a33*t56 + DOS*a11*a21*a32*a33*t61 + 
     1  DOS*a12*a22*a32*a33*t62 + DOS*a13*a23*a32*a33*t63 + 
     1  DOS*a12*a21*a32*a33*t64 + DOS*a11*a22*a32*a33*t64 + 
     1  DOS*a13*a21*a32*a33*t65 + DOS*a11*a23*a32*a33*t65 + 
     1  DOS*a13*a22*a32*a33*t66 + DOS*a12*a23*a32*a33*t66
       TENSOR_R(3,3,1,3)=a11*a31**3*t11 + a12*a31**2*a32*t12 + 
     1  a13*a31**2*a33*t13 + a12*a31**3*t14 + 
     1  a11*a31**2*a32*t14 + a13*a31**3*t15 + 
     1  a11*a31**2*a33*t15 + a13*a31**2*a32*t16 + 
     1  a12*a31**2*a33*t16 + a11*a31*a32**2*t21 + 
     1  a12*a32**3*t22 + a13*a32**2*a33*t23 + 
     1  a12*a31*a32**2*t24 + a11*a32**3*t24 + 
     1  a13*a31*a32**2*t25 + a11*a32**2*a33*t25 + 
     1  a13*a32**3*t26 + a12*a32**2*a33*t26 + 
     1  a11*a31*a33**2*t31 + a12*a32*a33**2*t32 + 
     1  a13*a33**3*t33 + a12*a31*a33**2*t34 + 
     1  a11*a32*a33**2*t34 + a13*a31*a33**2*t35 + 
     1  a11*a33**3*t35 + a13*a32*a33**2*t36 + 
     1  a12*a33**3*t36 + DOS*a11*a31**2*a32*t41 + 
     1  DOS*a12*a31*a32**2*t42 + DOS*a13*a31*a32*a33*t43 + 
     1  DOS*a12*a31**2*a32*t44 + DOS*a11*a31*a32**2*t44 + 
     1  DOS*a13*a31**2*a32*t45 + DOS*a11*a31*a32*a33*t45 + 
     1  DOS*a13*a31*a32**2*t46 + DOS*a12*a31*a32*a33*t46 + 
     1  DOS*a11*a31**2*a33*t51 + DOS*a12*a31*a32*a33*t52 + 
     1  DOS*a13*a31*a33**2*t53 + DOS*a12*a31**2*a33*t54 + 
     1  DOS*a11*a31*a32*a33*t54 + DOS*a13*a31**2*a33*t55 + 
     1  DOS*a11*a31*a33**2*t55 + DOS*a13*a31*a32*a33*t56 + 
     1  DOS*a12*a31*a33**2*t56 + DOS*a11*a31*a32*a33*t61 + 
     1  DOS*a12*a32**2*a33*t62 + DOS*a13*a32*a33**2*t63 + 
     1  DOS*a12*a31*a32*a33*t64 + DOS*a11*a32**2*a33*t64 + 
     1  DOS*a13*a31*a32*a33*t65 + DOS*a11*a32*a33**2*t65 + 
     1  DOS*a13*a32**2*a33*t66 + DOS*a12*a32*a33**2*t66
       TENSOR_R(3,3,2,2)=a21**2*a31**2*t11 + a22**2*a31**2*t12 + 
     1  a23**2*a31**2*t13 + DOS*a21*a22*a31**2*t14 + 
     1  DOS*a21*a23*a31**2*t15 + DOS*a22*a23*a31**2*t16 + 
     1  a21**2*a32**2*t21 + a22**2*a32**2*t22 + 
     1  a23**2*a32**2*t23 + DOS*a21*a22*a32**2*t24 + 
     1  DOS*a21*a23*a32**2*t25 + DOS*a22*a23*a32**2*t26 + 
     1  a21**2*a33**2*t31 + a22**2*a33**2*t32 + 
     1  a23**2*a33**2*t33 + DOS*a21*a22*a33**2*t34 + 
     1  DOS*a21*a23*a33**2*t35 + DOS*a22*a23*a33**2*t36 + 
     1  DOS*a21**2*a31*a32*t41 + DOS*a22**2*a31*a32*t42 + 
     1  DOS*a23**2*a31*a32*t43 + CUA*a21*a22*a31*a32*t44 + 
     1  CUA*a21*a23*a31*a32*t45 + CUA*a22*a23*a31*a32*t46 + 
     1  DOS*a21**2*a31*a33*t51 + DOS*a22**2*a31*a33*t52 + 
     1  DOS*a23**2*a31*a33*t53 + CUA*a21*a22*a31*a33*t54 + 
     1  CUA*a21*a23*a31*a33*t55 + CUA*a22*a23*a31*a33*t56 + 
     1  DOS*a21**2*a32*a33*t61 + DOS*a22**2*a32*a33*t62 + 
     1  DOS*a23**2*a32*a33*t63 + CUA*a21*a22*a32*a33*t64 + 
     1  CUA*a21*a23*a32*a33*t65 + CUA*a22*a23*a32*a33*t66
       TENSOR_R(3,3,2,3)=a21*a31**3*t11 + a22*a31**2*a32*t12 + 
     1  a23*a31**2*a33*t13 + a22*a31**3*t14 + 
     1  a21*a31**2*a32*t14 + a23*a31**3*t15 + 
     1  a21*a31**2*a33*t15 + a23*a31**2*a32*t16 + 
     1  a22*a31**2*a33*t16 + a21*a31*a32**2*t21 + 
     1  a22*a32**3*t22 + a23*a32**2*a33*t23 + 
     1  a22*a31*a32**2*t24 + a21*a32**3*t24 + 
     1  a23*a31*a32**2*t25 + a21*a32**2*a33*t25 + 
     1  a23*a32**3*t26 + a22*a32**2*a33*t26 + 
     1  a21*a31*a33**2*t31 + a22*a32*a33**2*t32 + 
     1  a23*a33**3*t33 + a22*a31*a33**2*t34 + 
     1  a21*a32*a33**2*t34 + a23*a31*a33**2*t35 + 
     1  a21*a33**3*t35 + a23*a32*a33**2*t36 + 
     1  a22*a33**3*t36 + DOS*a21*a31**2*a32*t41 + 
     1  DOS*a22*a31*a32**2*t42 + DOS*a23*a31*a32*a33*t43 + 
     1  DOS*a22*a31**2*a32*t44 + DOS*a21*a31*a32**2*t44 + 
     1  DOS*a23*a31**2*a32*t45 + DOS*a21*a31*a32*a33*t45 + 
     1  DOS*a23*a31*a32**2*t46 + DOS*a22*a31*a32*a33*t46 + 
     1  DOS*a21*a31**2*a33*t51 + DOS*a22*a31*a32*a33*t52 + 
     1  DOS*a23*a31*a33**2*t53 + DOS*a22*a31**2*a33*t54 + 
     1  DOS*a21*a31*a32*a33*t54 + DOS*a23*a31**2*a33*t55 + 
     1  DOS*a21*a31*a33**2*t55 + DOS*a23*a31*a32*a33*t56 + 
     1  DOS*a22*a31*a33**2*t56 + DOS*a21*a31*a32*a33*t61 + 
     1  DOS*a22*a32**2*a33*t62 + DOS*a23*a32*a33**2*t63 + 
     1  DOS*a22*a31*a32*a33*t64 + DOS*a21*a32**2*a33*t64 + 
     1  DOS*a23*a31*a32*a33*t65 + DOS*a21*a32*a33**2*t65 + 
     1  DOS*a23*a32**2*a33*t66 + DOS*a22*a32*a33**2*t66
       TENSOR_R(3,3,3,3)=a31**4*t11 + a31**2*a32**2*t12 + 
     1  a31**2*a33**2*t13 + DOS*a31**3*a32*t14 + 
     1  DOS*a31**3*a33*t15 + DOS*a31**2*a32*a33*t16 + 
     1  a31**2*a32**2*t21 + a32**4*t22 + 
     1  a32**2*a33**2*t23 + DOS*a31*a32**3*t24 + 
     1  DOS*a31*a32**2*a33*t25 + DOS*a32**3*a33*t26 + 
     1  a31**2*a33**2*t31 + a32**2*a33**2*t32 + 
     1  a33**4*t33 + DOS*a31*a32*a33**2*t34 + 
     1  DOS*a31*a33**3*t35 + DOS*a32*a33**3*t36 + 
     1  DOS*a31**3*a32*t41 + DOS*a31*a32**3*t42 + 
     1  DOS*a31*a32*a33**2*t43 + CUA*a31**2*a32**2*t44 + 
     1  CUA*a31**2*a32*a33*t45 + CUA*a31*a32**2*a33*t46 + 
     1  DOS*a31**3*a33*t51 + DOS*a31*a32**2*a33*t52 + 
     1  DOS*a31*a33**3*t53 + CUA*a31**2*a32*a33*t54 + 
     1  CUA*a31**2*a33**2*t55 + CUA*a31*a32*a33**2*t56 + 
     1  DOS*a31**2*a32*a33*t61 + DOS*a32**3*a33*t62 + 
     1  DOS*a32*a33**3*t63 + CUA*a31*a32**2*a33*t64 + 
     1  CUA*a31*a32*a33**2*t65 + CUA*a32**2*a33**2*t66

	TENSOR_R(1,1,2,1)=TENSOR_R(1,1,1,2)
	TENSOR_R(1,1,3,1)=TENSOR_R(1,1,1,3)
	TENSOR_R(1,1,3,2)=TENSOR_R(1,1,2,3)
	TENSOR_R(1,2,2,1)=TENSOR_R(1,2,1,2)
	TENSOR_R(1,2,3,1)=TENSOR_R(1,2,1,3)
	TENSOR_R(1,2,3,2)=TENSOR_R(1,2,2,3)
	TENSOR_R(1,3,2,1)=TENSOR_R(1,3,1,2)
	TENSOR_R(1,3,3,1)=TENSOR_R(1,3,1,3)
	TENSOR_R(1,3,3,2)=TENSOR_R(1,3,2,3)
	TENSOR_R(2,1,1,1)=TENSOR_R(1,2,1,1)
	TENSOR_R(2,1,1,2)=TENSOR_R(1,2,1,2)
	TENSOR_R(2,1,1,3)=TENSOR_R(1,2,1,3)
	TENSOR_R(2,1,2,1)=TENSOR_R(1,2,1,2)
	TENSOR_R(2,1,2,2)=TENSOR_R(1,2,2,2)
	TENSOR_R(2,1,2,3)=TENSOR_R(1,2,2,3)
	TENSOR_R(2,1,3,1)=TENSOR_R(1,2,1,3)
	TENSOR_R(2,1,3,2)=TENSOR_R(1,2,2,3)
	TENSOR_R(2,1,3,3)=TENSOR_R(1,2,3,3)
	TENSOR_R(2,2,2,1)=TENSOR_R(2,2,1,2)
	TENSOR_R(2,2,3,1)=TENSOR_R(2,2,1,3)
	TENSOR_R(2,2,3,2)=TENSOR_R(2,2,2,3)
	TENSOR_R(2,3,2,1)=TENSOR_R(2,3,1,2)
	TENSOR_R(2,3,3,1)=TENSOR_R(2,3,1,3)
	TENSOR_R(2,3,3,2)=TENSOR_R(2,3,2,3)
	TENSOR_R(3,1,1,1)=TENSOR_R(1,3,1,1)
	TENSOR_R(3,1,1,2)=TENSOR_R(1,3,1,2)
	TENSOR_R(3,1,1,3)=TENSOR_R(1,3,1,3)
	TENSOR_R(3,1,2,1)=TENSOR_R(1,3,1,2)
	TENSOR_R(3,1,2,2)=TENSOR_R(1,3,2,2)
	TENSOR_R(3,1,2,3)=TENSOR_R(1,3,2,3)
	TENSOR_R(3,1,3,1)=TENSOR_R(1,3,1,3)
	TENSOR_R(3,1,3,2)=TENSOR_R(1,3,2,3)
	TENSOR_R(3,1,3,3)=TENSOR_R(1,3,3,3)
	TENSOR_R(3,2,1,1)=TENSOR_R(2,3,1,1)
	TENSOR_R(3,2,1,2)=TENSOR_R(2,3,1,2)
	TENSOR_R(3,2,1,3)=TENSOR_R(2,3,1,3)
	TENSOR_R(3,2,2,1)=TENSOR_R(2,3,1,2)
	TENSOR_R(3,2,2,2)=TENSOR_R(2,3,2,2)
	TENSOR_R(3,2,2,3)=TENSOR_R(2,3,2,3)
	TENSOR_R(3,2,3,1)=TENSOR_R(2,3,1,3)
	TENSOR_R(3,2,3,2)=TENSOR_R(2,3,2,3)
	TENSOR_R(3,2,3,3)=TENSOR_R(2,3,3,3)
	TENSOR_R(3,3,2,1)=TENSOR_R(3,3,1,2)
	TENSOR_R(3,3,3,1)=TENSOR_R(3,3,1,3)
	TENSOR_R(3,3,3,2)=TENSOR_R(3,3,2,3)
 
	 RETURN
      END

c
c      SUBROUTINE STRAIN_G(STRAIN_GREEN,DEFGRAD_ELAS)
c      implicit none
c      REAL*8 STRAIN_GREEN(3,3),DEFGRAD_ELAS(3,3)
c      CALL TRANSPO



C     Polar decomposition by Simo 1991

      SUBROUTINE POLAR_DECOMP(UU,RR,CC,DEFGRAD)
      implicit none
      REAL*8 UU(3,3),UUINV(3,3),RR(3,3),DEFGRAD(3,3),DEFGRAD_T(3,3)
      REAL*8 CC(3,3),CC2(3,3),unit(3,3)
      REAL*8 I1,I2,I3,iu1,iu2,iu3
      REAL*8 x(3),lambda(3),lambda2(3)
      REAL*8 b,c,m,n,t,D,dd
      real*8 a,p,q,phi,temp1,y1,y2,y3,temp2,u,v,y2r,y2i
      INTEGER l,ii,jj,ndif
 
      real*8 pi


      PI=3.141592653589793238D0

      do,ii=1,3
         do,jj=1,3
            unit(ii,jj)=0d0
         enddo
         unit(ii,ii)=1d0
      enddo

c     Transpose F
      CALL TRANSPOSE(DEFGRAD_T,DEFGRAD,3,3)

c     CC=Ft*F
      CALL PMAT(CC,DEFGRAD_T,DEFGRAD,3,3,3)
c      write(*,*) 'cc',CC

C     Diagonalization of C (obtention of lambda(i))
      CALL INVARIANTS(I1,I2,I3,CC)

      b = I2 - I1*I1/3d0
      c = (-2d0/27d0)*I1*I1*I1 + I1*I2/3d0 - I3

       if(abs(b).LE.1D-15) then
         
         do,l=1,3
            x(l)=-c**(1D0/3d0) 
         enddo

c         write(*,*)'b<1e-12', x(1),x(2),x(3)
         
      else
         m  = 2d0*dsqrt(-b/3D0)
         n  = 3d0*c/(m*b)
         t  = datan2(dsqrt(abs(1d0-n*n)),n)/3D0
         do,l=1,3
            x(l) = m*dcos(t+2d0*(l-1)*PI/3d0)
         enddo
      endif

      do,l=1,3
         lambda(l)=dsqrt(x(l)+I1/3D0)
      enddo


c     Invariants of U

      iu1 = lambda(1) + lambda(2) + lambda(3)
      iu2 =  lambda(1)*lambda(2) + lambda(1)*lambda(3)
     1      +lambda(2)*lambda(3)
      iu3 = lambda(1) * lambda(2) * lambda(3)

      D= iu1*iu2 - iu3
      IF(ABS(D).le.1d-15) write(*,*)'ERROR D',D

      CALL PMAT(CC2,CC,CC,3,3,3)

      do,ii=1,3
         do,jj=1,3
            UU(ii,jj)=(1d0/D)*(-CC2(ii,jj)+(iu1*iu1-iu2)*CC(ii,jj)+
     1           iu1*iu3*unit(ii,jj))
            UUINV(ii,jj)=(1d0/iu3)*(CC(ii,jj)-iu1*UU(ii,jj)+
     1           iu2*unit(ii,jj))
         enddo
      enddo
      
C     ROTATION

      CALL PMAT(RR,DEFGRAD,UUINV,3,3,3)
            
c      write(*,*)'uu',uu
c      write(*,*)'rr',rr
      

      RETURN
      END


      SUBROUTINE DECOMP(Matrix,lambda,vec1,vec2,vec3)
      implicit none
      REAL*8 Matrix(3,3),lambda(3),vec1(3),vec2(3),vec3(3)
      REAL*8 I1,I2,I3
      real*8 a,d
      real*8 p,q,dd,phi,temp1,y1,y2,y3,temp2,u,v,y2r,y2i
      REAL*8 Matrix2(3,3)
      REAL*8 dnorm,det(3),detmax
      real*8 b,c,x(3),n,t,m,pi
      INTEGER ii,jj,ndif,nroot,i,ni,l



      PI=3.141592653589793238D0
C     Diagonalization of matrix (sym) (obtention of lambda(i))
      CALL INVARIANTS(I1,I2,I3,matrix)


c      write(*,*)'invariants',i1,i2,i3
      
      b = I2 - I1*I1/3.d0
      c = (-2/27.d0)*I1*I1*I1 + I1*I2/3.d0 - I3

c      write(*,*)'b,c',b,c

      if(abs(b).LE.1D-15) then
         
         do,l=1,3
            x(l)=-c**(1.D0/3.D0)
         enddo

         ndif=1
c         write(*,*) x(1),x(2),x(3)
         
      else
         
         m  = 2*dsqrt(-b/3.D0)
         n  = 3*c/(m*b)
         if(n.GT.1d0) THEN 
c            write(*,*)'RU error, SQRT() neg',1d0-n*n
            n=0d0
         endif
         t  = datan2(dsqrt(1d0-n*n),n)/3.D0
 
         if(abs(t).LE.1D-15) then
            ndif=2
         else
            ndif=3
         endif
c         write(*,*) m,n,t

         do,l=1,3
            x(l) = m*dcos(t+2*(l-1)*PI/3.d0)
         enddo

c         write(*,*) x(1),x(2),x(3)
         
      endif

      do,l=1,3
         lambda(l)=dsqrt(x(l)+I1/3D0)
      enddo

      write(*,*) 'lambdas',lambda(1),lambda(2),lambda(3)

      
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$      CALL INVARIANTS(I1,I2,I3,Matrix)
c$$$
c$$$c      write(*,*)'invariants',i1,i2,i3
c$$$      
c$$$      a=1d0
c$$$      b=-i1
c$$$      c=i2
c$$$      d=-i3
c$$$c      PI=3.141592653589793238D0
c$$$      pi=4d0*datan(1.d0)
c$$$
c$$$c Step 1: Calculate p and q --------------------------------------------
c$$$      p  = c/a - b*b/a/a/3D0
c$$$      q  = (2D0*b*b*b/a/a/a - 9D0*b*c/a/a + 27D0*d/a) / 27D0
c$$$
c$$$c Step 2: Calculate DD (discriminant) ----------------------------------
c$$$      DD = p*p*p/27D0 + q*q/4D0
c$$$c      write(*,*)'DD',dd
c$$$
c$$$c Step 3: Branch to different algorithms based on DD -------------------
c$$$      if(DD.lt.0.and.abs(dd/i3).GT.1D-15) then
c$$$c       Step 3b:
c$$$c       3 real unequal roots -- use the trigonometric formulation
c$$$        phi = dacos(-q/2D0/dsqrt(dabs(p*p*p)/27D0))
c$$$        temp1=2D0*dsqrt(dabs(p)/3D0)
c$$$        y1 =  temp1*dcos(phi/3D0)
c$$$        y2 = -temp1*dcos((phi+pi)/3D0)
c$$$        y3 = -temp1*dcos((phi-pi)/3D0)
c$$$      else
c$$$        if(abs(DD/i3).LE.1D-15) DD=0d0
c$$$c       Step 3a:
c$$$c       1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
c$$$        temp1 = -q/2D0 + dsqrt(DD)
c$$$        temp2 = -q/2D0 - dsqrt(DD)
c$$$        u = dabs(temp1)**(1D0/3D0)
c$$$        v = dabs(temp2)**(1D0/3D0)
c$$$        if(temp1 .lt. 0D0) u=-u
c$$$        if(temp2 .lt. 0D0) v=-v
c$$$        y1  = u + v
c$$$        y2r = -(u+v)/2D0
c$$$        y2i =  (u-v)*dsqrt(3D0)/2D0
c$$$      endif
c$$$
c$$$c Step 4: Final transformation -----------------------------------------
c$$$      temp1 = b/a/3D0
c$$$      y1 = y1-temp1
c$$$      y2 = y2-temp1
c$$$      y3 = y3-temp1
c$$$      y2r=y2r-temp1
c$$$
c$$$c Assign answers and number of lambdas differents-------------------------------------------------------
c$$$      if(DD.lt.0.and.abs(dd/i3).GT.1D-15) then
c$$$
c$$$        lambda(1) = y1
c$$$        lambda(2) = y2
c$$$        lambda(3) = y3
c$$$        ndif=3
c$$$       
c$$$      elseif(abs(DD/i3).LE.1D-15) then
c$$$
c$$$         if(abs(y1-y2r).LE.1D-15) then
c$$$            ndif=1 
c$$$            lambda(1)=y1
c$$$            lambda(2)=y2r
c$$$            lambda(3)=y2r
c$$$         else
c$$$            ndif=2
c$$$            if(y1.GT.y2r) THEN
c$$$               lambda(1)=y1
c$$$               lambda(2)=y2r
c$$$               lambda(3)=y2r
c$$$            else
c$$$               lambda(1)=y2r
c$$$               lambda(2)=y2r
c$$$               lambda(3)=y1
c$$$            endif
c$$$
c$$$         endif
c$$$
c$$$      else
c$$$      
c$$$         write(*,*)'error:COMPLEx'
c$$$      
c$$$      endif 
c$$$
c$$$c      write(*,*)'NDIF',ndif
c$$$c      write(*,*)'LAMBDAS',lambda(1),lambda(2),lambda(3)
c$$$
C     EIGENVECTORS

C     three different
      IF(ndif.EQ.3) THEN

c     FIRST
      do,ii=1,3
         do,jj=1,3
            Matrix2(ii,jj)=Matrix(ii,jj)
            if(ii.EQ.jj) MATRIX2(ii,jj)=MATRIX2(ii,jj)-lambda(1)
         enddo
      enddo

      det(1)=Matrix2(1,1)*Matrix2(2,2)-Matrix2(1,2)*Matrix2(2,1)
      det(2)=Matrix2(2,2)*Matrix2(3,3)-Matrix2(2,3)*Matrix2(3,2)
      det(3)=Matrix2(1,1)*Matrix2(3,3)-Matrix2(1,3)*Matrix2(3,1)
      detmax=0d0
      do,i=1,3
         if(abs(det(i)).GT.detmax) THEN
            ni=i
            detmax=abs(det(i))
         endif
      enddo
c      write(*,*)'dets',det
      if(ni.EQ.1) THEN
c         write(*,*)'supuesto1',detmax
         a=(1/det(ni))*
     1        (-Matrix2(2,2)*Matrix2(1,3)+Matrix2(2,1)*Matrix2(2,3))
         b=(1/det(ni))*
     1        (Matrix2(1,2)*Matrix2(1,3)-Matrix2(1,1)*Matrix2(2,3))
         c=1d0
      elseif(ni.EQ.2) THEN  
c         write(*,*)'supuesto2',detmax
         a=1d0
         b=(1/det(ni))*
     1        (-Matrix2(3,3)*Matrix2(2,1)+Matrix2(3,2)*Matrix2(3,1))
         c=(1/det(ni))*
     1        (Matrix2(2,3)*Matrix2(2,1)-Matrix2(2,2)*Matrix2(3,1))
      else
c         write(*,*)'supuesto3',detmax 
         b=1d0
         a=(1/det(ni))*
     1        (-Matrix2(3,3)*Matrix2(1,2)+Matrix2(3,1)*Matrix2(2,3))
         c=(1/det(ni))*
     1        (Matrix2(1,3)*Matrix2(1,2)-Matrix2(1,1)*Matrix2(2,3))
      endif
      
      dnorm=dsqrt(a*a+b*b+c*c)
      vec1(1)=a/dnorm
      vec1(2)=b/dnorm
      vec1(3)=c/dnorm
c     SECOND

      do,ii=1,3
         do,jj=1,3
            Matrix2(ii,jj)=Matrix(ii,jj)
            if(ii.EQ.jj) MATRIX2(ii,jj)=MATRIX2(ii,jj)-lambda(2)
         enddo
      enddo
      detmax=0d0
      det(1)=Matrix2(1,1)*Matrix2(2,2)-Matrix2(1,2)*Matrix2(2,1)
      det(2)=Matrix2(2,2)*Matrix2(3,3)-Matrix2(2,3)*Matrix2(3,2)
      det(3)=Matrix2(1,1)*Matrix2(3,3)-Matrix2(1,3)*Matrix2(3,1)
c      write(*,*)'dets',det
      do,i=1,3
         if(abs(det(i)).GT.detmax) THEN
            ni=i
            detmax=abs(det(i))
         endif
      enddo
      if(ni.EQ.1) THEN
c         write(*,*)'supuesto1',detmax
         a=(1/det(ni))*
     1        (-Matrix2(2,2)*Matrix2(1,3)+Matrix2(2,1)*Matrix2(2,3))
         b=(1/det(ni))*
     1        (Matrix2(1,2)*Matrix2(1,3)-Matrix2(1,1)*Matrix2(2,3))
         c=1d0
      elseif(ni.EQ.2) THEN  
c         write(*,*)'supuesto2',detmax
         a=1d0
         b=(1/det(ni))*
     1        (-Matrix2(3,3)*Matrix2(2,1)+Matrix2(3,2)*Matrix2(3,1))
         c=(1/det(ni))*
     1        (Matrix2(2,3)*Matrix2(2,1)-Matrix2(2,2)*Matrix2(3,1))
      else
c         write(*,*)'supuesto3',detmax 
         b=1d0
         a=(1/det(ni))*
     1        (-Matrix2(3,3)*Matrix2(1,2)+Matrix2(3,1)*Matrix2(2,3))
         c=(1/det(ni))*
     1        (Matrix2(1,3)*Matrix2(1,2)-Matrix2(1,1)*Matrix2(2,3))
      endif
      
      dnorm=dsqrt(a*a+b*b+c*c)
      vec2(1)=a/dnorm
      vec2(2)=b/dnorm
      vec2(3)=c/dnorm
         
c     THIRD
      CALL PVECT(vec3,vec1,vec2)
c      write(*,*)'v1:',vec1
c      write(*,*)'v2:',vec2
c      write(*,*)'v3:',vec3
      
C     TWO DIFFERENT EIGENVALUES
      ELSE IF(ndif.EQ.2) THEN

         if(y2r.GT.y1) then
            do,ii=1,3
               do,jj=1,3
                  Matrix2(ii,jj)=Matrix(ii,jj)
                  if(ii.EQ.jj) MATRIX2(ii,jj)=MATRIX2(ii,jj)-lambda(3)
               enddo
            enddo

         else
            
            do,ii=1,3
               do,jj=1,3
                  Matrix2(ii,jj)=Matrix(ii,jj)
                  if(ii.EQ.jj) MATRIX2(ii,jj)=MATRIX2(ii,jj)-lambda(1)
               enddo
            enddo

         endif

         detmax=0d0
         det(1)=Matrix2(1,1)*Matrix2(2,2)-Matrix2(1,2)*Matrix2(2,1)
         det(2)=Matrix2(2,2)*Matrix2(3,3)-Matrix2(2,3)*Matrix2(3,2)
         det(3)=Matrix2(1,1)*Matrix2(3,3)-Matrix2(1,3)*Matrix2(3,1)
         do,i=1,3
            if(abs(det(i)).GT.detmax) THEN
               ni=i
               detmax=abs(det(i))
            endif
         enddo
         if(ni.EQ.1) THEN
c            write(*,*)'supuesto1',detmax
            a=(1/det(ni))*
     1           (-Matrix2(2,2)*Matrix2(1,3)+Matrix2(2,1)*Matrix2(2,3))
            b=(1/det(ni))*
     1           (Matrix2(1,2)*Matrix2(1,3)-Matrix2(1,1)*Matrix2(2,3))
            c=1d0
         elseif(ni.EQ.2) THEN  
c            write(*,*)'supuesto2',detmax
            a=1d0
            b=(1/det(ni))*
     1           (-Matrix2(3,3)*Matrix2(2,1)+Matrix2(3,2)*Matrix2(3,1))
            c=(1/det(ni))*
     1           (Matrix2(2,3)*Matrix2(2,1)-Matrix2(2,2)*Matrix2(3,1))
         else
c            write(*,*)'supuesto3',detmax 
            b=1d0
            a=(1/det(ni))*
     1           (-Matrix2(3,3)*Matrix2(1,2)+Matrix2(3,1)*Matrix2(2,3))
            c=(1/det(ni))*
     1           (Matrix2(1,3)*Matrix2(1,2)-Matrix2(1,1)*Matrix2(2,3))
         endif
         
         dnorm=dsqrt(a*a+b*b+c*c)

         if(y2r.GT.y1) then
            vec3(1)=a/dnorm
            vec3(2)=b/dnorm
            vec3(3)=c/dnorm
            if(abs(vec3(3)).GT.1D-15) THEN
               
               vec1(1)=-vec3(3)/dsqrt(vec3(1)**2+vec3(3)**2)
               vec1(2)=0d0       
               vec1(3)=vec3(1)/dsqrt(vec3(1)**2+vec3(3)**2)
               
            else 
               vec1(1)=0d0
               vec1(2)=0d0
               vec1(3)=1d0
            endif
            call pvect(vec2,vec3,vec1)
           

         else

            vec1(1)=a/dnorm
            vec1(2)=b/dnorm
            vec1(3)=c/dnorm
            if(abs(vec1(3)).GT.1D-15) THEN
               vec2(1)=-vec1(3)/dsqrt(vec1(1)**2+vec1(3)**2)
               vec2(2)=0d0
               vec2(3)=vec1(1)/dsqrt(vec1(1)**2+vec1(3)**2)
            else 
               vec2(1)=0d0
               vec2(2)=0d0
               vec2(3)=1d0
            endif
            call pvect(vec3,vec1,vec2)
           

            
         endif
         
c     lambda1=lambda2=lambda3
      ELSE
         vec1(1)=1d0
         vec1(2)=0d0
         vec1(3)=0d0
         vec2(1)=0d0
         vec2(2)=1d0
         vec2(3)=0d0
         vec3(1)=0d0
         vec3(2)=0d0
         vec3(3)=1d0
         
      ENDIF
      
      RETURN 
      END
      

      SUBROUTINE LOGtens(matrix,lambda,vec1,vec2,vec3)
      REAL*8 matrix(3,3),lambda(3),vec1(3),vec2(3),vec3(3)
      REAL*8 TENS1(3,3),TENS2(3,3),TENS3(3,3)
      do,ii=1,3
         lambda(ii)=log(lambda(ii))
      enddo
      CALL PTENS(TENS1,vec1,vec1)
      CALL PTENS(TENS2,vec2,vec2)
      CALL PTENS(TENS3,vec3,vec3)
      do,ii=1,3
         do,jj=1,3
            matrix(ii,jj)=lambda(1)*TENS1(ii,jj)+
     1           lambda(2)*TENS2(ii,jj)+
     1           lambda(3)*TENS3(ii,jj)
         enddo
      enddo
      RETURN
      END

      SUBROUTINE INVARIANTS(I1,I2,I3,CC)
c     Calculate invariants of a symmetric matrix
      implicit none
      REAL*8 I1,I2,I3
      REAL*8 CC(3,3),CC2(3,3)
      real*8 trac2
  

      I1=CC(1,1)+CC(2,2)+CC(3,3)
      
      CALL PMAT(CC2,CC,CC,3,3,3)
      
      trac2=CC2(1,1)+CC2(2,2)+CC2(3,3)
      
      I2=0.5d0*(I1*I1-trac2)

      I3=CC(1,1)*CC(2,2)*CC(3,3)+CC(1,2)*CC(2,3)*CC(3,1)
     1  +CC(1,3)*CC(2,1)*CC(3,2)
     2  -CC(1,3)*CC(2,2)*CC(3,1)-CC(2,3)*CC(3,2)*CC(1,1)
     3  -CC(3,3)*CC(1,2)*CC(2,1)

      RETURN
      END

      SUBROUTINE PMAT(matout,mat1,mat2,id,jd,kd)
C     Matrix product
      REAL*8 MAT1,MAT2,MATOUT
      DIMENSION mat1(id,kd),mat2(kd,jd),matout(id,jd)
      do,i=1,id
         Do,j=1,jd
            matout(i,j)=0D0
            do,k=1,kd
               matout(i,j)=matout(i,j)+mat1(i,k)*mat2(k,j)
            enddo
         enddo
      enddo
      RETURN
      END
c
      SUBROUTINE PTENS(matout,vectors1,vectors2)
c     Tensorial product of 2 vectors
       integer i,j
       real*8 vectors1(3),vectors2(3),matout(3,3)
       do,i=1,3
          do,j=1,3
             matout(i,j)=vectors1(i)*vectors2(j)
          enddo
       enddo
       RETURN
       END
      
      SUBROUTINE PTENS2(matout4,tens1,tens2)
      integer i,j,k,l
      real*8 tens1(3,3),tens2(3,3),matout4(3,3,3,3)
      do,i=1,3
         do,j=1,3
            do,k=1,3
               do,l=1,3
                  matout4(i,j,k,l)=tens1(i,j)*tens2(k,l)
               enddo
            enddo
         enddo
      enddo
      return
      end
      

      SUBROUTINE TRANSPOSE(Mat2,Mat1,ifil,icol)
      REAL*8 MAT1,MAT2
      DIMENSION MAT1(ifil,icol),MAT2(icol,ifil)
      DO,i=1,ifil
         do,j=1,icol
            MAT2(j,i)=0D0
            MAT2(j,i)=MAT1(i,j)
C     write(8,*) j,i,'-->',i,j,MAT2(i,j)
         enddo
      enddo
      RETURN
      END

      SUBROUTINE PCONTRACT(pcont,MAT1,MAT2,n)
      REAL*8 MAT1(n,n),MAT2(n,n),pcont
      INTEGER n,ii,jj
      pcont=0d0
      do,ii=1,n
         do,jj=1,n
            pcont=pcont+MAT1(ii,jj)*MAT2(ii,jj)
         enddo
      enddo
      RETURN 
      END

      SUBROUTINE PCONTRACT2(pcont,MAT1,MAT2,n)
      REAL*8 MAT1(n,n,n,n),MAT2(n,n),pcont(n,n)
      INTEGER n,ii,jj,kk,ll
      do,ii=1,n
         do,jj=1,n
            pcont(ii,jj)=0d0
            do,kk=1,n
               do,ll=1,n
                  pcont(ii,jj)=pcont(ii,jj)+MAT1(ii,jj,kk,ll)*
     1                 MAT2(kk,ll)
               enddo
            enddo
         enddo
      enddo


      RETURN 
      END


      SUBROUTINE PVECT(vecout,vectors1,vectors2)
C     Vectorial product
      REAL*8 vectors1,vectors2,vecout
      DIMENSION vectors1(3),vectors2(3),vecout(3)
      Dimension ii(3,2)
      ii(1,1)=2
      ii(1,2)=3
      ii(2,1)=1
      ii(2,2)=3
      ii(3,1)=1
      ii(3,2)=2
      Do,i=1,3
         vecout(i)=0d0
         vecout(i)=(-1)**(i+1)*( vectors1(ii(i,1))*vectors2(ii(i,2))-
     1   vectors2(ii(i,1))*vectors1(ii(i,2)) )
      ENDDO
C      write(8,*) vecout(1),vecout(2),vecout(3)
      RETURN
      END




c
      SUBROUTINE MATINV3(CC2,CC,iflag)
      REAL*8 CC(3,3),CC2(3,3)
      REAL*8 I3
      integer iflag

      I3=CC(1,1)*CC(2,2)*CC(3,3)+CC(1,2)*CC(2,3)*CC(3,1)
     1     +CC(1,3)*CC(2,1)*CC(3,2)
     2     -CC(1,3)*CC(2,2)*CC(3,1)-CC(2,3)*CC(3,2)*CC(1,1)
     3     -CC(3,3)*CC(1,2)*CC(2,1)

      if(abs(I3).GT.1D-50) THEN
         iflag=0
         I3=1D0/I3
      else
         iflag=1
      endif

      CC2(1,1)=I3*(cc(2,2)*cc(3,3)-cc(2,3)*cc(3,2))
      CC2(2,2)=I3*(cc(1,1)*cc(3,3)-cc(1,3)*cc(3,1))
      CC2(3,3)=I3*(cc(1,1)*cc(2,2)-cc(1,2)*cc(2,1))

      CC2(1,2)=-I3*(cc(1,2)*cc(3,3)-cc(3,2)*cc(1,3))
      CC2(2,1)=-I3*(cc(2,1)*cc(3,3)-cc(2,3)*cc(3,1))

      CC2(1,3)=I3*(cc(1,2)*cc(2,3)-cc(1,3)*cc(2,2))
      CC2(3,1)=I3*(cc(2,1)*cc(3,2)-cc(2,2)*cc(3,1))

      CC2(2,3)=-I3*(cc(1,1)*cc(2,3)-cc(1,3)*cc(2,1))
      CC2(3,2)=-I3*(cc(1,1)*cc(3,2)-cc(1,2)*cc(3,1))

      RETURN 
      END


      SUBROUTINE MATINV (AA,LDA,N,nmax,IFLAG) 
     
C
C-----------------------------------------------------------------------
C   MATINV   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
C            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
C            MARYLAND  20899
C
C   FOR: COMPUTING THE INVERSE OF A GENERAL N BY N MATRIX IN PLACE,
C        I.E., THE INVERSE OVERWRITES THE ORIGINAL MATRIX.  THE STEPS 
C        OF THE ALGORITHM ARE DESCRIBED BELOW AS THEY OCCUR.  ROW
C        INTERCHANGES ARE DONE AS NEEDED IN ORDER TO INCREASE THE
C        ACCURACY OF THE INVERSE MATRIX.  WITHOUT INTERCHANGES THIS
C        ALGORITHM WILL FAIL WHEN ANY OF THE LEADING PRINCIPAL
C        SUBMATRICES ARE SINGULAR OR WHEN THE MATRIX ITSELF IS
C        SINGULAR.  WITH INTERCHANGES THIS ALGORITHM WILL FAIL ONLY
C        WHEN THE MATRIX ITSELF IS SINGULAR.  THE LEADING PRINCIPAL
C
C                                   [A B C]
C        SUBMATRICES OF THE MATRIX  [D E F]  ARE  [A]  AND  [A B] .
C                                   [G H I]                 [D E]
C
C   SUBPROGRAMS CALLED: -NONE-
C
C   CURRENT VERSION COMPLETED JANUARY 15, 1987
C
C   REFERENCE: STEWART, G.W., 'INTRODUCTION TO MATRIX COMPUTATIONS',
C              ACADEMIC PRESS, INC., 1973
C-----------------------------------------------------------------------
C   DEFINITION OF PASSED PARAMETERS
C
C     * A = MATRIX (SIZE NXN) TO BE INVERTED (REAL)
C
C   * LDA = LEADING DIMENSION OF MATRIX A [LDA>=N] (INTEGER)
C
C     * N = NUMBER OF ROWS AND COLUMNS OF MATRIX A (INTEGER)
C
C   IFLAG = ERROR INDICATOR ON OUTPUT (INTEGER)   INTERPRETATION: 
C           -2 -> TOO MANY ROW INTERCHANGES NEEDED - INCREASE MX
C           -1 -> N>LDA
C            0 -> NO ERRORS DETECTED
C            K -> MATRIX A FOUND TO BE SINGULAR AT THE KTH STEP OF
C                 THE CROUT REDUCTION (1<=K<=N)
C
C   * INDICATES PARAMETERS REQUIRING INPUT VALUES 
C-----------------------------------------------------------------------
      REAL*8 AA,A,Q,S,R
      PARAMETER (MX=100)
      DIMENSION AA(nmax,nmax)
   
      DIMENSION A(LDA,N),IEX(MX,2)
      IFLAG = 0
C
C--- CHECK CONSISTENCY OF PASSED PARAMETERS
C
      IF (N.GT.LDA) THEN
         IFLAG = -1 
         RETURN
      ENDIF

      do,ii=1,lda
         do,jj=1,lda
            A(ii,jj)=AA(ii,jj)
         enddo
      enddo
C     enddo
C--- COMPUTE A = LU BY THE CROUT REDUCTION WHERE L IS LOWER TRIANGULAR
C--- AND U IS UNIT UPPER TRIANGULAR (ALGORITHM 3.4, P. 138 OF THE
C--- REFERENCE)
C
      NEX = 0
      DO 70 K = 1, N
         DO 20 I = K, N
            S = A(I,K)
            DO 10 L = 1, K-1
               S = S-A(I,L)*A(L,K)
   10       CONTINUE
            A(I,K) = S
   20    CONTINUE
C
C--- INTERCHANGE ROWS IF NECESSARY
C
         Q = 0.0
         L = 0
         DO 30 I = K, N
            R = ABS(A(I,K))
            IF (R.GT.Q) THEN
               Q = R
               L = I
            ENDIF
   30    CONTINUE
         IF (L.EQ.0) THEN
            IFLAG = K
            RETURN
         ENDIF
         IF (L.NE.K) THEN
            NEX = NEX+1
            IF (NEX.GT.MX) THEN
               IFLAG = -2
               RETURN
            ENDIF
            IEX(NEX,1) = K
            IEX(NEX,2) = L
            DO 40 J = 1, N
               Q = A(K,J)
               A(K,J) = A(L,J)
               A(L,J) = Q
   40       CONTINUE
         ENDIF
C
C--- END ROW INTERCHANGE SECTION
C
         DO 60 J = K+1, N
            S = A(K,J)
            DO 50 L = 1, K-1
               S = S-A(K,L)*A(L,J)
   50       CONTINUE
            A(K,J) = S/A(K,K) 
   60    CONTINUE
   70 CONTINUE
C
C--- INVERT THE LOWER TRIANGLE L IN PLACE (SIMILAR TO ALGORITHM 1.5,
C--- P. 110 OF THE REFERENCE) 
C
      DO 100 K = N, 1, -1
         A(K,K) = 1.0/A(K,K)
         DO 90 I = K-1, 1, -1 
            S = 0.0 
            DO 80 J = I+1, K
               S = S+A(J,I)*A(K,J)
   80       CONTINUE
            A(K,I) = -S/A(I,I)
   90    CONTINUE
  100 CONTINUE
C
C--- INVERT THE UPPER TRIANGLE U IN PLACE (ALGORITHM 1.5, P. 110 OF
C--- THE REFERENCE) 
C
      DO 130 K = N, 1, -1
         DO 120 I = K-1, 1, -1
            S = A(I,K)
            DO 110 J = I+1, K-1
               S = S+A(I,J)*A(J,K)
  110       CONTINUE
            A(I,K) = -S
  120    CONTINUE
  130 CONTINUE
C
C--- COMPUTE INV(A) = INV(U)*INV(L)
C
      DO 160 I = 1, N
         DO 150 J = 1, N
            IF (J.GT.I) THEN
               S = 0.0
               L = J
            ELSE
               S = A(I,J)
               L = I+1
            ENDIF
            DO 140 K = L, N
               S = S+A(I,K)*A(K,J)
  140       CONTINUE
            A(I,J) = S
  150    CONTINUE
  160 CONTINUE
C
C--- INTERCHANGE COLUMNS OF INV(A) TO REVERSE EFFECT OF ROW 
C--- INTERCHANGES OF A
C
      DO 180 I = NEX, 1, -1
         K = IEX(I,1)
         L = IEX(I,2)
         DO 170 J = 1, N
            Q = A(J,K)
            A(J,K) = A(J,L)
            A(J,L) = Q
  170    CONTINUE
  180 CONTINUE

       do,ii=1,lda
         do,jj=1,lda
            AA(ii,jj)=A(ii,jj)
         enddo
      enddo

      RETURN
      END 

      SUBROUTINE EXPONENTIAL3(EXPO,MATRIX)
      REAL*8 EXPO(3,3),MATRIX(3,3)
      REAL*8 AUX33(3,3),MATRIXn(3,3)
      REAL*8 dnorm,toler
      toler=1D-5
      do,ii=1,3
         do,jj=1,3
            EXPO(ii,jj)=0d0
            MATRIXn(ii,jj)=0d0
         enddo
         EXPO(ii,ii)=1d0 
         MATRIXn(ii,ii)=1d0
      enddo
      n=0
      nfact=1
      dnorm=1D0
      DO 100 WHILE(dnorm.GE.toler)
         n=n+1
         nfact=nfact*n
         CALL PMAT(AUX33,MATRIX,MATRIXn,3,3,3)
         do,ii=1,3
            do,jj=1,3
               MATRIXn(ii,jj)=AUX33(ii,jj)
               EXPO(ii,jj)=EXPO(ii,jj)+(1D0/nfact)*MATRIXn(ii,jj)
            enddo
         enddo
c         write(*,*)'n,fact',n,nfact
c         write(*,*)'MATRIXn'
c         write(*,102) MATRIXn
c         write(*,*)'EXPO'
c        write(*,102) EXPO
         dnorm=0d0
         do,ii=1,3
            do,jj=1,3
               dnorm=dnorm+MATRIXn(ii,jj)*MATRIXn(ii,jj)
            enddo
         enddo
         dnorm=dsqrt(dnorm)/nfact
c         write(*,*) dnorm
 100  CONTINUE
      RETURN
 102  FORMAT(3(3(E16.4,' '),/))
      
      END

      SUBROUTINE EXPMAP
     1     (   EXPX       ,NOCONV     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1     (   NDIM=3     ,NDIM2=9    )
C     Arguments
      LOGICAL  NOCONV
      DIMENSION
     1     EXPX(NDIM,NDIM)           ,X(NDIM,NDIM)
C     Local arrays and variables
      DIMENSION
     1     XN(NDIM,NDIM)     ,XNM1(NDIM,NDIM)     ,XNM2(NDIM,NDIM)     ,
     2     XNM3(NDIM,NDIM)   ,X2(NDIM,NDIM)
      DATA
     1     R0   ,RP5  ,R1   ,R2   ,TOL    ,OVER    ,UNDER   /
     2     0.0D0,0.5D0,1.0D0,2.0D0,1.0D-5,1.0D+100,1.0D-100/
      DATA
     1     NMAX / 100 /
C***********************************************************************
C     COMPUTES THE EXPONENTIAL OF A (GENERALLY UNSYMMETRIC) 3-D TENSOR. USES
C     THE SERIES REPRESENTATION OF THE TENSOR EXPONENTIAL
C     
C     REFERENCE: Section B.1
C     Box B.1
C***********************************************************************
C     Initialise series convergence flag
      NOCONV=.FALSE.
C     Compute X square
      CALL RVZERO(X2,NDIM2)
      DO 30 I=1,NDIM
         DO 20 J=1,NDIM
            DO 10 K=1,NDIM
               X2(I,J)= X2(I,J)+X(I,K)*X(K,J)
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
C     Compute principal invariants of X
      C1=X(1,1)+X(2,2)+X(3,3)
      C2=RP5*(C1*C1-(X2(1,1)+X2(2,2)+X2(3,3)))
      C3=X(1,1)*X(2,2)*X(3,3)+X(1,2)*X(2,3)*X(3,1)+
     1     X(1,3)*X(2,1)*X(3,2)-X(1,2)*X(2,1)*X(3,3)-
     2     X(1,1)*X(2,3)*X(3,2)-X(1,3)*X(2,2)*X(3,1)
C     Start computation of exponential using its series definition
C     ============================================================
      DO 50 I=1,NDIM
         DO 40 J=1,NDIM
            XNM1(I,J)=X2(I,J)
            XNM2(I,J)=X(I,J)
 40      CONTINUE
 50   CONTINUE
      XNM3(1,1)=R1
      XNM3(1,2)=R0
      XNM3(1,3)=R0
      XNM3(2,1)=R0
      XNM3(2,2)=R1
      XNM3(2,3)=R0
      XNM3(3,1)=R0
      XNM3(3,2)=R0
      XNM3(3,3)=R1
C     Add first three terms of series
C     -------------------------------
      DO 70 I=1,NDIM
         DO 60 J=1,NDIM
            EXPX(I,J)=RP5*XNM1(I,J)+XNM2(I,J)+XNM3(I,J)
 60      CONTINUE
 70   CONTINUE
C     Add remaining terms (with X to the powers 3 to NMAX)
C     ----------------------------------------------------
      FACTOR=R2
      DO 140 N=3,NMAX
C     Use recursive formula to obtain X to the power N
         DO 90 I=1,NDIM
            DO 80 J=1,NDIM
               XN(I,J)=C1*XNM1(I,J)-C2*XNM2(I,J)+C3*XNM3(I,J)
 80         CONTINUE
 90      CONTINUE
C     Update factorial
         FACTOR=DBLE(N)*FACTOR
         R1DFAC=R1/FACTOR
C     Add Nth term of the series
         DO 110 I=1,NDIM
            DO 100 J=1,NDIM
               EXPX(I,J)=EXPX(I,J)+R1DFAC*XN(I,J)
 100        CONTINUE
 110     CONTINUE
C     Check convergence of series
         XNNORM=SQRT(SCAPRD(XN(1,1),XN(1,1),NDIM2))
         IF(XNNORM.GT.OVER.OR.(XNNORM.LT.UNDER.AND.XNNORM.GT.R0)
     1        .OR.R1DFAC.LT.UNDER)THEN
C...  first check possibility of overflow or underflow.
C...  numbers are to small or too big: Break (unconverged) loop and exit
            NOCONV=.TRUE.
            GOTO 999
         ELSEIF(XNNORM*R1DFAC.LT.TOL)THEN
C...  converged: Break series summation loop and exit with success
            GOTO 999
         ENDIF 
         DO 130 I=1,NDIM
            DO 120 J=1,NDIM
               XNM3(I,J)=XNM2(I,J)
               XNM2(I,J)=XNM1(I,J)
               XNM1(I,J)=XN(I,J)
 120        CONTINUE
 130     CONTINUE
 140  CONTINUE
C     Re-set convergence flag if series did not converge
      NOCONV=.TRUE.
C     
 999  CONTINUE
      RETURN
      END

 


      SUBROUTINE DEXPMP
     1(   DEXPX      ,NOCONV     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=3     ,NDIM2=9    ,NDIM4=81   ,MAXN=100   )
C Arguments
      LOGICAL  NOCONV
      DIMENSION
     1    DEXPX(NDIM,NDIM,NDIM,NDIM),X(NDIM,NDIM)
C Local arrays and variables
C...matrix of powers of X
      DIMENSION
     1    R1DFAC(MAXN)       ,XMATX(NDIM,NDIM,0:MAXN)
C...initialise identity matrix: X to the power 0
      DATA
     1    XMATX(1,1,0)  ,XMATX(1,2,0)  ,XMATX(1,3,0)  /
     2    1.D0          ,0.D0          ,0.D0          /
     3    XMATX(2,1,0)  ,XMATX(2,2,0)  ,XMATX(2,3,0)  /
     4    0.D0          ,1.D0          ,0.D0          /
     5    XMATX(3,1,0)  ,XMATX(3,2,0)  ,XMATX(3,3,0)  /
     6    0.D0          ,0.D0          ,1.D0          /
      DATA
     1    R0   ,RP5  ,R1   ,TOL    ,OVER    ,UNDER   /
     2    0.0D0,0.5D0,1.0D0,1.0D-10,1.0D+100,1.0D-100/
C***********************************************************************
C COMPUTES THE DERIVATIVE OF THE EXPONENTIAL OF A (GENERALLY
C UNSYMMETRIC) 3-D TENSOR. USES THE SERIES REPRESENTATION OF THE TENSOR
C EXPONENTIAL.
C
C REFERENCE: Section B.2
C            Box B.2
C***********************************************************************
C Initialise convergence flag
      NOCONV=.FALSE.
C X to the power 1
      DO 20 I=1,NDIM
        DO 10 J=1,NDIM
          XMATX(I,J,1)=X(I,J)
   10   CONTINUE
   20 CONTINUE
C Zero remaining powers of X
      CALL RVZERO(XMATX(1,1,2),NDIM*NDIM*(MAXN-1))
C Compute X square
      DO 50 I=1,NDIM
        DO 40 J=1,NDIM
          DO 30 K=1,NDIM
            XMATX(I,J,2)=XMATX(I,J,2)+X(I,K)*X(K,J)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
C Compute principal invariants of X
      C1=X(1,1)+X(2,2)+X(3,3)
      C2=RP5*(C1*C1-(XMATX(1,1,2)+XMATX(2,2,2)+XMATX(3,3,2)))
      C3=X(1,1)*X(2,2)*X(3,3)+X(1,2)*X(2,3)*X(3,1)+
     1   X(1,3)*X(2,1)*X(3,2)-X(1,2)*X(2,1)*X(3,3)-
     2   X(1,1)*X(2,3)*X(3,2)-X(1,3)*X(2,2)*X(3,1)
C Compute X to the powers 3,4,...,NMAX using recursive formula
      R1DFAC(1)=R1
      R1DFAC(2)=RP5
      DO 80 N=3,MAXN
        R1DFAC(N)=R1DFAC(N-1)/DBLE(N)
        DO 70 I=1,NDIM
          DO 60 J=1,NDIM
            XMATX(I,J,N)=C1*XMATX(I,J,N-1)-C2*XMATX(I,J,N-2)+
     1                   C3*XMATX(I,J,N-3)
   60     CONTINUE
   70   CONTINUE
        XNNORM=SQRT(SCAPRD(XMATX(1,1,N),XMATX(1,1,N),NDIM2))
C...check number of terms required for series convergence
        IF(XNNORM.GT.OVER.OR.(XNNORM.LT.UNDER.AND.XNNORM.GT.R0)
     1                                  .OR.R1DFAC(N).LT.UNDER)THEN
C...numbers are to small or too big: Exit without computing derivative
          NOCONV=.TRUE.
          GOTO 999
        ELSEIF(XNNORM*R1DFAC(N).LT.TOL)THEN
C...series will converge with NMAX terms:
C   Carry on to derivative computation
          NMAX=N
          GOTO 90
        ENDIF
   80 CONTINUE
C...series will not converge for the currently prescribed tolerance
C   with the currently prescribed maximum number of terms MAXN:
C   Exit without computing derivative
      NOCONV=.TRUE.
      GOTO 999
   90 CONTINUE
C Compute the derivative of exponential map
      CALL RVZERO(DEXPX,NDIM4)
      DO 150 I=1,NDIM
        DO 140 J=1,NDIM
          DO 130 K=1,NDIM
            DO 120 L=1,NDIM
              DO 110 N=1,NMAX
                DO 100 M=1,N
                  DEXPX(I,J,K,L)=DEXPX(I,J,K,L)+
     1                           R1DFAC(N)*XMATX(I,K,M-1)*XMATX(L,J,N-M)
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
C
  999 CONTINUE
      RETURN
      END


      

      SUBROUTINE RVZERO
     1     (   V          ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N)
      DATA R0/0.0D0/
C***********************************************************************
C     INITIALISES TO ZERO A DOUBLE PRECISION ARRAY OF DIMENSION N
C***********************************************************************
      DO 10 I=1,N
         V(I)=R0
 10   CONTINUE
      RETURN
      END


      DOUBLE PRECISION FUNCTION SCAPRD(U  ,V  ,N  ) 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N), V(N)
      DATA  R0  / 0.0D0 /
C***********************************************************************
C SCALAR PRODUCT OF DOUBLE PRECISION VECTORS U AND V OF DIMENSION N
C***********************************************************************
      SCAPRD=R0
      DO 10 I=1,N
        SCAPRD=SCAPRD+U(I)*V(I)
  10  CONTINUE
      RETURN
      END




      SUBROUTINE norm(dnorm,MAT,m,n)
      REAL*8 MAT(m,n),dnorm
      integer ii,jj,m,n

      dnorm=0d0
      DO,ii=1,m
         do,jj=1,n
            dnorm=dnorm+MAT(ii,jj)*MAT(ii,jj)
         enddo
      enddo
      dnorm=dsqrt(dnorm)
      RETURN 
      END

      subroutine TENS3333_sym(A,A66)
      real*8 A(3,3,3,3)
      real*8 a66(6,6)
      integer indexi(6),indexj(6),iflag
      integer indexii(3,3)
      integer ii,jj
  

      do,ii=1,3
         indexi(ii)=ii
         indexj(ii)=ii
      enddo
      
      indexi(4)=1
      indexj(4)=2
      indexi(5)=1
      indexj(5)=3
      indexi(6)=2
      indexj(6)=3
      indexii(1,1)=1
      indexii(2,2)=2
      indexii(3,3)=3
      indexii(1,2)=4
      indexii(1,3)=5
      indexii(2,3)=6
      indexii(2,1)=4
      indexii(3,1)=5
      indexii(3,2)=6
     
      do,ii=1,6
         do,jj=1,6
c            write(*,*) 'de ij a ijkl',ii,jj
c            write(*,*) '+,',indexi(ii),indexj(ii),indexi(jj),indexj(jj)
c            write(*,*) '+,',indexi(ii),indexj(ii),indexj(jj),indexi(jj)

c            write(*,*)  A(indexi(ii),indexj(ii),indexi(jj),indexj(jj))
c            a66(ii,jj) = .5d0*
c     1      (A( indexi(ii),indexj(ii),indexi(jj),indexj(jj))+
c     1       A( indexi(ii),indexj(ii),indexj(jj),indexi(jj)) )
            a66(ii,jj) = A( indexi(ii),indexj(ii),indexi(jj),indexj(jj))
c$$$            a66(ii,jj) = .25d0*
c$$$     1           (A( indexi(ii),indexj(ii),indexi(jj),indexj(jj))+
c$$$     1           A( indexi(ii),indexj(ii),indexj(jj),indexi(jj))+
c$$$     1           A( indexj(ii),indexi(ii),indexi(jj),indexj(jj))+
c$$$     1           A( indexj(ii),indexi(ii),indexj(jj),indexi(jj)))
         
         enddo
c         write(*,100) (a66(ii,j),j=1,6)
      enddo
      RETURN 
 100  FORMAT( 6(E16.8,' ')) 
      END
      


      subroutine TENS3333(A,A99)
      real*8 A(3,3,3,3)
      real*8 a99(9,9)
      integer indexi(9),indexj(9),iflag
      integer indexii(3,3)
      
       do,ii=1,3
         indexi(ii)=ii
         indexj(ii)=ii
      enddo
      
      indexi(4)=1
      indexj(4)=2
      indexi(5)=1
      indexj(5)=3
      indexi(6)=2
      indexj(6)=3
      indexi(7)=2
      indexj(7)=1
      indexi(8)=3
      indexj(8)=1
      indexi(9)=3
      indexj(9)=2
      
      indexii(1,1)=1
      indexii(2,2)=2
      indexii(3,3)=3
      indexii(1,2)=4
      indexii(2,1)=7
      indexii(1,3)=5
      indexii(3,1)=8
      indexii(2,3)=6
      indexii(3,2)=9
      do,ii=1,9
         do,jj=1,9
            a99(ii,jj)=A(indexi(ii),indexj(ii),indexi(jj),indexj(jj))
         enddo
      enddo
c$$$      do,ii=1,9
c$$$         write(*,100) (A99(ii,j), j=1,9)
c$$$      enddo
      return
 100  FORMAT( 9(E16.8,' ')) 
      end

      subroutine INVERSE3333_sym(A_I,A,iflag)
      real*8 A_I(3,3,3,3),A(3,3,3,3)
      real*8 a66(6,6)
      integer indexi(6),indexj(6),iflag
      integer indexii(3,3)
       do,ii=1,3
         indexi(ii)=ii
         indexj(ii)=ii
      enddo
      
      indexi(4)=1
      indexj(4)=2
      indexi(5)=1
      indexj(5)=3
      indexi(6)=2
      indexj(6)=3
      indexii(1,1)=1
      indexii(2,2)=2
      indexii(3,3)=3
      indexii(1,2)=4
      indexii(1,3)=5
      indexii(2,3)=6
      do,ii=1,6
         do,jj=ii,6
            a66(ii,jj) = .25d0*
     1           (A( indexi(ii),indexj(ii),indexi(jj),indexj(jj))+
     1           A( indexi(ii),indexj(ii),indexj(jj),indexi(jj))+
     1           A( indexj(ii),indexi(ii),indexi(jj),indexj(jj))+
     1           A( indexj(ii),indexi(ii),indexj(jj),indexi(jj)))
            if(ii.NE.jj) a66(jj,ii)=a66(ii,jj)
         enddo
         
      enddo

      call MATINV(A66,6,6,6,IFLAG) 

      if(iflag.NE.0)  write(*,*)'iflag',iflag
c$$$      do,ii=1,6
c$$$         write(*,100) (a66(ii,j),j=1,6)
c$$$      enddo
      do,ii=1,3
         do,jj=1,3
            do,kk=1,3
               do,ll=1,3
                  A_I(ii,jj,kk,ll)=A66(indexii(ii,jj),indexii(kk,ll))
                  A_I(jj,ii,kk,ll)=A66(indexii(ii,jj),indexii(kk,ll))
                  A_I(jj,ii,ll,kk)=A66(indexii(ii,jj),indexii(kk,ll))
                  A_I(ii,jj,ll,kk)=A66(indexii(ii,jj),indexii(kk,ll))
               enddo
            enddo
         enddo
      enddo
 100  FORMAT( 6(e16.8,' '))
      RETURN
      END

      subroutine INVERSE3333(A_I,A,iflag)
      real*8 A_I(3,3,3,3),A(3,3,3,3)
      real*8 a99(9,9)
      integer indexi(9),indexj(9),iflag
      integer indexii(3,3)
      do,ii=1,3
         indexi(ii)=ii
         indexj(ii)=ii
      enddo
      
      indexi(4)=1
      indexj(4)=2
      indexi(5)=1
      indexj(5)=3
      indexi(6)=2
      indexj(6)=3
      indexi(7)=2
      indexj(7)=1
      indexi(8)=3
      indexj(8)=1
      indexi(9)=3
      indexj(9)=2
      
      indexii(1,1)=1
      indexii(2,2)=2
      indexii(3,3)=3
      indexii(1,2)=4
      indexii(2,1)=7
      indexii(1,3)=5
      indexii(3,1)=8
      indexii(2,3)=6
      indexii(3,2)=9
      do,ii=1,9
         do,jj=1,9
            a99(ii,jj)=A(indexi(ii),indexj(ii),indexi(jj),indexj(jj))
c            write(*,*) 'de ij a ijkl',ii,jj,a99(ii,jj)
c            write(*,*) indexi(ii),indexj(ii),indexi(jj),indexj(jj)
         enddo
      enddo
c$$$      do,ii=1,9
c$$$         write(*,100) (A99(ii,j), j=1,9)
c$$$      enddo
c$$$      write(*,*)'**INV**'
      call MATINV(A99,9,9,9,IFLAG)
      if(iflag.NE.0)  write(*,*)'iflag',iflag
c$$$      do,ii=1,9
c$$$         write(*,100) (A99(ii,j), j=1,9)
c$$$      enddo
      do,ii=1,3
         do,jj=1,3
            do,kk=1,3
               do,ll=1,3
                  A_I(ii,jj,kk,ll)=A99(indexii(ii,jj),indexii(kk,ll))
               enddo
            enddo
         enddo
      enddo
 100  FORMAT( 9(e16.8,' '))
      RETURN
      END
      
      integer function kroneker(i,j)
      integer i,j
      kroneker=0
      if(i.EQ.j) kroneker=1
      return 
      end


      SUBROUTINE ROTATE_TENS3(AROT,A,R)
c     A given in system ei
C     ei*=Aei
C     Arot=R^T A R
      REAL*8 A(3,3),R(3,3),AROT(3,3)
      REAL*8 aux33(3,3)
      
      do,i=1,3
         do,j=1,3
            aux33(i,j)=0d0
            do,k=1,3
               aux33(i,j)=R(k,i)*A(k,j)+aux33(i,j)
            enddo
         enddo
      enddo
      do,i=1,3
         do,j=1,3
            arot(i,j)=0d0
            do,k=1,3
               arot(i,j)=aux33(i,k)*R(k,j)+arot(i,j)
            enddo
         enddo
      enddo
      return
      end

      Subroutine LUDCMP(A,N,INDX,D,CODE)
      PARAMETER(NMAX=100,TINY=1d-16)
      REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX)
      INTEGER CODE, D, INDX(N)

      D=1; CODE=0

      DO I=1,N
         AMAX=0.d0
         DO J=1,N
             IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
         END DO                 ! j loop
         IF(AMAX.LT.TINY) THEN
             CODE = 1
             RETURN
         END IF
         VV(I) = 1.d0 / AMAX
      END DO                    ! i loop
      
      DO J=1,N
         DO I=1,J-1
            SUM = A(I,J)
            DO K=1,I-1
               SUM = SUM - A(I,K)*A(K,J) 
                   END DO       ! k loop
                   A(I,J) = SUM
                END DO          ! i loop
                AMAX = 0.d0
                DO I=J,N
                   SUM = A(I,J)
                   DO K=1,J-1
                      SUM = SUM - A(I,K)*A(K,J) 
                   END DO       ! k loop
                   A(I,J) = SUM
                   DUM = VV(I)*DABS(SUM)
                   IF(DUM.GE.AMAX) THEN
                      IMAX = I
                      AMAX = DUM
                   END IF
                END DO          ! i loop  
   
                IF(J.NE.IMAX) THEN
                    DO K=1,N
                       DUM = A(IMAX,K)
                       A(IMAX,K) = A(J,K)
                       A(J,K) = DUM
                    END DO      ! k loop
                    D  = -D
                    VV(IMAX) = VV(J)
                 END IF

                 INDX(J) = IMAX
                 IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

                 IF(J.NE.N) THEN
                    DUM = 1.d0 / A(J,J)
                        DO I=J+1,N
                           A(I,J) = A(I,J)*DUM
                        END DO  ! i loop
                     END IF 
                  END DO        ! j loop

                  RETURN
                  END subroutine LUDCMP

c
      SUBROUTINE gauss_3(AA,BB,xx,nmax,n,det,iflag)
      integer iflag
      real*8 AA(nmax,nmax),bb(nmax),xx(nmax)
      real*8 A(n,n),b(n),det
      integer n,nmax,i,j, d, indx(n)
      det=0.
     
      do,i=1,n
         b(i)=bb(i)
         do,j=1,n
            A(i,j)=AA(i,j)
         enddo
      enddo
      call ludcmp(a,n,indx,d,iflag)
      call LUBKSB(A,N,INDX,B)
      do,i=1,n
         xx(i)=b(i)
      enddo
      return
      end

      Subroutine LUBKSB(A,N,INDX,B)
      REAL*8  SUM, A(N,N),B(N)
      INTEGER INDX(N)
      
      II = 0
      

      DO I=1,N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF(II.NE.0) THEN
            DO J=II,I-1
               SUM = SUM - A(I,J)*B(J)
            END DO              ! j loop
         ELSE IF(SUM.NE.0.d0) THEN
            II = I
         END IF
         B(I) = SUM
      END DO                    ! i loop
      
      DO I=N,1,-1
         SUM = B(I)
         IF(I < N) THEN
            DO J=I+1,N
               SUM = SUM - A(I,J)*B(J)
            END DO              ! j loop
         END IF
         B(I) = SUM / A(I,I)
      END DO                    ! i loop
      
     

      RETURN
      END subroutine LUBKSB


! end of file lu.f90




!      subroutine determinant(A,n,d)
!      REAL*8 A(n,n),b(n,n)
!      integer indx(n)
!      real*8 d
!      do,i=1,n
!         do,j=1,n
!            b(i,j)=a(i,j)
!         enddo
!      enddo
!      CALL ludcmp(b,n,n,indx,d)
!      do  j=1,n
!         d=d*b(j,j)
!      enddo
!      return
!      end

      subroutine dzeros_vec(A,i)
      real*8 A(i)
      integer ii,i
      do,ii=1,i
         A(ii)=0d0
      enddo
      return 
      end

      subroutine dzeros(A,i,j)
      real*8 A(i,j)
      integer i,j,ii,jj
      do,ii=1,i
         do,jj=1,j
            A(ii,jj)=0d0
         enddo
      enddo
      return
      end

      subroutine dzeros_tensor4(A,i,j,k,l)
      real*8 A(i,j,k,l)
      integer i,j,ii,jj
      do,ii=1,i
         do,jj=1,j
            do,kk=1,k
               do,ll=1,l
                  A(ii,jj,kk,ll)=0d0
               enddo
            enddo
         enddo
      enddo
      return

      end

c     Residual Subroutine
      subroutine residual(RES_CURR,DFG_ELAS_CURR,DFGRD_ELAS_pred0,
     1 tau_pred,schmidt,STIFF_ORIENT,nsystems,
     1 nsystems_max,gamma_0,mm,dtime,gamma_dot)
      
      implicit none
      integer ii,jj,isystem,nsystems,kroneker,nsystems_max
      real*8 DFG_ELAS_CURR(3,3),DFGRD_ELAS_pred0(3,3),
     1 L_p(3,3),EXP_L_P(3,3),DFGRD_ELAS_1(3,3),
     1 UU(3,3),RR(3,3),CC(3,3),EPSILON_EL(3,3),
     1 stress_k(3,3),AUX33_2(3,3),RES_CURR(3,3),
     1 temp_lp(3,3)
      real*8 STIFF_ORIENT(3,3,3,3)
      real*8 schmidt(nsystems_max,3,3)
      real*8 tau_resolved(nsystems_max),gamma_dot(nsystems_max),
     1 tau_pred(nsystems_max)
      real*8 gamma_0,mm,dtime
      logical noconv
      call POLAR_DECOMP(UU,RR,CC,DFG_ELAS_CURR)
      call LAGRANGIAN_DEF(CC,EPSILON_EL)
      call PCONTRACT2(stress_k,STIFF_ORIENT,epsilon_el,3)
      do,isystem=1,nsystems
         do,ii=1,3
            do,jj=1,3
               AUX33_2(ii,jj)=schmidt(isystem,ii,jj)
            enddo
         enddo
         call PCONTRACT(tau_resolved(isystem),stress_k,AUX33_2,3)
      enddo
      do,isystem=1,nsystems
         call VISCOLAW(gamma_dot(isystem),gamma_0,mm,
     1        tau_resolved(isystem),tau_pred(isystem),noconv)
         if(noconv .EQV. .TRUE.) THEN
         write(*,*), 'Convergency Error in Viscolaw'
         endif
      enddo
      call dzeros(L_p,3,3)
      do,isystem=1,nsystems
         do,ii=1,3
            do,jj=1,3
               L_p(ii,jj)=L_p(ii,jj)+gamma_dot(isystem)* 
     1              schmidt(isystem,ii,jj)
            enddo
         enddo
      enddo
      do,ii=1,3
         do,jj=1,3
            L_P(ii,jj)=L_P(ii,jj)*(-1d0*dtime)
         enddo
      enddo

      do,ii=1,3
         do,jj=1,3
            EXP_L_P(ii,jj)=kroneker(ii,jj)+L_p(ii,jj)
         enddo

      enddo
      call PMAT(DFGRD_ELAS_1,DFGRD_ELAS_PRED0,exp_l_p,3,3,3)

       do,ii=1,3
         do,jj=1,3
          RES_CURR(ii,jj)=DFG_ELAS_CURR(ii,jj)-DFGRD_ELAS_1(ii,jj)
         enddo
       enddo
      return
      end
c Hossein
      subroutine min_location(min_value,min_loc,input_array,n)
      implicit none

      integer ii, n, min_loc
      real*8 min_value, input_array(n)

      min_value=input_array(1)
      min_loc=1

      do, ii=2,n
       if(min_value.gt.input_array(ii)) then
       min_value=input_array(ii)
       min_loc=ii
       endif
      enddo

      return
      end

           
