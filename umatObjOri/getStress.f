subroutine ABQMAIN
!====================================================================
! This program must be compiled and linked with the command:
!     abaqus make job=fprin
! Run the program using the command:
!     abaqus fprin
!====================================================================
!
!  Purpose: 
! 
!    This program computes the principal stresses and strains and their
!    directions from stress and strain values stored in an ABAQUS 
!    results file (.fil).
!
!  Input File names: `FNAME.fil', where FNAME is the root file name of 
!                               the input file.
!
!  Output File name:  pvalue.dat
!
!====================================================================
!
!  Variables used by this program and ABAQUS subroutine SPRIND :
!
!   NDI    -- Number of direct components in stress/strain tensor.
!   NSHR   -- Number of shear components in stress/strain tensor.
!   NDIP1  -- NDI + 1
!   ARRAY  -- Real array containing values read from results file 
!              (.fil). Equivalenced to JRRAY.
!   JRRAY  -- Integer array containing values read from results file 
!              (.fil). Equivalenced to ARRAY.
!   FNAME  -- Root file name of input file (w/o .fil extension).
!   NRU    -- Number of results files (.fil) to be read.
!   LRUNIT -- Array containing unit number and format of results files:
!               LRUNIT(1,*) --> Unit number of input file.
!               LRUNIT(2,*) --> Format of input file.
!   LOUTF  -- Format of output file:
!               0 --> Standard ASCII format.
!               1 --> ABAQUS results file ASCII format.
!               2 --> ABAQUS results file binary format.
!   JUNIT  -- Unit number of file to be opened.
!   JRCD   -- Error check return code.
!               .EQ. 0 --> No errors.
!               .NE. 0 --> Errors detected.
!   KEY    -- Current record key identifier.
!   JELNUM -- Current element number.
!   INTPN  -- Integration point number.
!   LSTR   -- Indicates type of principal value (stress/strain) and
!               ordering used:
!                 For calculation of principal value (stress/strain):
!                   1 -->  stress.
!                   2 -->  strain. 
!                 For calculation of directions:
!                   1 -->  stress.
!                   2 -->  strain.
!   S      -- Array containing stress tensor.
!   PS     -- Array containing principal stresses.
!   ANPS   -- Array containing directions of principal stresses.
!   E      -- Array containing strain tensor.
!   PE     -- Array containing principal strains.
!   ANPE   -- Array containing directions of principal strains.
!
!====================================================================
!
!  The use of ABA_PARAM.INC eliminates the need to have different  
!  versions of the code for single and double precision.  
!  ABA_PARAM.INC defines an appropriate IMPLICIT REAL statement 
!  and sets the value of NPRECD to 1 or 2, depending on whether 
!  the machine uses single or double precision.
!
!====================================================================
!
!      INCLUDE 'aba_param.inc'
      
      real(kind = 8)    :: array(513)
      real(kind = 8)    :: S(6), E(6), Ps(3), PE(3), anps(3, 3), anpe(3, 3)
      integer(kind = 4) :: NRU, LOUTF, LRUNIT(2, 1), KEY, JRCD, JUNIT
      integer(kind = 4) :: jElnum, intPn, locate
      integer(kind = 4) :: K100, K1, IXX, IYY, IZZ, NDI, nshr, NDIP1
      character(len=80) :: fname

!====================================================================
!  Get the name of the results file.
!
!====================================================================
      WRITE(6,*) 'Enter the name of the input file (w/o .fil):'
      READ(5,'(A)') FNAME

!====================================================================
!  Open the output file.
!
!====================================================================
      OPEN(UNIT=9,FILE='pvalue.dat',STATUS='NEW')

      NRU = 1
      LOUTF = 0
      LRUNIT(1,1) = 8
      LRUNIT(2,1) = 2

      CALL INITPF(FNAME,NRU,LRUNIT,LOUTF)

      JUNIT = 8

      CALL DBRNU(JUNIT)

!====================================================================
!  Read records from the results (.fil) file and process the data.  
!  Cover a maximum of 10 million records in the file. 
!
!====================================================================
      DO 1000 K100 = 1, 100
      DO 1000 K1 = 1, 99999
         CALL DBFILE(0,ARRAY,JRCD)
         IF (JRCD .NE. 0) GO TO 1001
!         KEY = JRRAY(1,2)
         key = nint(array(2))
!
!====================================================================
!  Get the heading (title) record.
!
!====================================================================
         IF (KEY .EQ. 1922) THEN
            WRITE(9,1100) (ARRAY(IXX),IXX=3,12)
 1100       FORMAT(1X,10A8)

!====================================================================
!  Get the current step and increment number.
!
!====================================================================
         ELSE IF (KEY .EQ. 2000) THEN
!            WRITE(9,1200) JRRAY(1,8), JRRAY(1,9)
            WRITE(9,1200) nint(array(8)), nint(array(9))
 1200       FORMAT(1X,'** STEP ',I2,'    INCREMENT ',I3)

!====================================================================
!  Get the element and integration point numbers, JELNUM and INTPN, 
!  and the location of INTPN (0--at int.pt., 1--at centroid, 
!  4--nodal average) and the number of direct and shear components 
!  in the analysis.
!
!====================================================================
         ELSE IF (KEY .EQ. 1) THEN
          !  JELNUM = JRRAY(1,3)
          !  INTPN  = JRRAY(1,4)
          !  LOCATE = JRRAY(1,6)
          !  NDI    = JRRAY(1,8)
          !  NSHR   = JRRAY(1,9)
            JELNUM = nint(array(3))
            INTPN  = nint(array(4))
            LOCATE = nint(array(6))
            NDI    = nint(array(8))
            NSHR   = nint(array(9))
            NDIP1  = NDI + 1
            IF(LOCATE.LE.1) THEN
              WRITE(9,1201) JELNUM, INTPN ,NDI,NSHR
 1201         FORMAT(2X,'ELEMENT NUMBER = ',I8,5X, 'INT. PT. NUMBER = ',I2,5X, 'NDI/HSHR = ',2I2)
            ELSEIF(LOCATE.EQ.4) THEN
              WRITE(9,1191) JELNUM, NDI,NSHR
 1191         FORMAT(2X,'NODE NUMBER = ',I8,5X,'NDI/HSHR = ',2I2)
            END IF

!====================================================================
!  Get the stress tensor. 
!
!====================================================================
         ELSE  IF (KEY .EQ. 11) THEN
            WRITE(9,1202)
 1202       FORMAT(3X,'STRESSES:')

            DO 10 IXX = 1, NDI
               S(IXX) = ARRAY(IXX+2)
   10       CONTINUE
            WRITE(9,1203) (S(IZZ), IZZ = 1, NDI)
 1203       FORMAT(4X,'S11 = ',E12.5,'  S22 = ',E12.5,'  S33 = ',E12.5)
            DO 20 IYY = NDI + 1, NSHR + NDI
               S(IYY) = ARRAY(IYY+2)             
   20       CONTINUE
            WRITE(9,1204) (S(IZZ), IZZ = NDI + 1, NSHR + NDI)
 1204       FORMAT(4X,'S12 = ',E12.5,'  S13 = ',E12.5,'  S23 = ',E12.5)


!====================================================================
!  Calculate the principal stresses and corresponding principal 
!  directions in unsorted order.
!====================================================================
            LSTR = 1
            CALL SPRIND(S,PS,ANPS,LSTR,NDI,NSHR)
            WRITE(9,1205) PS(1), ANPS(1,1), ANPS(1,2), ANPS(1,3)
 1205       FORMAT(4X,'PS1 = ',E12.5,/, 5X,'PD11 =',F8.3,2X,'PD12 =',F8.3,2X,'PD13 =',F8.3)
            WRITE(9,1206) PS(2), ANPS(2,1), ANPS(2,2), ANPS(2,3)
 1206       FORMAT(4X,'PS2 = ',E12.5,/, 5X,'PD21 =',F8.3,2X,'PD22 =',F8.3,2X,'PD23 =',F8.3)
            WRITE(9,1207) PS(3), ANPS(3,1), ANPS(3,2), ANPS(3,3)
 1207       FORMAT(4X,'PS3 = ',E12.5,/, 5X,'PD31 =',F8.3,2X,'PD32 =',F8.3,2X,'PD33 =',F8.3)


!====================================================================
!  Get the strain tensor. 
!
!====================================================================
         ELSE  IF (KEY .EQ. 21) THEN
            WRITE(9,2202)
 2202       FORMAT(3X,'STRAINS:')

            DO 30 IXX = 1, NDI
               E(IXX) = ARRAY(IXX+2)
   30       CONTINUE
            WRITE(9,2203) (E(IZZ), IZZ = 1, NDI)
 2203       FORMAT(4X,'E11 = ',E12.5,'  E22 = ',E12.5,'  E33 = ',E12.5)
            DO 40 IYY = NDI + 1, NSHR + NDI
               E(IYY) = ARRAY(IYY+2)             
   40       CONTINUE
            WRITE(9,2204) (E(IZZ), IZZ = NDI + 1, NSHR + NDI)
 2204       FORMAT(4X,'E12 = ',E12.5,'  E13 = ',E12.5,'  E23 = ',E12.5)


!====================================================================
!  Calculate the principal strains and corresponding principal 
!  directions in unsorted order.
!====================================================================
            LSTR = 2
            CALL SPRIND(E,PE,ANPE,LSTR,NDI,NSHR)
            WRITE(9,2205) PE(1), ANPE(1,1), ANPE(1,2), ANPE(1,3)
 2205       FORMAT(4X,'PE1 = ',E12.5,/, 5X,'PD11 =',F8.3,2X,'PD12 =',F8.3,2X,'PD13 =',F8.3)
            WRITE(9,2206) PE(2), ANPE(2,1), ANPE(2,2), ANPE(2,3)
 2206       FORMAT(4X,'PE2 = ',E12.5,/, 5X,'PD21 =',F8.3,2X,'PD22 =',F8.3,2X,'PD23 =',F8.3)
            WRITE(9,2207) PE(3), ANPE(3,1), ANPE(3,2), ANPE(3,3)
 2207       FORMAT(4X,'PE3 = ',E12.5,/, 5X,'PD31 =',F8.3,2X,'PD32 =',F8.3,2X,'PD33 =',F8.3)

         END IF

 1000 CONTINUE
 1001 CONTINUE

      CLOSE (UNIT=9)

      RETURN
      END
