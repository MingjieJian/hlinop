C
C   Table of contents:
C
C     SUBROUTINE HLINPROF(WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,
C    *                   H1FRC,HE1FRC,DOP,PROF)
C     SUBROUTINE HLINPROF_SET_FILES(STEFILE,SELFFILE,VCSFILE)
C     INTEGER FUNCTION HLINPROF_SOURCE()
C     SUBROUTINE HLINPROF_SET_MESSAGES(MODE)
C     SUBROUTINE HLINOP_SET_TB09(TBFILE,MODE,DLO,DHI)
C     BLOCK DATA HLPRFD
C     REAL*8 FUNCTION HLYAPRF(HFILE,WAVE,WAVE0,TEMP,XXNE,H1FRC,HE1FRC,...)
C     REAL*8 FUNCTION HLPCONV(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,...)
C     SUBROUTINE STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
C     SUBROUTINE STEHLS(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN,HLINM,IBOTH)
C     SUBROUTINE STEHL2(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN,HLINM,IBOTH)
C     SUBROUTINE HTB09P(LFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,PROF,VALID)
C     REAL*8 FUNCTION STHPCH(XL,K,NA,AA,SS), STHPD(J,NA,AA,SS)
C     SUBROUTINE HGRID(HFILE,WAVE,WAVE0,TEMP,NHYD,NLOW,NUP,HLIN)
C     SUBROUTINE HSELF_PDWIDTH(NLOW,NUP,TEMP,H1FRC,WIDTH)
C     SUBROUTINE HE_PDWIDTH(NLOW,NUP,TEMP,HE1FRC,HEWID)
C     SUBROUTINE RAD_WIDTH(NLOW,NUP,WIDTH)
C     SUBROUTINE ALI_GRIEM(NLOW,NUP,WAVE0,H1FRC,WIDTH)
C     SUBROUTINE VCS(PR,XNE,T,DEL,II,N,M)                                      
C
C
C  Hydrogen Line Routines
C  All routines by Paul Barklem unless otherwise stated
C
      SUBROUTINE HLINPROF(WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,
     *                   H1FRC,HE1FRC,DOP,PROF)
C
C  Returns H line Profile at specific point
C
C  Inputs: 
C          WAVE   - wavelength in AA
C          WAVE0  - line centre in AA
C          TEMP   - temperature in K
C          XXNE   - electron number density in cm^-3
C          NLOW   - n quantum number of lower level
C          NUP    - n quantum number of upper level
C          H1FRC  - neutral hydrogen number density in cm^-3
C          HE1FRC - neutral helium number density in cm^-3
C          DOP    - the doppler width delta lambda / lambda
C                   (only used if HLINOP is called instead!)
C  Outputs:
C          PROF   - area normalised (in AA) profile for WAVE
C  Common Inputs: 
C     (CURRENTLY SET HARD!!!! see below)
C          ISTARK - switch for different recipes
C                  0=none
C                  1=Stehle
C                  2=Vidal, Cooper & Smith
C          ISELF  - switch for self broadening
C                  0=none
C                  1=Barklem, Piskunov and O'Mara (grid)
C                  2=Barklem, Piskunov and O'Mara (p-d)
C                  3=Ali & Griem
C          ICONV  - switch for convolution method
C                  1=adding
C                  2=trapezoidal
C          IHE    - switch for broadening by He
C                  0=none
C                  1=estimated from Barklem et al (1997)
C          NB: ICONV=2 is incompatible with either 0 option 
C  
C  The trapezoidal convolution is only for testing.  Use at own risk!
C  It's left here in as it might be useful to someone.
C      
      IMPLICIT NONE

      INTEGER NKMAX
      PARAMETER (NKMAX = 10001)
      INTEGER NLOW,NUP,LINE
      INTEGER ISTARK,ISELF,ICONV
      INTEGER IHE
      INTEGER NLOWW,NHIGH,NKRN,IKRN,J,NMID,I1,I2
      REAL*4 TEMP,XXNE,PROF,H1FRC,HE1FRC
      REAL*4 HLIN,STARK
      REAL*4 HEWID,RAD,DOP
      REAL*4 HLINOP
      REAL*8 TD,XNED,H1FRCD
      REAL*8 WAVE,WAVE0,DEL,HHLIN
      REAL*8 KRNSTEP,WLKRN(NKMAX),KERN(NKMAX),KMAX,WV2
      REAL*8 DEL1,DDEL,XNEXP,F1,HLIND
      REAL*8 STARKD,PRA(1),DELA(1)
      REAL*8 LOR,HWID,RES
      REAL*8 GAML,AEXT,DOPTH,HLPCONV,HLYAPRF
      INTEGER CONVOLUTION_MODE
      COMMON /HLCONV/ CONVOLUTION_MODE
      EXTERNAL HLCONVD
      REAL*8 C,PI
      CHARACTER*256 HFILE,HVCSFILE,HSELFFILE
      LOGICAL FIRST
C  table files in use: from DATA.FILES, or set with HLINPROF_SET_FILES
      CHARACTER*592 HF,HSF,HVF
      CHARACTER*592 FSTE,FSELF,FVCS
      COMMON /HLFILN/ FSTE,FSELF,FVCS
      INTEGER NFSET,ISRC,LFSET
      COMMON /HLFILC/ NFSET,ISRC
      INTEGER IMSG
      COMMON /HLMSG/ IMSG
      EXTERNAL HLPRFD
C
      SAVE FIRST,HF,HSF,HVF,LFSET
C
      DATA FIRST/.TRUE./, LFSET/-1/
      PARAMETER (C = 2.997925E+18, PI = 3.14159265)
      INCLUDE 'DATA.FILES'
C
C  Set the switches here 
C
      COMMON /HSWITCH/ISTARK,ISELF,ICONV,IHE
      ISTARK = 1
      ISELF = 2
      ICONV = 1
      IHE = 1
C
C Table file names (only when they change)
C
      IF (LFSET.NE.NFSET) THEN
        HF = HFILE
        HSF = HSELFFILE
        HVF = HVCSFILE
        IF (FSTE.NE.' ') HF = FSTE
        IF (FSELF.NE.' ') HSF = FSELF
        IF (FVCS.NE.' ') HVF = FVCS
        LFSET = NFSET
      ENDIF
C
C We check we have a line which we can handle
C  if not then pass to HLINOP
C
C  Lines in the Stehle & Hutcheon (1999) table: Lyman, Balmer and 
C  Paschen to n_up = 30.  Lyman alpha uses the table for its core and
C  the Allard et al. quasi-molecular profiles for its far red wing 
C  (HLYAPRF) with the default settings; otherwise it is left to HLINOP,
C  which has the earlier treatment.
C
      IF ((NLOW .GT. 3) .OR. (NUP .GT. 30) .OR. 
     ;    ((NLOW .EQ. 1) .AND. (NUP .EQ. 2) .AND. 
     ;     .NOT.((CONVOLUTION_MODE.EQ.1).AND.(ICONV.EQ.1).AND.
     ;           (ISTARK.EQ.1).AND.(ISELF.NE.1)))) THEN
        IF (FIRST) THEN
          FIRST = .FALSE.
          IF (IMSG.EQ.1) PRINT*,' Using Hlinop code!'
        ENDIF
        ISRC = 0
        PROF = HLINOP(WAVE,NLOW,NUP,WAVE0,TEMP,XXNE,H1FRC,
     ;                                          HE1FRC,DOP)
        PROF = PROF * C/WAVE/WAVE
      RETURN
      END IF 
C
C For historical reasons copy some variables to real*8
C
      XNED = XXNE
      TD = TEMP
      H1FRCD = H1FRC
C      
C Precomputed things here
C
      DEL = DABS(WAVE0 - WAVE)
      LINE = NUP - NLOW
C
C Nov 2026: by default (HLINOP_SET_CONVOLUTION(1)) the tabulated 
C profile (Stark convolved with the thermal Doppler profile of H) is 
C convolved with the remaining Doppler broadening (DOP beyond the 
C thermal width, e.g. microturbulence; Gaussian of 1/e half width 
C WAVE0*sqrt(DOP**2-DOPTH**2)) and with the Lorentz profile of the 
C radiative, self (ISELF = 2 or 3) and helium widths (HLPCONV).  The
C additions below (Lorentzian added beyond 0.2 A, no microturbulence)
C are used with HLINOP_SET_CONVOLUTION(0), and for ISELF = 1 (self-
C broadening grid profiles) or ICONV = 2.
C
      IF ((CONVOLUTION_MODE.EQ.1).AND.(ICONV.EQ.1).AND.(ISTARK.EQ.1)
     ;    .AND.(ISELF.NE.1)) THEN
        CALL STEHLE(HF,WAVE0,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
        IF (HLIN.LT.0.) THEN
          ISRC = 0
          PROF = HLINOP(WAVE,NLOW,NUP,WAVE0,TEMP,XXNE,H1FRC,
     ;                                          HE1FRC,DOP)
          PROF = PROF * C/WAVE/WAVE
          RETURN
        ENDIF
        ISRC = 1
C  Lorentz half width (A): self + radiative + helium
        HWID = 0.D0
        IF (H1FRC.GT.0.) THEN
          IF (ISELF.EQ.2) THEN
            CALL HSELF_PDWIDTH(NLOW,NUP,TEMP,H1FRC,HWID)
            HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
          ELSE IF (ISELF.EQ.3) THEN
            CALL ALI_GRIEM(NLOW,NUP,WAVE0,H1FRC,RES)
            HWID = RES
          ENDIF
        ENDIF
        HEWID = 0.
        IF (IHE .EQ. 1 ) CALL HE_PDWIDTH(NLOW,NUP,TEMP,HE1FRC,HEWID)
        RAD = 0.
        CALL RAD_WIDTH(NLOW,NUP,RAD)
        GAML = HWID + (RAD + HEWID) * WAVE0 * WAVE0 / C / 2. / PI
C  Doppler broadening beyond the thermal width of H in the tables
        DOPTH = SQRT(2.D0*1.380649D-16*TEMP/1.6735575D-24)/2.99792458D10
        AEXT = 0.D0
        IF (DOP.GT.DOPTH) AEXT = WAVE0*SQRT(DBLE(DOP)**2-DOPTH**2)
        IF ((NLOW.EQ.1).AND.(NUP.EQ.2)) THEN
          PROF = HLYAPRF(HF,WAVE,WAVE0,TEMP,XXNE,H1FRC,HE1FRC,AEXT,
     ;                   DBLE(RAD)*WAVE0*WAVE0/C/2./PI)
          RETURN
        ENDIF
        PROF = HLPCONV(HF,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,AEXT,GAML)
        RETURN
      ENDIF
C
C Adding for convolution of Stark and self
C
      IF (ICONV .EQ. 1) THEN
        HLIN = 0. 
        HHLIN = 0.
        IF (ISTARK .EQ. 1) THEN
          CALL STEHLE(HF,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
          ISRC = 1
C
C If outside the table (temperature, or density below the table) use
C HLINOP, which includes all broadening mechanisms itself
C
          IF (HLIN.LT.0.) THEN
            ISRC = 0
            PROF = HLINOP(WAVE,NLOW,NUP,WAVE0,TEMP,XXNE,H1FRC,
     ;                                          HE1FRC,DOP)
            PROF = PROF * C/WAVE/WAVE
            RETURN
          ENDIF
        END IF
        IF (ISTARK .EQ. 2) THEN
C  (VCS takes arrays of detunings and profile values)
          DELA(1) = DEL
          CALL VCS(PRA,XNED,TD,DELA,1,NLOW,NUP,HVF)
          HLIND = PRA(1)
          ISRC = 2
          IF (HLIND.LT.0.) THEN
            ISRC = 0
            PROF = HLINOP(WAVE,NLOW,NUP,WAVE0,TEMP,XXNE,H1FRC,
     ;                                          HE1FRC,DOP)
            PROF = PROF * C/WAVE/WAVE
            RETURN
          ENDIF
          HLIN = HLIND
        END IF
        IF ((DEL .GT. .2d0) .AND. (H1FRC .GT. 0.)) THEN
          IF (ISELF .EQ. 1) THEN
             CALL HGRID(HSF,
     ;            WAVE,WAVE0,TD,H1FRCD,NLOW,NUP,HHLIN)
C
C If outside table range use p-d approx
C
            IF (HHLIN.LT.0.) THEN
              CALL HSELF_PDWIDTH(NLOW,NUP,TEMP,H1FRC,HWID)
              HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
              HHLIN = HWID / (HWID * HWID + DEL * DEL) / PI
            ENDIF
          ENDIF
          IF (ISELF .EQ. 2) THEN
            CALL HSELF_PDWIDTH(NLOW,NUP,TEMP,H1FRC,HWID)
            HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
            HHLIN = HWID / (HWID * HWID + DEL * DEL) / PI
          END IF
          IF (ISELF .EQ. 3) THEN
            CALL ALI_GRIEM(NLOW,NUP,WAVE0,H1FRC,RES)
            HHLIN = RES / (RES * RES + DEL * DEL) / PI
          END IF 
          HLIN = HLIN + HHLIN
        END IF
      END IF
C
C Convolution by simple trapezoidal method
C
C   Notes:  Should have Stark profile > Lorentzian
C   If ever this is broken you should reverse the two or 
C   increase NKRN points significantly.
C
      IF (ICONV.EQ.2) THEN
        IF (H1FRC .EQ. 0.) THEN
          PRINT*,'STOPPING : ZERO H1FRC, CONVOLVING WILL FAIL'
          STOP
        END IF
C
C Use a logarithmic spread of kernel points 
C
        KMAX = 100.     !we go out KMAX AA into wings
        NKRN = 1001     !use an odd number, testing may be 
C                       !needed to see what this number should 
C                       !for a given case
        NMID = (NKRN + 1) / 2
        NLOWW = NMID - 1
        NHIGH = NMID + 1
        XNEXP = DLOG(1.d-4)   !this is the low limit of detuning 
C                             !(the log of it in AA)
C                             !should be smaller than smallest half width
C                             !which are typically ~0.01AA for solar
        KRNSTEP = (LOG(KMAX) - XNEXP) / DBLE(NLOWW - 1)
        DO J = NKRN, NHIGH, -1
          WLKRN(J) = EXP(LOG(KMAX) - KRNSTEP * DBLE(NKRN - J))
        END DO
        WLKRN(NMID) = 0.d0
        DO J = 1, NLOWW, 1
          WLKRN(J) = - EXP(LOG(KMAX) - KRNSTEP * DBLE(J - 1))
        END DO
C
        DO IKRN = 1, NKRN
          WV2 = WAVE0 + WLKRN(IKRN)
            HHLIN = 0.
            IF (ISELF .EQ. 1) CALL HGRID(HSF,
     ;              WV2,WAVE0,TD,H1FRCD,NLOW,NUP,HHLIN)
            IF (ISELF .EQ. 2) THEN
              CALL HSELF_PDWIDTH(NLOW,NUP,TEMP,H1FRC,HWID)
              HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
              HHLIN = HWID / (HWID * HWID + WLKRN(IKRN) * WLKRN(IKRN)) 
     ,                     / PI
            END IF
            IF (ISELF .EQ. 3) THEN
              CALL ALI_GRIEM(NLOW,NUP,WAVE0,H1FRC,RES)
              HHLIN = RES / (RES * RES + WLKRN(IKRN) * WLKRN(IKRN)) / PI
            END IF 
          KERN(IKRN) = HHLIN
        END DO
C
C  Do Convolution -- Integration is simple Trapezoidal here
C
        HLIND = 0.d0
        DO IKRN = 1, NKRN
          DEL1 = DEL + WLKRN(IKRN)
          I1 = MAX0(   1, IKRN - 1)
          I2 = MIN0(NKRN, IKRN + 1)
		  DDEL = WLKRN(I2) - WLKRN(I1)
          IF (ISTARK .EQ. 1) THEN
            CALL STEHLE(HF,WAVE0+DEL1,WAVE0,TEMP,XXNE,NLOW,NUP,STARK)
          END IF
          IF (ISTARK .EQ. 2) THEN
            DELA(1) = DEL1
            CALL VCS(PRA,XNED,TD,DELA,1,NLOW,NUP,HVF)
            STARKD = PRA(1)
            STARK = STARKD
          END IF
          F1 = KERN(IKRN) * STARK
          HLIND = HLIND + DDEL * F1
        END DO
        HLIND= HLIND * .5d0
        HLIN = HLIND
        ISRC = ISTARK
      END IF
C
C The He and radiative broadening are always assumed small enough to 
C convolve by adding 
C
      HEWID = 0.
      IF (IHE .EQ. 1 ) CALL HE_PDWIDTH(NLOW,NUP,TEMP,HE1FRC,HEWID)
      CALL RAD_WIDTH(NLOW,NUP,RAD)
C
C  Compute the width in lambda units (AA)
C
      LOR = RAD + HEWID
      LOR = LOR * WAVE0 * WAVE0 / C / 2. / PI
C
      IF (DEL .GT. .2)
     ;  HLIN = HLIN + LOR / (LOR * LOR + DEL * DEL) / PI
C
      PROF = HLIN
      RETURN
C
      END

C***********************************************************************
      SUBROUTINE HLINPROF_SET_FILES(STEFILE,SELFFILE,VCSFILE)
C
C  Table files for HLINPROF at run time (full paths, at most 592 
C  characters): the Stehle & Hutcheon table (stehle_long.dat), the 
C  self-broadening grid (ISELF = 1) and the VCS table (ISTARK = 2).  A
C  blank name keeps the file compiled in from DATA.FILES.  Call before 
C  the first profile: each table is read once.
C
      CHARACTER*(*) STEFILE,SELFFILE,VCSFILE
      CHARACTER*592 FSTE,FSELF,FVCS
      COMMON /HLFILN/ FSTE,FSELF,FVCS
      INTEGER NFSET,ISRC
      COMMON /HLFILC/ NFSET,ISRC
      EXTERNAL HLPRFD
      IF (LEN(STEFILE).GT.592 .OR. LEN(SELFFILE).GT.592 .OR.
     *    LEN(VCSFILE).GT.592)
     *  STOP 'HLINPROF_SET_FILES: file name longer than 592 characters'
      FSTE = STEFILE
      FSELF = SELFFILE
      FVCS = VCSFILE
      NFSET = NFSET + 1
      RETURN
      END

C***********************************************************************
      INTEGER FUNCTION HLINPROF_SOURCE()
C
C  Which profile the last call of HLINPROF returned: 1 = Stehle & 
C  Hutcheon table (convolved as described in HLINPROF), 2 = VCS table
C  (ISTARK = 2), 0 = HLINOP (no table for the line, or temperature or 
C  density outside it), -1 = no call yet.  For codes that report the 
C  lines that fall back to HLINOP.
C
      INTEGER NFSET,ISRC
      COMMON /HLFILC/ NFSET,ISRC
      EXTERNAL HLPRFD
      HLINPROF_SOURCE = ISRC
      RETURN
      END

C***********************************************************************
      SUBROUTINE HLINPROF_SET_MESSAGES(MODE)
C
C  MODE = 1 (default): HLINPROF prints which tables it uses ("Will use
C  Stehle's tables ...", " Using Hlinop code!") at the first use; 
C  MODE = 0: no such messages (error messages are always printed).
C
      INTEGER MODE
      INTEGER IMSG
      COMMON /HLMSG/ IMSG
      EXTERNAL HLPRFD
      IMSG = MODE
      RETURN
      END

C***********************************************************************
      SUBROUTINE HLINOP_SET_TB09(TBFILE,MODE,DLO,DHI)
C
C  Balmer lines H10-H22: blend the Stehle & Hutcheon profile with the
C  TB09 profile of Tremblay & Bergeron (2009, ApJ 696, 1755; table of
C  November 2015 as distributed with SYNSPEC, data/tremblay.dat, read
C  unchanged by HTB09P).  TB09 includes the non-ideal (Hummer &
C  Mihalas critical field) cut of the proton microfield and is
C  tabulated ~0.5-1.5 dex beyond Stehle's last density (Inglis-Teller
C  limit); it is meant to be used with the occupation probability of
C  the charged perturbers.
C  MODE = 0 (default): Stehle only.  MODE = 1: weight of TB09
C     w = 3t**2 - 2t**3,  t = (log n_e - (L - DLO))/(DLO - DHI) in [0,1]
C  with L = log10 of the line's last tabulated Stehle density, i.e.
C  Stehle below L - DLO, TB09 above L - DHI (DLO > DHI >= 0, in dex).
C  The profile is (1-w) Stehle + w TB09 (both of unit area).  Where
C  TB09 has no value (outside its T, n_e grid) Stehle is used.
C  TBFILE: SYNSPEC's tremblay.dat.  The table is read once, at the
C  first blended profile (a later call with another file name does not
C  reload it); a file that cannot be read stops the run.  The call
C  empties the profile caches.
C  Fit tests (Nov 2026; UVES metal-poor dwarfs and DESI turnoff stars,
C  H10-H22): DLO = 2.0, DHI = 1.5 fitted best.
C
      CHARACTER*(*) TBFILE
      INTEGER MODE
      REAL DLO,DHI
      INTEGER ITBMOD
      COMMON /HLTB9N/ ITBMOD,TBLO,TBHI
      CHARACTER*592 FTB09
      COMMON /HLTB9F/ FTB09
      EXTERNAL HLPRFD
      IF (LEN(TBFILE).GT.592)
     *  STOP 'HLINOP_SET_TB09: file name longer than 592 characters'
      IF (MODE.NE.0 .AND. MODE.NE.1)
     *  STOP 'HLINOP_SET_TB09: MODE must be 0 or 1'
      IF (MODE.EQ.1 .AND. .NOT.(DLO.GT.DHI .AND. DHI.GE.0.))
     *  STOP 'HLINOP_SET_TB09: need DLO > DHI >= 0'
      ITBMOD = MODE
      TBLO = DLO
      TBHI = DHI
      FTB09 = TBFILE
C  the profile caches (HLPCONV) do not know the setting: empty them
      CALL HLINOP_CLEAR_CACHE
      RETURN
      END

C***********************************************************************
      BLOCK DATA HLPRFD
C  No table files set at run time; no HLINPROF call yet; messages on;
C  no TB09 blending.
      CHARACTER*592 FSTE,FSELF,FVCS
      COMMON /HLFILN/ FSTE,FSELF,FVCS
      INTEGER NFSET,ISRC
      COMMON /HLFILC/ NFSET,ISRC
      INTEGER IMSG
      COMMON /HLMSG/ IMSG
      INTEGER ITBMOD
      COMMON /HLTB9N/ ITBMOD,TBLO,TBHI
      CHARACTER*592 FTB09
      COMMON /HLTB9F/ FTB09
      DATA FSTE,FSELF,FVCS /3*' '/
      DATA NFSET,ISRC /0,-1/
      DATA IMSG /1/
      DATA ITBMOD,TBLO,TBHI /0,2.0,1.5/
      DATA FTB09 /' '/
      END

      REAL*8 FUNCTION HLYAPRF(HFILE,WAVE,WAVE0,TEMP,XXNE,H1FRC,HE1FRC,
     *                        AEXT,GNAT)
C
C  Lyman alpha profile (per A):
C  - core: the Stehle & Hutcheon table convolved with the extra Doppler
C    width AEXT and the Lorentz profile of the natural (GNAT, A), 
C    resonance (H-H, Ali & Griem) and van der Waals (He) widths, as in 
C    HPROFL (HLPCONV).  Redward of 1240 A the natural width is left out
C    (RAYLCUT in HPROFL: the far red natural wing is Rayleigh 
C    scattering, treated elsewhere), as a second cached state;
C  - red side: times the Boltzmann factor exp(-hc dsigma/kT) for the
C    absorption wing (Stehle 1994, A&AS 104, 509, eqn 7), as STARK1;
C  - red wing beyond -4000 cm^-1 (1278 A): the Allard et al. H-H and 
C    H-H+ quasi-molecular profiles (HLYAFW) plus half the tabulated 
C    quasistatic wing (times the Boltzmann factor) for the electrons:
C    the H-H+ profiles are for protons only, the tabulated far wing 
C    counts ions and electrons; 0 redward of 5000 A, as in HPROFL;
C  - in between, the red core profile times a factor going linearly 
C    from 1 at the line centre to the value that joins the far wing at
C    -4000 cm^-1 (as HPROFL joins its impact and Holtsmark parts there).
C
      IMPLICIT NONE
      CHARACTER HFILE*(*)
      REAL*8 WAVE,WAVE0,AEXT,GNAT
      REAL TEMP,XXNE,H1FRC,HE1FRC
      REAL*8 CA,RAYLCUT,HCK
      PARAMETER (CA=2.997925D18, RAYLCUT=1240.D0, HCK=1.438777D0)
      REAL*8 GNM,RESONT,VDW,T3NHE,GS,D,SIG0,W4,D4,C4,F4,CORE,BZ
      REAL*8 BZ4,HLPCONV,HLYAFW
      REAL HFNM,PV
C  resonance (Ali & Griem) and van der Waals widths of Lyman alpha as 
C  in HPROFL (dnu/nu), converted to A
      GNM = 0.75D0
      RESONT = DBLE(HFNM(1,2))/(1.D0-0.25D0)*2.07D-24/GNM
      VDW = 4.45D-26/GNM*(4.D0*(7.D0*4.D0+5.D0))**0.4D0
      T3NHE = (DBLE(TEMP)/1.D4)**0.3D0*DBLE(HE1FRC)
      GS = (RESONT*DBLE(H1FRC) + VDW*T3NHE)*WAVE0
      D = WAVE - WAVE0
      IF (D.LE.0.D0) THEN
        HLYAPRF = HLPCONV(HFILE,WAVE,WAVE0,TEMP,XXNE,1,2,AEXT,GS+GNAT)
        RETURN
      ENDIF
      IF (WAVE.GE.5000.D0) THEN
        HLYAPRF = 0.D0
        RETURN
      ENDIF
      SIG0 = 1.D8/WAVE0
      BZ = EXP(-HCK*(SIG0-1.D8/WAVE)/DBLE(TEMP))
C  -4000 cm^-1 from 82259.105 cm^-1, the anchor of the Allard tables
      W4 = 1.D8/(82259.105D0-4000.D0)
      D4 = W4 - WAVE0
      IF (WAVE.LT.W4) THEN
        IF (WAVE.LT.RAYLCUT) THEN
          CORE = HLPCONV(HFILE,WAVE,WAVE0,TEMP,XXNE,1,2,AEXT,GS+GNAT)
        ELSE
          CORE = HLPCONV(HFILE,WAVE,WAVE0,TEMP,XXNE,1,2,AEXT,GS)
        ENDIF
        BZ4 = EXP(-HCK*4000.D0/DBLE(TEMP))
        C4 = HLPCONV(HFILE,W4,WAVE0,TEMP,XXNE,1,2,AEXT,GS)*BZ4
        CALL STEHLE(HFILE,W4,WAVE0,TEMP,XXNE,1,2,PV)
        F4 = HLYAFW(W4,XXNE,H1FRC)*CA/W4**2 + 0.5D0*DBLE(PV)*BZ4
        HLYAPRF = CORE*BZ*(1.D0 + (F4/C4-1.D0)*D/D4)
      ELSE
        CALL STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,1,2,PV)
        HLYAPRF = HLYAFW(WAVE,XXNE,H1FRC)*CA/WAVE**2 + 0.5D0*DBLE(PV)*BZ
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION HLPCONV(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,
     *                        AEXT,GAML)
C
C  Tabulated profile (STEHLE, per A) convolved with a Gaussian of 1/e 
C  half width AEXT (A) and a Lorentzian of half width GAML (A).  The 
C  Gaussian and Lorentzian enter through their analytic Fourier 
C  transforms; the tabulated profile is sampled on a uniform grid and
C  scaled to its own area within the grid (log-grid quadrature; the 
C  tables are not renormalised to 1).  Done once per (NLOW, NUP, WAVE0,
C  TEMP, XXNE, AEXT, GAML) and cached on a logarithmic grid for 
C  |dlambda| <= RC (40 widths, at most 0.15 WAVE0); beyond RC the 
C  tabulated profile and the Lorentzian are added.  The cache holds up
C  to NENT = 3/4 NSLOT states (hash index HLCKEY, stored densely), 
C  enough for every line and depth of a model atmosphere; it is 
C  emptied when full.
C
      IMPLICIT NONE
      CHARACTER HFILE*(*)
      INTEGER NLOW,NUP
      REAL*8 WAVE,WAVE0,AEXT,GAML
      REAL TEMP,XXNE
      INTEGER NSLOT,NOUT,NMAXG,NFFT,NLG,NK
      INTEGER NENT
      PARAMETER (NSLOT=8192, NOUT=1200, NMAXG=8191, NFFT=32768, NLG=600,
     *           NK=7, NENT=3*NSLOT/4)
      REAL*8 PI
      PARAMETER (PI=3.14159265358979D0)
      INTEGER NUSED,IGEN,IS,I,K,NG,NH,NF,IDX(NSLOT)
      LOGICAL HIT
      REAL*8 KEY(NK),KEYS(NK,NENT),LHS(NENT),LRCS(NENT),RCS(NENT)
      REAL*8 Y0(NENT),CV(NMAXG)
      REAL YB(NOUT,NENT),YR(NOUT,NENT)
      REAL*8 RA(NFFT)
      COMPLEX*16 CS(NFFT/2+1),ZW(NFFT/2)
      REAL P0,PV,PV2,PR1,PB1
      REAL*8 WT,W,WMIN,RC,G,H,X,AK,SUMS,AREA,XMIN,X1,X2,R,WW,YL,YH,PHI
      REAL*8 AD,D,XO
      SAVE
      DATA IDX/NSLOT*0/, NUSED/0/, IGEN/0/
C
      KEY(1) = DBLE(NLOW)
      KEY(2) = DBLE(NUP)
      KEY(3) = WAVE0
      KEY(4) = DBLE(TEMP)
      KEY(5) = DBLE(XXNE)
      KEY(6) = AEXT
      KEY(7) = GAML
      CALL HLCKEY(KEY,NK,KEYS,NSLOT,IDX,NENT,NUSED,IGEN,IS,HIT)
      IF (HIT) GOTO 50
C
C  Build.  Width of the tabulated profile from its peak.
C
      CALL STEHLE(HFILE,WAVE0,WAVE0,TEMP,XXNE,NLOW,NUP,P0)
      WT = 0.5D0/DBLE(P0)
      W = MAX(WT,AEXT,GAML)
      WMIN = WT
      IF (AEXT.GT.0.D0) WMIN = MIN(WMIN,AEXT)
      IF (GAML.GT.0.D0) WMIN = MIN(WMIN,GAML)
      RC = MIN(40.D0*W, 0.15D0*WAVE0)
      G = 2.D0*RC
      H = WMIN/8.D0
C  (compare in floating point: G/H can exceed the integer range when a
C  component is extremely narrow)
      IF (G/H.GT.DBLE((NMAXG-1)/2)) THEN
        NH = (NMAXG-1)/2
        H = G/DBLE(NH)
      ELSE
        NH = INT(G/H)
      ENDIF
      NG = 2*NH + 1
      NF = 1
   15 NF = 2*NF
      IF (NF.LT.4*NG) GOTO 15
      IF (NF.GT.NFFT) STOP 'HLPCONV: NFFT too small'
      DO 20 I = 1, NF
        RA(I) = 0.D0
   20 CONTINUE
C  (both sides at once, STEHL2)
      SUMS = 0.D0
      DO 21 I = NH+1, NG
        X = DBLE(I-NH-1)*H
        CALL STEHL2(HFILE,WAVE0+X,WAVE0,TEMP,XXNE,NLOW,NUP,PV,PV2,1)
        RA(I) = DBLE(MAX(PV,0.))
        RA(2*NH+2-I) = DBLE(MAX(PV2,0.))
   21 CONTINUE
      DO 19 I = 1, NG
        SUMS = SUMS + RA(I)
   19 CONTINUE
      SUMS = SUMS*H
C  area of the tabulated profile within +-G (log grid, both sides)
      XMIN = MIN(1.D-4, H/4.D0)
      AREA = 2.D0*DBLE(P0)*XMIN
      X1 = XMIN
      CALL STEHL2(HFILE,WAVE0+X1,WAVE0,TEMP,XXNE,NLOW,NUP,PR1,PB1,1)
      DO 22 K = 1, NLG-1
        X1 = XMIN*(G/XMIN)**(DBLE(K-1)/DBLE(NLG-1))
        X2 = XMIN*(G/XMIN)**(DBLE(K)/DBLE(NLG-1))
        CALL STEHL2(HFILE,WAVE0+X2,WAVE0,TEMP,XXNE,NLOW,NUP,PV,PV2,1)
        AREA = AREA + 0.5D0*(DBLE(PR1)+DBLE(PV))*(X2-X1)
        PR1 = PV
        AREA = AREA + 0.5D0*(DBLE(PB1)+DBLE(PV2))*(X2-X1)
        PB1 = PV2
   22 CONTINUE
      IF (SUMS.GT.0.D0) THEN
        DO 23 I = 1, NG
          RA(I) = RA(I)*AREA/SUMS
   23   CONTINUE
      ENDIF
      CALL HLRFFT(RA,CS,NF,1,ZW)
      DO 25 I = 1, NF/2+1
        K = I - 1
        AK = DBLE(K)/(DBLE(NF)*H)
        CS(I) = CS(I)*EXP(-(PI*AEXT*AK)**2 - 2.D0*PI*GAML*ABS(AK))
   25 CONTINUE
      CALL HLRFFT(RA,CS,NF,-1,ZW)
      DO 30 I = 1, NG
        CV(I) = MAX(RA(I)/DBLE(NF), 0.D0)
   30 CONTINUE
C  store on a log grid of |dlambda|; YR: dlambda > 0, YB: dlambda < 0
      Y0(IS) = CV(NH+1)
      DO 40 K = 1, NOUT
        XO = H*(RC/H)**(DBLE(K-1)/DBLE(NOUT-1))
        R = XO/H
        I = INT(R)
        WW = R - DBLE(I)
        YR(K,IS) = REAL((1.D0-WW)*LOG(MAX(CV(NH+1+I),1.D-300))+WW*
     *            LOG(MAX(CV(NH+2+I),1.D-300)))
        YB(K,IS) = REAL((1.D0-WW)*LOG(MAX(CV(NH+1-I),1.D-300))+WW*
     *            LOG(MAX(CV(NH-I),1.D-300)))
   40 CONTINUE
      LHS(IS) = LOG(H)
      LRCS(IS) = LOG(RC)
      RCS(IS) = RC
C
C  Evaluate
C
   50 D = WAVE - WAVE0
      AD = ABS(D)
      IF (AD.LE.RCS(IS)) THEN
        IF (AD.LT.EXP(LHS(IS))) THEN
          IF (D.GE.0.D0) THEN
            YH = EXP(DBLE(YR(1,IS)))
          ELSE
            YH = EXP(DBLE(YB(1,IS)))
          ENDIF
          PHI = Y0(IS) + (YH-Y0(IS))*AD/EXP(LHS(IS))
        ELSE
          R = 1.D0 + DBLE(NOUT-1)*(LOG(AD)-LHS(IS))/(LRCS(IS)-LHS(IS))
          K = MIN(INT(R),NOUT-1)
          WW = R - DBLE(K)
          IF (D.GE.0.D0) THEN
            YL = DBLE(YR(K,IS))
            YH = DBLE(YR(K+1,IS))
          ELSE
            YL = DBLE(YB(K,IS))
            YH = DBLE(YB(K+1,IS))
          ENDIF
          PHI = EXP((1.D0-WW)*YL + WW*YH)
        ENDIF
        HLPCONV = PHI
      ELSE
        CALL STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,PV)
        HLPCONV = DBLE(PV) + GAML/PI/(GAML*GAML+D*D)
      ENDIF
      RETURN
      END

      SUBROUTINE STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
C
C  Tabulated profile at WAVE (see STEHL2).
C
      CHARACTER HFILE*(*)
      DOUBLE PRECISION WAVE,WAVE0
      REAL TEMP,XXNE,HLIN,HDUM
      INTEGER NLOW,NUP
      CALL STEHL2(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN,HDUM,0)
      RETURN
      END

C***********************************************************************
      SUBROUTINE STEHLS(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN,
     *                  HLINM,IBOTH)
C
C  N Piskunov's code to return Stehle Stark profiles
C
C  Nov 2026: called through STEHL2, which can blend the result with
C  the TB09 profile (HLINOP_SET_TB09).  Sets ANESL (COMMON /HLSLIM/),
C  log10 of the line's last tabulated density at the bracketing
C  temperatures.
C
C  Nov 2026: reads the table built by tools/make_stehle_table.py from 
C  the CDS catalogue J/A+AS/140/93 (format version 2: Lyman, Balmer 
C  and Paschen to n_up = 30, 10 temperatures, log n_e 10 - 19.5).  
C  NALPHA = 0 marks profiles that are not tabulated (above a line's 
C  Inglis-Teller limit, and the lowest temperatures at high density).
C  Electron densities:
C    - inside the tabulated range: interpolation as before;
C    - above the last tabulated density of the line (at both bracketing
C      temperatures): the profile of that density is used in alpha 
C      units, with the normal Holtsmark field F0 scaled to the actual 
C      density as n_e**(2/3), i.e. quasistatic scaling of the width; 
C      the profile stays normalised;
C    - below the first density (1e10): the first profile is used in
C      wavelength units (Doppler-dominated cores).
C  Outside the temperature range HLIN = -1 is returned (caller falls 
C  back to HLINOP).
C  Nov 2026: with IBOTH = 1 also HLINM, the profile at the mirror point
C  2*WAVE0 - WAVE (the tabulated profiles depend on |dlambda| only; the
C  extrapolated far wing has a factor sqrt(lambda/lambda0)), so that a
C  profile on a symmetric grid takes half the work.
C
      DOUBLE PRECISION WAVE,WAVE0
      PARAMETER (NWL=60,NTEMP=10,NNE=20,NLINE=90,SQRTPI=1.77245385)
      REAL T(NTEMP),XNE(NNE)
      REAL WING(NNE,NLINE),F0(NNE,NLINE),WIDTH(NTEMP,NNE,NLINE),
     *     ALPHA(NWL,NTEMP,NNE,NLINE),SPROF(NWL,NTEMP,NNE,NLINE)
      INTEGER NALPHA(NTEMP,NNE,NLINE),NWL1,NTEMP1,NNE1,NLINE1,
     *     N1(NLINE),M1(NLINE)
      INTEGER*2 LINDEX(100,100)
      LOGICAL FIRST
      CHARACTER HFILE*(*)
      INTEGER CONVOLUTION_MODE
      COMMON /HLCONV/ CONVOLUTION_MODE
      REAL*8 STHPCH,STHPD
C  log10 alpha, log10 S and the monotone-cubic slopes, computed once 
C  when the table is read (the interpolation then needs no logarithms
C  of table values and finds the interval by bisection)
      REAL LALPHA(NWL,NTEMP,NNE,NLINE),LSPROF(NWL,NTEMP,NNE,NLINE)
      REAL*8 SLOPE(NWL,NTEMP,NNE,NLINE),LALD(NWL,NTEMP,NNE,NLINE),
     *       LSPD(NWL,NTEMP,NNE,NLINE)
C  the temperature and density part of the last call
      INTEGER LLINE,LITEMP,LINE_,LJMAX
      REAL LTEMP,LXXNE,LDT1,LDT2,LDNE1,LDNE2,LF0A,LF0B
      INTEGER LINE2
      INTEGER ICOR,ITC(4),INC(4),KLO,KHI,KMID,NA
      REAL WCOR(4),F0C(4),SC(4),SCM(4),SM
      INTEGER IBOTH
      REAL HLINM
      SAVE T,XNE,WING,F0,WIDTH,ALPHA,SPROF,NALPHA,
     *     NWL1,NTEMP1,NNE1,NLINE1,FIRST,LINDEX,LALPHA,LSPROF,SLOPE,
     *     LALD,LSPD,
     *     LLINE,LITEMP,LINE_,LJMAX,LTEMP,LXXNE,LDT1,LDT2,LDNE1,LDNE2,
     *     LF0A,LF0B
      INTEGER IMSG
      COMMON /HLMSG/ IMSG
      COMMON /HLSLIM/ ANESL
      EXTERNAL HLPRFD
      DATA FIRST/.TRUE./,LINDEX/10000*0/,LLINE/0/
      EXP10(X)=EXP(2.30258509299405E0*X)
C
C  Read in the table (this is done only once)
C
      IF(FIRST) THEN
        FIRST=.FALSE.
C
C  Open file and process errors
C
        OPEN(1,file=HFILE,IOSTAT=IERR,FORM='UNFORMATTED',STATUS='OLD')
        IF(IERR.NE.0) THEN
          WRITE(*,*) 'ERROR: STEHLE did not find Hydrogen line file'
          WRITE(*,*) HFILE
          STOP
        ENDIF
C
C  Read table dimensions and process errors
C
        READ(1,ERR=11,END=11) NWL1,NTEMP1,NNE1,NLINE1,IVERS
        IF(IVERS.NE.2) GO TO 11
        IF(NWL1.GT.NWL) THEN
          WRITE(*,*) 'ERROR: NWL must be at least',NWL1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
        IF(NTEMP1.GT.NTEMP) THEN
          WRITE(*,*) 'ERROR: NTEMP must be at least',NTEMP1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
        IF(NNE1.GT.NNE) THEN
          WRITE(*,*) 'ERROR: NNE must be at least',NNE1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
        IF(NLINE1.GT.NLINE) THEN
          WRITE(*,*) 'ERROR: NLINE must be at least',NLINE1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
C
C  Read line ID and construct line index
C
        READ(1,ERR=9,END=10) (N1(LINE),M1(LINE),LINE=1,NLINE1)
        DO 1 LINE=1,NLINE1
   1    LINDEX(N1(LINE),M1(LINE))=LINE
C
C  Read temperatures and electron densities
C
        READ(1,ERR=9,END=10) (T(I),I=1,NTEMP1)
        READ(1,ERR=9,END=10) (XNE(J),J=1,NNE1)
C
C  Read wavelength conversion factor, wing extrapolation parameter
C  and line widths
C
        READ(1,ERR=9,END=10) ((F0(J,L),J=1,NNE1),L=1,NLINE1)
        READ(1,ERR=9,END=10) ((WING(J,L),J=1,NNE1),L=1,NLINE1)
        READ(1,ERR=9,END=10) (((WIDTH(I,J,L),I=1,NTEMP1),J=1,NNE1),
     *                       L=1,NLINE1)
C
C  Read the number of detuning points and the two big arrays
C  for detuning (ALPHA) and line profiles (SPROF)
C
        READ(1,ERR=9,END=10) (((NALPHA(I,J,L),I=1,NTEMP1),J=1,NNE1), 
     *                       L=1,NLINE1)
        DO 2 L=1,NLINE1
        READ(1,ERR=9,END=10) (((ALPHA(K,I,J,L),K=1,NWL1),
     *                       I=1,NTEMP1),J=1,NNE1)
   2    READ(1,ERR=9,END=10) (((SPROF(K,I,J,L),K=1,NWL1),
     *                       I=1,NTEMP1),J=1,NNE1)
        CLOSE(1)
C  Guard against trailing padding (alpha = 0 or zero intensity) in a
C  profile; the table built from CDS has none, older tables did.
        DO 23 L=1,NLINE1
        DO 23 J=1,NNE1
        DO 23 I=1,NTEMP1
        DO 22 K=1,NALPHA(I,J,L)
        IF(SPROF(K,I,J,L).LE.0.) THEN
          NALPHA(I,J,L)=K-1
          GO TO 23
        END IF
C  (separate test: .AND. need not short-circuit, ALPHA(0) out of bounds)
        IF(K.GT.1) THEN
          IF(ALPHA(K,I,J,L).LE.ALPHA(K-1,I,J,L)) THEN
            NALPHA(I,J,L)=K-1
            GO TO 23
          END IF
        END IF
   22   CONTINUE
   23   CONTINUE
        DO 26 L=1,NLINE1
        DO 26 J=1,NNE1
        DO 26 I=1,NTEMP1
        NA=NALPHA(I,J,L)
        DO 24 K=2,NA
        LALPHA(K,I,J,L)=LOG10(ALPHA(K,I,J,L))
        LSPROF(K,I,J,L)=LOG10(SPROF(K,I,J,L))
        LALD(K,I,J,L)=LOG10(DBLE(ALPHA(K,I,J,L)))
        LSPD(K,I,J,L)=LOG10(DBLE(SPROF(K,I,J,L)))
   24   CONTINUE
        DO 25 K=2,NA
        SLOPE(K,I,J,L)=STHPD(K,NA,ALPHA(1,I,J,L),SPROF(1,I,J,L))
   25   CONTINUE
   26   CONTINUE
        IF(LINDEX(NLOW,NUP).GT.0.AND.IMSG.EQ.1) WRITE(*,*) 
     *    'Will use Stehle''s tables for Hydrogen line profiles'
      END IF
C
C  Interpolation section
C
      LINE=LINDEX(NLOW,NUP)
C
C  Check if line is in the table
C
      IF(LINE.EQ.0) THEN
        HLIN=-1.
        HLINM=-1.
        RETURN
      END IF
C
C  Temperature and density part (kept from the last call when the line,
C  temperature and density are the same)
C
      IF((LINE.EQ.LLINE).AND.(TEMP.EQ.LTEMP).AND.(XXNE.EQ.LXXNE)) THEN
        IF(LITEMP.EQ.0) THEN
          HLIN=-1.
          HLINM=-1.
          RETURN
        END IF
        ITEMP=LITEMP
        INE=LINE_
        DT1=LDT1
        DT2=LDT2
        DNE1=LDNE1
        DNE2=LDNE2
        F0A=LF0A
        F0B=LF0B
        GO TO 30
      END IF
      LLINE=LINE
      LTEMP=TEMP
      LXXNE=XXNE
      LITEMP=0
C
C  Temperature interpolation (no extrapolation)
C
      IF(TEMP.LT.T(1).OR.TEMP.GT.T(NTEMP1)) THEN
        HLIN=-1.
        HLINM=-1.
        RETURN
      END IF
      ITEMP=2
   3  IF(T(ITEMP).LT.TEMP.AND.ITEMP.LT.NTEMP1) THEN
        ITEMP=ITEMP+1
        GO TO 3
      END IF
      DT=LOG10(T(ITEMP))-LOG10(T(ITEMP-1))
      DT1=(LOG10(TEMP)-LOG10(T(ITEMP-1)))/DT
      DT2=(LOG10(T(ITEMP))-LOG10(TEMP))/DT
C
C  Electron density.  JMAX = number of leading densities tabulated at
C  both bracketing temperatures.
      JMAX=0
   21 IF(JMAX.LT.NNE1) THEN
        IF(NALPHA(ITEMP-1,JMAX+1,LINE).GT.0.AND.
     *     NALPHA(ITEMP,JMAX+1,LINE).GT.0) THEN
          JMAX=JMAX+1
          GO TO 21
        END IF
      END IF
      IF(JMAX.EQ.0) THEN
        HLIN=-1.
        HLINM=-1.
        RETURN
      END IF
      ANESL=XNE(JMAX)
      ANE=LOG10(XXNE)
      IF(ANE.LE.XNE(1).OR.JMAX.EQ.1) THEN
C  Below the table (or only one density): first profile, in A
        INE=2
        DNE1=0.
        DNE2=1.
        F0A=F0(1,LINE)
        F0B=F0A
        IF(ANE.GT.XNE(1)) F0A=F0(1,LINE)*EXP10((ANE-XNE(1))/1.5)
      ELSE IF(ANE.GE.XNE(JMAX)) THEN
C  Above the line's last tabulated density: quasistatic scaling
        INE=JMAX
        DNE1=1.
        DNE2=0.
        F0B=F0(JMAX,LINE)*EXP10((ANE-XNE(JMAX))/1.5)
        F0A=F0B
      ELSE
        INE=2
    4   IF(XNE(INE).LT.ANE.AND.INE.LT.JMAX) THEN
          INE=INE+1
          GO TO 4
        END IF
        DNE=XNE(INE)-XNE(INE-1)
        DNE1=(ANE-XNE(INE-1))/DNE
        DNE2=(XNE(INE)-ANE)/DNE
        F0A=F0(INE-1,LINE)
        F0B=F0(INE,LINE)
      END IF
      LITEMP=ITEMP
      LINE_=INE
      LDT1=DT1
      LDT2=DT2
      LDNE1=DNE1
      LDNE2=DNE2
      LF0A=F0A
      LF0B=F0B
C  Interpolate over alpha. For each of the four (T, n_e) corners with
C  non-zero density weight, 3 cases:
C  1. alpha is in the first interval - simple linear interpolation
C     (.LE.: at alpha = ALPHA(2) the log branch would use ALPHA(1) = 0)
C  2. alpha is outside the table     - extrapolate using wing parameter
C  3. alpha is in between            - log interpolation: linear in 
C     (log alpha, log S) with HLINOP_SET_CONVOLUTION(0); by default a 
C     monotone cubic (STHPCH), because the linear form lies below the 
C     profile in Doppler-dominated cores (half widths 1-6 per cent too 
C     small, areas 0.96-0.99 instead of 1)
C
   30 ITC(1)=ITEMP-1
      INC(1)=INE-1
      F0C(1)=F0A
      WCOR(1)=DNE2
      ITC(2)=ITEMP
      INC(2)=INE-1
      F0C(2)=F0A
      WCOR(2)=DNE2
      ITC(3)=ITEMP-1
      INC(3)=INE
      F0C(3)=F0B
      WCOR(3)=DNE1
      ITC(4)=ITEMP
      INC(4)=INE
      F0C(4)=F0B
      WCOR(4)=DNE1
      DO 40 ICOR=1,4
      SC(ICOR)=0.
      SCM(ICOR)=0.
      IF(WCOR(ICOR).LE.0.) GO TO 40
      IT=ITC(ICOR)
      IN=INC(ICOR)
      NA=NALPHA(IT,IN,LINE)
      ALP=ABS(WAVE-WAVE0)/F0C(ICOR)
      IF(ALP.LE.ALPHA(2,IT,IN,LINE)) THEN
        DAL=ALPHA(2,IT,IN,LINE)-ALPHA(1,IT,IN,LINE)
        DAL1=(ALP-ALPHA(1,IT,IN,LINE))/DAL
        DAL2=(ALPHA(2,IT,IN,LINE)-ALP)/DAL
        S=DAL2*SPROF(1,IT,IN,LINE)+DAL1*SPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.ALPHA(NA,IT,IN,LINE)) THEN
        S=WING(IN,LINE)/(ALP*ALP*SQRT(ALP))*SQRT(WAVE/WAVE0)
        IF(IBOTH.EQ.1) THEN
          SM=WING(IN,LINE)/(ALP*ALP*SQRT(ALP))*
     *       SQRT((2.D0*WAVE0-WAVE)/WAVE0)
          SCM(ICOR)=LOG10(SM/F0C(ICOR))
        END IF
        SC(ICOR)=LOG10(S/F0C(ICOR))
        GO TO 40
      ELSE
        ALP=LOG10(ALP)
C  first K in 2..NA with log10 alpha(K) >= ALP (NA if none)
        KLO=2
        KHI=NA
   31   IF(KLO.LT.KHI) THEN
          KMID=(KLO+KHI)/2
          IF(LALPHA(KMID,IT,IN,LINE).LT.ALP) THEN
            KLO=KMID+1
          ELSE
            KHI=KMID
          END IF
          GO TO 31
        END IF
        K=MAX(KLO,3)
        IF (CONVOLUTION_MODE.EQ.1) THEN
          S=STHPCH(DBLE(ALP),K,LALD(1,IT,IN,LINE),
     *             LSPD(1,IT,IN,LINE),SLOPE(1,IT,IN,LINE))
        ELSE
          DAL=LALPHA(K,IT,IN,LINE)-LALPHA(K-1,IT,IN,LINE)
          DAL1=(ALP-LALPHA(K-1,IT,IN,LINE))/DAL
          DAL2=(LALPHA(K,IT,IN,LINE)-ALP)/DAL
          S=DAL2*LSPROF(K-1,IT,IN,LINE)+DAL1*LSPROF(K,IT,IN,LINE)
          S=EXP10(S)
        END IF
      END IF
      SC(ICOR)=LOG10(S/F0C(ICOR))
      SCM(ICOR)=SC(ICOR)
   40 CONTINUE
      S11=SC(1)
      S21=SC(2)
      S12=SC(3)
      S22=SC(4)
C  Interpolate over T and Ne
      HLIN=EXP10(DT2*DNE2*S11+DT1*DNE2*S21+DT2*DNE1*S12+DT1*DNE1*S22)
      IF(IBOTH.EQ.1) HLINM=EXP10(DT2*DNE2*SCM(1)+DT1*DNE2*SCM(2)+
     *                     DT2*DNE1*SCM(3)+DT1*DNE1*SCM(4))
      RETURN
C
C  I/O error processing
C
   9  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'ERROR reading binary file '//HFILE(1:I)//' in HTABLE'
      STOP
  10  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'EOF found while reading binary file '//HFILE(1:I)//
     *           ' in HTABLE'
      STOP
  11  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'ERROR: '//HFILE(1:I)//' is not a version-2 Stehle '//
     *           'table; build it with tools/make_stehle_table.py'
      STOP
      END

      SUBROUTINE STEHL2(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN,
     *                  HLINM,IBOTH)
C
C  Tabulated Stark profile (with thermal Doppler) of hydrogen lines:
C  Stehle & Hutcheon (STEHLS), for the Balmer lines of the TB09 table
C  (H10-H22) optionally blended with TB09 (HTB09P) from DLO dex below
C  the line's last tabulated Stehle density upwards (HLINOP_SET_TB09;
C  off by default).  Arguments as STEHLS.
C
      DOUBLE PRECISION WAVE,WAVE0
      CHARACTER HFILE*(*)
      INTEGER NLOW,NUP,IBOTH,IV
      REAL TEMP,XXNE,HLIN,HLINM,HT,W,T
      INTEGER ITBMOD
      COMMON /HLTB9N/ ITBMOD,TBLO,TBHI
      CHARACTER*592 FTB09
      COMMON /HLTB9F/ FTB09
      COMMON /HLSLIM/ ANESL
      LOGICAL EX
      INTEGER ITBCHK
      SAVE ITBCHK
      DATA ITBCHK /0/
      EXTERNAL HLPRFD
      CALL STEHLS(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN,HLINM,IBOTH)
      IF (ITBMOD.EQ.0) RETURN
      IF (ITBCHK.EQ.0) THEN
        ITBCHK = 1
        INQUIRE(FILE=FTB09,EXIST=EX)
        IF (.NOT.EX) THEN
          WRITE(*,*) 'ERROR: HLINOP_SET_TB09: TB09 table not found: ',
     *               FTB09(1:MAX(1,LEN_TRIM(FTB09)))
          STOP
        END IF
      END IF
      IF (HLIN.LT.0.) RETURN
      IF (NLOW.NE.2 .OR. NUP.LT.10 .OR. NUP.GT.22) RETURN
      T = (LOG10(XXNE) - (ANESL - TBLO))/(TBLO - TBHI)
      T = MAX(0., MIN(1., T))
      W = T*T*(3. - 2.*T)
      IF (W.LE.0.) RETURN
      CALL HTB09P(FTB09,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HT,IV)
      IF (IV.NE.1) RETURN
C  TB09 depends on |WAVE - WAVE0| only: the same value at the mirror
      HLIN = (1. - W)*HLIN + W*HT
      IF (IBOTH.EQ.1) HLINM = (1. - W)*HLINM + W*HT
      RETURN
      END

      SUBROUTINE HTB09P(LFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,PROF,VALID)
C
C  TB09 Stark profile (Tremblay & Bergeron 2009) of the line NLOW-NUP,
C  per Angstrom, unit area, at |WAVE - WAVE0|; Doppler-convolved at the
C  table temperature like the Stehle table.  VALID = 0 if the line is
C  not in the table or (T, n_e) is outside its grid.  LFILE: the table
C  as distributed with SYNSPEC (data/tremblay.dat; November 2015
C  version, Lemke 1997 format; Lyman 2-21, Balmer 3-22, Paschen 4-22,
C  Brackett 5-14), read unchanged: number of blocks; per block the
C  number of lines, a header per line (NL NU log alpha0 log ne0 log T0
C  dlog alpha dlog ne dlog T MP MNE MT), then per line a label line
C  'nl= NL nu= NU' and MNE*MT records (n_e outer, T inner), each the
C  TB09 flag (Tremblay's README: 0 computed, non-zero not converged or
C  copied from a neighbouring T or n_e) and MP values of log10 S(alpha)
C  with alpha = dlambda/F0, F0 = 1.25e-9 n_e**(2/3) A (as SYNSPEC).
C  The table is read once; a file that cannot be read stops the run.
C  Interpolation: bilinear in log T, log n_e, linear in log alpha,
C  alpha**-2.5 beyond the last point; the area is computed exactly per
C  (T, n_e) and cached.  As the analysis version of HLINPROF
C  (TB09_BALMER), which read a compact copy of the H10-H22 part.
      IMPLICIT NONE
      INTEGER MAXP,MAXNE,MAXT,MAXL,MAXC
      PARAMETER (MAXP=60,MAXNE=17,MAXT=7,MAXL=80,MAXC=100)
      CHARACTER*(*) LFILE
      REAL*8 WAVE,WAVE0
      REAL*4 TEMP,XXNE,PROF
      INTEGER NLOW,NUP,VALID
      INTEGER NLINE,NL(MAXL),NU(MAXL),MP(MAXL),MNE(MAXL),MT(MAXL)
      INTEGER QFLAG(MAXNE,MAXT,MAXL)
      REAL*8 LA0(MAXL),LNE0(MAXL),LT0(MAXL)
      REAL*8 DLA(MAXL),DNE(MAXL),DT(MAXL)
      REAL*8 LOGS(MAXP,MAXNE,MAXT,MAXL)
      REAL*8 CTEMP(MAXC,MAXL),CNE(MAXC,MAXL),CNORM(MAXC,MAXL)
      INTEGER NCACHE(MAXL)
      LOGICAL FIRST,LOADED
      INTEGER I,J,K,L,Q,IERR,IC,NBLK,NB,IB,LL,IN1,IN2,NL1,NU1
      REAL*8 LOGT,LOGNE,WT,WN,X,XLAST,WA,Y1,Y2,Y
      REAL*8 F0,ALPHA,NORM,SUM,A1,S1,SLOPE,RATIO,TERM,FLAG
      CHARACTER*80 LABEL
      SAVE NLINE,NL,NU,MP,MNE,MT,QFLAG,LA0,LNE0,LT0
      SAVE DLA,DNE,DT,LOGS,FIRST,LOADED
      SAVE CTEMP,CNE,CNORM,NCACHE
      DATA FIRST/.TRUE./,LOADED/.FALSE./,NCACHE/MAXL*0/

      VALID = 0
      PROF = 0.
      IF (FIRST) THEN
        FIRST = .FALSE.
        OPEN(96,FILE=LFILE,STATUS='OLD',FORM='FORMATTED',
     *       ACTION='READ',IOSTAT=IERR)
        IF (IERR.NE.0) GOTO 90
        READ(96,*,IOSTAT=IERR) NBLK
        IF ((IERR.NE.0).OR.(NBLK.LT.1)) GOTO 90
        NLINE = 0
        DO IB=1,NBLK
          READ(96,*,IOSTAT=IERR) NB
          IF ((IERR.NE.0).OR.(NB.LT.1).OR.(NLINE+NB.GT.MAXL)) GOTO 90
          DO L=NLINE+1,NLINE+NB
            READ(96,*,IOSTAT=IERR) NL(L),NU(L),LA0(L),LNE0(L),LT0(L),
     *        DLA(L),DNE(L),DT(L),MP(L),MNE(L),MT(L)
            IF (IERR.NE.0) GOTO 90
            IF ((MP(L).LT.2).OR.(MP(L).GT.MAXP).OR.(MNE(L).LT.2).OR.
     *          (MNE(L).GT.MAXNE).OR.(MT(L).LT.2).OR.
     *          (MT(L).GT.MAXT)) GOTO 90
          ENDDO
          DO L=NLINE+1,NLINE+NB
C  label line 'nl= NL nu= NU' (blank lines skipped): must match header
            LABEL = ' '
            DO LL=1,1000
              READ(96,'(A)',IOSTAT=IERR) LABEL
              IF (IERR.NE.0) GOTO 90
              IF (LEN_TRIM(LABEL).GT.0) GOTO 20
            ENDDO
            GOTO 90
   20       IN1 = INDEX(LABEL,'nl=')
            IN2 = INDEX(LABEL,'nu=')
            IF ((IN1.EQ.0).OR.(IN2.LE.IN1)) GOTO 90
            READ(LABEL(IN1+3:IN2-1),*,IOSTAT=IERR) NL1
            IF (IERR.NE.0) GOTO 90
            READ(LABEL(IN2+3:),*,IOSTAT=IERR) NU1
            IF ((IERR.NE.0).OR.(NL1.NE.NL(L)).OR.(NU1.NE.NU(L)))
     *        GOTO 90
            DO J=1,MNE(L)
              DO I=1,MT(L)
                READ(96,*,IOSTAT=IERR) FLAG,(LOGS(K,J,I,L),K=1,MP(L))
                IF (IERR.NE.0) GOTO 90
                QFLAG(J,I,L) = NINT(FLAG)
              ENDDO
            ENDDO
          ENDDO
          NLINE = NLINE + NB
        ENDDO
        CLOSE(96)
        LOADED = .TRUE.
      ENDIF
      IF (.NOT.LOADED) RETURN

      L = 0
      DO Q=1,NLINE
        IF ((NL(Q).EQ.NLOW).AND.(NU(Q).EQ.NUP)) L=Q
      ENDDO
      IF (L.EQ.0) RETURN
      IF ((TEMP.LE.0.).OR.(XXNE.LE.0.)) RETURN
      LOGT = DLOG10(DBLE(TEMP))
      LOGNE = DLOG10(DBLE(XXNE))
      IF ((LOGT.LT.LT0(L)).OR.
     *    (LOGT.GT.LT0(L)+DT(L)*DBLE(MT(L)-1))) RETURN
      IF ((LOGNE.LT.LNE0(L)).OR.
     *    (LOGNE.GT.LNE0(L)+DNE(L)*DBLE(MNE(L)-1))) RETURN

      I = INT((LOGT-LT0(L))/DT(L))+1
      J = INT((LOGNE-LNE0(L))/DNE(L))+1
      I = MAX(1,MIN(I,MT(L)-1))
      J = MAX(1,MIN(J,MNE(L)-1))
      WT = (LOGT-(LT0(L)+DT(L)*DBLE(I-1)))/DT(L)
      WN = (LOGNE-(LNE0(L)+DNE(L)*DBLE(J-1)))/DNE(L)

C  Cache exact area of the interpolated reduced profile per layer and line.
      IC = 0
      DO Q=1,NCACHE(L)
        IF ((CTEMP(Q,L).EQ.DBLE(TEMP)).AND.
     *      (CNE(Q,L).EQ.DBLE(XXNE))) IC=Q
      ENDDO
      IF (IC.GT.0) THEN
        NORM = CNORM(IC,L)
      ELSE
        Y1 = (1.D0-WN)*((1.D0-WT)*LOGS(1,J,I,L)+
     *       WT*LOGS(1,J,I+1,L))+WN*((1.D0-WT)*LOGS(1,J+1,I,L)+
     *       WT*LOGS(1,J+1,I+1,L))
        A1 = 10.D0**LA0(L)
        S1 = 10.D0**Y1
        SUM = S1*A1
        RATIO = 10.D0**DLA(L)
        DO K=1,MP(L)-1
          Y1 = (1.D0-WN)*((1.D0-WT)*LOGS(K,J,I,L)+
     *         WT*LOGS(K,J,I+1,L))+WN*
     *         ((1.D0-WT)*LOGS(K,J+1,I,L)+
     *         WT*LOGS(K,J+1,I+1,L))
          Y2 = (1.D0-WN)*((1.D0-WT)*LOGS(K+1,J,I,L)+
     *         WT*LOGS(K+1,J,I+1,L))+WN*
     *         ((1.D0-WT)*LOGS(K+1,J+1,I,L)+
     *         WT*LOGS(K+1,J+1,I+1,L))
          A1 = 10.D0**(LA0(L)+DLA(L)*DBLE(K-1))
          S1 = 10.D0**Y1
          SLOPE = (Y2-Y1)/DLA(L)
          IF (DABS(SLOPE+1.D0).LT.1.D-12) THEN
            TERM = S1*A1*DLOG(RATIO)
          ELSE
            TERM = S1*A1*(RATIO**(SLOPE+1.D0)-1.D0)/
     *             (SLOPE+1.D0)
          ENDIF
          SUM = SUM+TERM
        ENDDO
        Y2 = (1.D0-WN)*((1.D0-WT)*LOGS(MP(L),J,I,L)+
     *       WT*LOGS(MP(L),J,I+1,L))+WN*
     *       ((1.D0-WT)*LOGS(MP(L),J+1,I,L)+
     *       WT*LOGS(MP(L),J+1,I+1,L))
        A1 = 10.D0**(LA0(L)+DLA(L)*DBLE(MP(L)-1))
        SUM = SUM+(10.D0**Y2)*A1/1.5D0
        NORM = 2.D0*SUM
        IF ((NORM.LE.0.D0).OR.(NORM.NE.NORM)) RETURN
        IF (NCACHE(L).LT.MAXC) THEN
          NCACHE(L)=NCACHE(L)+1
          CTEMP(NCACHE(L),L)=DBLE(TEMP)
          CNE(NCACHE(L),L)=DBLE(XXNE)
          CNORM(NCACHE(L),L)=NORM
        ENDIF
      ENDIF

      F0 = 1.25D-9*DBLE(XXNE)**(2.D0/3.D0)
      ALPHA = DABS(WAVE-WAVE0)/F0
      IF (ALPHA.LE.0.D0) THEN
        X = LA0(L)
      ELSE
        X = DLOG10(ALPHA)
      ENDIF
      XLAST = LA0(L)+DLA(L)*DBLE(MP(L)-1)
      IF (X.LE.LA0(L)) THEN
        K=1
        WA=0.D0
      ELSE IF (X.GE.XLAST) THEN
        K=MP(L)
        WA=0.D0
      ELSE
        K=INT((X-LA0(L))/DLA(L))+1
        K=MAX(1,MIN(K,MP(L)-1))
        WA=(X-(LA0(L)+DLA(L)*DBLE(K-1)))/DLA(L)
      ENDIF
      Y1 = (1.D0-WN)*((1.D0-WT)*LOGS(K,J,I,L)+
     *     WT*LOGS(K,J,I+1,L))+WN*((1.D0-WT)*LOGS(K,J+1,I,L)+
     *     WT*LOGS(K,J+1,I+1,L))
      IF ((K.LT.MP(L)).AND.(X.GT.LA0(L)).AND.(X.LT.XLAST)) THEN
        Y2 = (1.D0-WN)*((1.D0-WT)*LOGS(K+1,J,I,L)+
     *       WT*LOGS(K+1,J,I+1,L))+WN*
     *       ((1.D0-WT)*LOGS(K+1,J+1,I,L)+
     *       WT*LOGS(K+1,J+1,I+1,L))
        Y = (1.D0-WA)*Y1+WA*Y2
      ELSE IF (X.GE.XLAST) THEN
        Y = Y1-2.5D0*(X-XLAST)
      ELSE
        Y = Y1
      ENDIF
      PROF = REAL(10.D0**Y/F0/NORM)
      VALID = 1
      RETURN
   90 CLOSE(96)
      WRITE(*,*) 'ERROR: HLINOP_SET_TB09: cannot read the TB09 table ',
     *  '(SYNSPEC tremblay.dat format): ',
     *  LFILE(1:MAX(1,LEN_TRIM(LFILE)))
      STOP
      END

      REAL*8 FUNCTION STHPCH(XL,K,XA,YA,DD)
C
C  Monotone cubic Hermite interpolation (Fritsch & Carlson 1980; the 
C  PCHIP form) of log10 S against log10 alpha at XL = log10 alpha, in 
C  the interval of nodes K-1 and K (K-1 >= 2: node 1 is alpha = 0).  
C  XA, YA: log10 alpha and log10 S at the nodes of one profile; DD: the
C  slopes at the nodes (STHPD).  Passes through the nodes and keeps the
C  profile monotonic between them.  Returns S.
C
      IMPLICIT NONE
      INTEGER K
      REAL*8 XL,DD(*),XA(*),YA(*)
      REAL*8 X0,X1,Y0,Y1,H,T,D0,D1
      X0 = XA(K-1)
      X1 = XA(K)
      Y0 = YA(K-1)
      Y1 = YA(K)
      H = X1 - X0
      T = (XL - X0)/H
      D0 = DD(K-1)
      D1 = DD(K)
      STHPCH = 10.D0**((2.D0*T**3-3.D0*T**2+1.D0)*Y0
     *       + (T**3-2.D0*T**2+T)*H*D0 + (-2.D0*T**3+3.D0*T**2)*Y1
     *       + (T**3-T**2)*H*D1)
      RETURN
      END

C***********************************************************************
      REAL*8 FUNCTION STHPD(J,NA,AA,SS)
C
C  Slope d(log10 S)/d(log10 alpha) at node J for STHPCH, from the nodes
C  2..NA: weighted harmonic mean of the neighbouring secants inside, 
C  zero at a local extremum, the three-point form at the ends (as in
C  PCHIP).
C
      IMPLICIT NONE
      INTEGER J,NA
      REAL AA(*),SS(*)
      REAL*8 HM,HP,DM,DP,W1,W2,D
      IF (NA-1.LT.2) THEN
        STHPD = 0.D0
        RETURN
      ENDIF
      IF (NA-1.EQ.2) THEN
        STHPD = (LOG10(DBLE(SS(3)))-LOG10(DBLE(SS(2))))/
     *          (LOG10(DBLE(AA(3)))-LOG10(DBLE(AA(2))))
        RETURN
      ENDIF
      IF ((J.GT.2).AND.(J.LT.NA)) THEN
        HM = LOG10(DBLE(AA(J)))-LOG10(DBLE(AA(J-1)))
        HP = LOG10(DBLE(AA(J+1)))-LOG10(DBLE(AA(J)))
        DM = (LOG10(DBLE(SS(J)))-LOG10(DBLE(SS(J-1))))/HM
        DP = (LOG10(DBLE(SS(J+1)))-LOG10(DBLE(SS(J))))/HP
        IF (DM*DP.LE.0.D0) THEN
          STHPD = 0.D0
        ELSE
          W1 = 2.D0*HP + HM
          W2 = HP + 2.D0*HM
          STHPD = (W1+W2)/(W1/DM + W2/DP)
        ENDIF
        RETURN
      ENDIF
C  end points: J = 2 (first node with alpha > 0) or J = NA
      IF (J.LE.2) THEN
        HM = LOG10(DBLE(AA(3)))-LOG10(DBLE(AA(2)))
        HP = LOG10(DBLE(AA(4)))-LOG10(DBLE(AA(3)))
        DM = (LOG10(DBLE(SS(3)))-LOG10(DBLE(SS(2))))/HM
        DP = (LOG10(DBLE(SS(4)))-LOG10(DBLE(SS(3))))/HP
      ELSE
        HM = LOG10(DBLE(AA(NA)))-LOG10(DBLE(AA(NA-1)))
        HP = LOG10(DBLE(AA(NA-1)))-LOG10(DBLE(AA(NA-2)))
        DM = (LOG10(DBLE(SS(NA)))-LOG10(DBLE(SS(NA-1))))/HM
        DP = (LOG10(DBLE(SS(NA-1)))-LOG10(DBLE(SS(NA-2))))/HP
      ENDIF
      D = ((2.D0*HM+HP)*DM - HM*DP)/(HM+HP)
      IF (D*DM.LE.0.D0) THEN
        D = 0.D0
      ELSE IF ((DM*DP.LE.0.D0).AND.(ABS(D).GT.3.D0*ABS(DM))) THEN
        D = 3.D0*DM
      ENDIF
      STHPD = D
      RETURN
      END

C***********************************************************************
      SUBROUTINE HGRID(HFILE,WAVE,WAVE0,TEMP,NHYD,NLOW,NUP,HLIN)
C
C  N Piskunov's code to return Stehle Stark profiles
C  ADAPTED BY PAUL BARKLEM FOR MY GRIDS OF HYDROGEN PROFILES
C  ALL DOUBLE-PRECISION NOW ALSO
C  MAIN THING IS PROFILES ARE NOT ASSUMED SYMMETRIC
C
C  THERE IS NO WING APPROX YET, JUST SET TO A VERY SMALL NUMBER
C
      IMPLICIT NONE
      REAL*8 WAVE,WAVE0
      INTEGER MNPROFS,MNTEMP,MNNH,MNLINE
      PARAMETER (MNTEMP=10,MNNH=30,MNLINE=10,MNPROFS=200)
      REAL*8 T(MNTEMP),NH(MNNH),NHYD,TEMP,HLIN
      REAL*8 F0(MNNH,MNLINE),
     *     MALPHA(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     PALPHA(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     FALPHA(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     MSPROF(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     PSPROF(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     FSPROF(MNPROFS,MNTEMP,MNNH,MNLINE)
      INTEGER MNALPHA(MNTEMP,MNNH,MNLINE),
     *     PNALPHA(MNTEMP,MNNH,MNLINE),
     *     FNALPHA(MNTEMP,MNNH,MNLINE)
      REAL*8 ALP,DAL,DAL1,DAL2,S11,S12,S21,S22
      REAL*8 DT,DT1,DT2,DNH,DNH1,DNH2
      REAL*8 SQRTPI
      INTEGER NPROFS,NTEMP,NNH,NLINE
      INTEGER NU(MNLINE),NL(MNLINE),NLOW,NUP
      INTEGER IERR,I,K,LINE,IN,IT,ITEMP,INH
      INTEGER MIDDLE
      INTEGER*2 LINDEX(100,100)
      LOGICAL FIRST
      CHARACTER HFILE*(*)
      SAVE T,NH,F0,MALPHA,MSPROF,PALPHA,PSPROF,PNALPHA,MNALPHA,
     *     NPROFS,NTEMP,NNH,NLINE,FIRST,LINDEX
      INTEGER IMSG
      COMMON /HLMSG/ IMSG
      EXTERNAL HLPRFD
      DATA FIRST/.TRUE./,LINDEX/10000*0/
      PARAMETER (SQRTPI=1.77245385)
      REAL*8 DEXP10,X
      DEXP10(X)=EXP(2.30258509299405D0*X)
C
C  Read in the table (this is done only once)
C
      IF(FIRST) THEN
        FIRST=.FALSE.
C
C  Open file and process errors
C
        OPEN(1,file=HFILE,IOSTAT=IERR,FORM='UNFORMATTED',STATUS='OLD')
        IF(IERR.NE.0) THEN
          WRITE(*,*) 'ERROR: HGRID did not find Hydrogen line file'
          WRITE(*,*) HFILE
          STOP
        ENDIF
C
C  Read table dimensions and process errors
C
        READ(1,ERR=9,END=10) NLINE
C
C  Read line ID and construct line index, read grid parameters
C
        READ(1,ERR=9,END=10) (NL(LINE),NU(LINE),LINE=1,NLINE)
        DO 401 I=1,NLINE
 401       LINDEX(NL(I),NU(I))=I
        READ(1,ERR=9,END=10) NNH
        READ(1,ERR=9,END=10) (NH(I),I=1,NNH)
        READ(1,ERR=9,END=10) NTEMP
        READ(1,ERR=9,END=10) (T(I),I=1,NTEMP)
C
C  Read the number of detuning points and the two big arrays
C  for detuning (ALPHA) and line profiles (SPROF)
C  Right Now ALPHA=LAMBDA, SPROF=INTENSITY (normalised)
C  hence F0=1 everywhere 
C
        READ(1,ERR=9,END=10) NPROFS
        MIDDLE=(NPROFS-1)/2+1
        DO LINE=1,NLINE
           DO IN=1,NNH
              DO IT=1,NTEMP
                 DO I=1,NPROFS
          READ(1,ERR=9,END=10) FALPHA(I,IT,IN,LINE),FSPROF(I,IT,IN,LINE)
                 END DO
              END DO
           END DO
        END DO
        DO LINE=1,NLINE
           DO IN=1,NNH
              F0(IN,LINE)=1.d0
           END DO
        END DO
        DO LINE=1,NLINE
           DO IN=1,NNH
              DO ITEMP=1,NTEMP
                 FNALPHA(ITEMP,IN,LINE)=NPROFS
                 PNALPHA(ITEMP,IN,LINE)=MIDDLE
                 MNALPHA(ITEMP,IN,LINE)=MIDDLE
              END DO
           END DO
        END DO
C
C  SET UP SEPARATE ARRAYS OF EACH SIDE OF THE PROFILE
C
        DO LINE=1,NLINE
           DO IN=1,NNH
              DO IT=1,NTEMP
                 DO K=1,MIDDLE
                    PALPHA(K,IT,IN,LINE)=FALPHA(MIDDLE+K-1,IT,IN,LINE)
              MALPHA(K,IT,IN,LINE)=DABS(FALPHA(MIDDLE-K+1,IT,IN,LINE))
                    PSPROF(K,IT,IN,LINE)=FSPROF(MIDDLE+K-1,IT,IN,LINE)
                    MSPROF(K,IT,IN,LINE)=FSPROF(MIDDLE-K+1,IT,IN,LINE)
                 END DO
              END DO
           END DO
        END DO                    
C
C  IF OUR LINE EXIST THEN WRITE OUT FIRST TIME
C
        IF(LINDEX(NLOW,NUP).GT.0.AND.IMSG.EQ.1) WRITE(*,*) 
     *    'Will use PB''s tables for Hydrogen line profiles'
C
      END IF  ! ENDS PART ONLY DONE FIRST TIME
C
C  Interpolation section
C
      LINE=LINDEX(NLOW,NUP)
C
C  Check if line is in the table
C
      IF(LINE.EQ.0) THEN
        HLIN=-1.
        RETURN
      END IF
C
C  Temperature and perturber density interpolation (no extrapolation)
C
      IF(TEMP.LT.T(1).OR.TEMP.GT.T(NTEMP)) THEN
        HLIN=-1.
        RETURN
      END IF
      ITEMP=2
   3  IF(T(ITEMP).LT.TEMP.AND.ITEMP.LT.NTEMP) THEN
        ITEMP=ITEMP+1
        GO TO 3
      END IF
      DT=DLOG10(T(ITEMP))-DLOG10(T(ITEMP-1))
      DT1=(DLOG10(TEMP)-DLOG10(T(ITEMP-1)))/DT
      DT2=(DLOG10(T(ITEMP))-DLOG10(TEMP))/DT
C
      IF(NHYD.LT.NH(1).OR.NHYD.GT.NH(NNH)) THEN
        HLIN=-1.
        RETURN
      END IF
      INH=2
   4  IF(NH(INH).LT.NHYD.AND.INH.LT.NNH) THEN
        INH=INH+1
        GO TO 4
      END IF
      DNH=DLOG10(NH(INH))-DLOG10(NH(INH-1))
      DNH1=(DLOG10(NHYD)-DLOG10(NH(INH-1)))/DNH
      DNH2=(DLOG10(NH(INH))-DLOG10(NHYD))/DNH
C
C  WE MUST SET UP TO ALLOW FOR THE FACT THAT THE PROFILE IS NOT ASSUMED 
C  TO BE SYMMETRIC
C
C  Interpolate over alpha. For each IT and IN we consider 3 cases:
C  1. alpha is in the first interval - simple linear interpolation
C  2. alpha is outside the table     - extrapolate using wing parameter
C  3. alpha is in between            - log interpolation
C
c  no wing extrapolation (impact approx invalid), we just take it as zero.
c
      IF ((WAVE-WAVE0).GE.0.d0) THEN
C
C  CASE OF POSITIVE DETUNING
C
      IT=ITEMP-1
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S11=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S11=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   5    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 5
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S11=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S11=DEXP10(S11)
      END IF
      if (S11.ne.0.d0) S11=DLOG10(S11/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S21=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S21=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   6    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 6
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S21=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S21=DEXP10(S21)
      END IF
      if (S21.ne.0.d0) S21=DLOG10(S21/F0(IN,LINE))
C
      IT=ITEMP-1
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S12=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S12=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   7    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 7
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S12=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S12=DEXP10(S12)
      END IF
      if (S12.ne.0.d0) S12=DLOG10(S12/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S22=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S22=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   8    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 8
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S22=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S22=DEXP10(S22)
      END IF
      if (S22.ne.0.d0) S22=DLOG10(S22/F0(IN,LINE))
C
      ELSE
C
C  NOW THE CASE OF NEGATIVE ALPHA
C
      IT=ITEMP-1
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S11=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S11=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  15    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 15
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S11=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S11=DEXP10(S11)
      END IF
      if (S11.ne.0.d0) S11=DLOG10(S11/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S21=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S21=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  16    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 16
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S21=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S21=DEXP10(S21)
      END IF
      if (S21.ne.0.d0) S21=DLOG10(S21/F0(IN,LINE))
C
      IT=ITEMP-1
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S12=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S12=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  17    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 17
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S12=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S12=DEXP10(S12)
      END IF
      if (S12.ne.0.d0) S12=DLOG10(S12/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S22=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S22=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  18    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 18
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S22=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S22=DEXP10(S22)
      END IF
      if (S22.ne.0.d0) S22=DLOG10(S22/F0(IN,LINE))
C
      END IF
C
C  Interpolate over T and NH
C
      if (S11.eq.0.d0.and.S21.eq.0.d0.and.S12.eq.0.d0.and.S22.eq.0.d0)
     * then
         HLIN=0.d0
      else
      HLIN=DEXP10(DT2*DNH2*S11+DT1*DNH2*S21+DT2*DNH1*S12+DT1*DNH1*S22)
      end if
      RETURN
C
C  I/O error processing
C
   9  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'ERROR reading binary file '//HFILE(1:I)//' in HGRID'
      STOP
  10  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'EOF found while reading binary file '//HFILE(1:I)//
     *           ' in HGRID'
      STOP
      END


      SUBROUTINE HSELF_PDWIDTH(NLOW,NUPP,TEMP,H1FRC,WIDTH)
C
C   Computes the Lorentzian WIDTH self broadening 
C   of the Hydrogen Balmer Lines in the p-d approximation
C   Output in usual rad /s cm^3
C   
      IMPLICIT NONE
      REAL*4 H1FRC
      REAL*4 C,PI,KB,M0,A0,SIGMAH(4),ALPHAH(4),
     *       X,GAMMAF,GX,VBAR,SIGMA,TEMP
      REAL*8 GVW, WIDTH
      INTEGER LINE,NLOW,NUPP
C
      PARAMETER (C=2.997925E+18,PI=3.14159265)
      PARAMETER (KB=1.380658E-23,M0=1.660540E-27)
      PARAMETER (A0=5.29177249E-11)  
      DATA SIGMAH/ 1180., 2320., 4208., 0.0/
      DATA ALPHAH/ 0.677, 0.455, 0.380, 0.0/
C
      IF ((NLOW.NE.2) .OR. (NUPP.GT.5)) THEN
         WIDTH = 0.d0
         RETURN
      ENDIF  
C
      LINE = NUPP - NLOW
      X = 2. - ALPHAH(LINE) * .5
      GX = X - 1.0
      GAMMAF = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX
     ;           )*GX)*GX)*GX)*GX
      SIGMA = SIGMAH(LINE) * A0 * A0
      GVW = (4./PI)**(ALPHAH(LINE)*0.5)*GAMMAF*1.E4*SIGMA
      VBAR = SQRT(8.*KB*TEMP/PI/M0 * 2./1.008)
      GVW = GVW * ((VBAR/1.E4)**(1.-ALPHAH(LINE)))
      GVW = GVW * H1FRC * 1.E6
C 
      WIDTH=GVW
      RETURN
C
      END

      SUBROUTINE HE_PDWIDTH(NLOW,NUPP,TEMP,HE1FRC,HEWID)
C 
C  Compute approximate broadening by He collisions (p-d approx)
C
C  Output in usual rad /s cm^3
C   
C  27/4/00 PB
C
      IMPLICIT NONE
      REAL*4 PI,C,HE1FRC,GAMMAF
      REAL*4 SIGMAH(4),ALPHAH(4),X,GX,A0,SIGMA,VBAR,KB,M0,TEMP
      REAL*4 HEWID,GVW
      INTEGER NLOW,NUPP,LINE
      PARAMETER (C=2.997925E+18,PI=3.14159265,A0=5.29177249E-11)    
      PARAMETER (KB=1.380658E-23,M0=1.660540E-27)
      DATA SIGMAH/ 834., 1998., 3837., 5700. /
      DATA ALPHAH/ .280, .330, .260, .250 /
C
      IF ((NLOW.NE.2) .OR. (NUPP.GT.6)) THEN
         HEWID = 0.d0
         RETURN
      ENDIF  
C
      LINE = NUPP - NLOW
      X = 2. - ALPHAH(LINE) * .5
      GX = X - 1.0
      GAMMAF = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX
     ;           )*GX)*GX)*GX)*GX
      SIGMA = SIGMAH(LINE) * A0 * A0
      GVW = (4./PI)**(ALPHAH(LINE)*0.5)*GAMMAF*1.E4*SIGMA
      VBAR = SQRT(8.*KB*TEMP/PI/M0 * (1./1.008 + 1./4.003))
      HEWID = GVW * ((VBAR/1.E4)**(1.-ALPHAH(LINE))) 
      HEWID = HEWID * 0.625            !polarisability difference term
      HEWID = HEWID * HE1FRC * 1.E6    !in cgs
C
      RETURN
      END

      SUBROUTINE RAD_WIDTH(NLOW,NUPP,WIDTH)
C
C  Returns the radiation damping width (hwhm) in rad/s cm^3
C
C  12/6/2013 PB based on Kurucz and Peterson
C
      INTEGER NLOW,NUPP,LINE
      REAL*4 WIDTH,GAMMA,GL,GU,ASUM(100),ASUMLYMAN(100)
C
C  Einstein A-value sums for H lines
C
      DATA ASUM/
     1 0.000E+00, 4.696E+08, 9.980E+07, 3.017E+07, 1.155E+07, 5.189E+06,
     2 2.616E+06, 1.437E+06, 8.444E+05, 5.234E+05, 3.389E+05, 2.275E+05,
     3 1.575E+05, 1.120E+05, 8.142E+04, 6.040E+04, 4.560E+04, 3.496E+04,
     4 2.719E+04, 2.141E+04, 1.711E+04, 1.377E+04, 1.119E+04, 9.166E+03,
     5 7.572E+03, 6.341E+03, 5.338E+03, 4.523E+03, 3.854E+03, 3.302E+03,
     6 2.844E+03, 2.460E+03, 2.138E+03, 1.866E+03, 1.635E+03, 1.438E+03,
     7 1.269E+03, 1.124E+03, 9.983E+02, 8.894E+02, 7.947E+02, 7.120E+02,
     8 6.396E+02, 5.759E+02, 5.198E+02, 4.703E+02, 4.263E+02, 3.873E+02,
     9 3.526E+02, 3.215E+02, 2.938E+02, 2.689E+02, 2.465E+02, 2.264E+02,
     A 2.082E+02, 1.918E+02, 1.769E+02, 1.634E+02, 1.512E+02, 1.400E+02,
     1 1.298E+02, 1.206E+02, 1.121E+02, 1.043E+02, 9.720E+01, 9.066E+01,
     2 8.465E+01, 7.912E+01, 7.403E+01, 6.933E+01, 6.498E+01, 6.097E+01,
     3 5.725E+01, 5.381E+01, 5.061E+01, 4.765E+01, 4.489E+01, 4.232E+01,
     4 3.994E+01, 3.771E+01, 3.563E+01, 3.369E+01, 3.188E+01, 3.019E+01,
     5 2.860E+01, 2.712E+01, 2.572E+01, 2.442E+01, 2.319E+01, 2.204E+01,
     6 2.096E+01, 1.994E+01, 1.898E+01, 1.808E+01, 1.722E+01, 1.642E+01,
     7 1.566E+01, 1.495E+01, 1.427E+01, 1.363E+01/
C
C  For Lyman lines only the s-p transition is allowed.
C
      DATA ASUMLYMAN/
     1 0.000E+00, 6.265E+08, 1.897E+08, 8.126E+07, 4.203E+07, 2.450E+07,
     2 1.236E+07, 8.249E+06, 5.782E+06, 4.208E+06, 3.158E+06, 2.430E+06,
     3 1.910E+06, 1.567E+06, 1.274E+06, 1.050E+06, 8.752E+05, 7.373E+05,
     4 6.269E+05, 5.375E+05, 4.643E+05, 4.038E+05, 3.534E+05, 3.111E+05,
     5 2.752E+05, 2.447E+05, 2.185E+05, 1.959E+05, 1.763E+05, 1.593E+05,
     6 1.443E+05, 1.312E+05, 1.197E+05, 1.094E+05, 1.003E+05, 9.216E+04,
     7 8.489E+04, 7.836E+04, 7.249E+04, 6.719E+04, 6.239E+04, 5.804E+04,
     8 5.408E+04, 5.048E+04, 4.719E+04, 4.418E+04, 4.142E+04, 3.888E+04,
     9 3.655E+04, 3.440E+04, 3.242E+04, 3.058E+04, 2.888E+04, 2.731E+04,
     A 2.585E+04, 2.449E+04, 2.322E+04, 2.204E+04, 2.094E+04, 1.991E+04,
     1 1.894E+04, 1.804E+04, 1.720E+04, 1.640E+04, 1.566E+04, 1.496E+04,
     2 1.430E+04, 1.368E+04, 1.309E+04, 1.254E+04, 1.201E+04, 1.152E+04,
     3 1.105E+04, 1.061E+04, 1.019E+04, 9.796E+03, 9.419E+03, 9.061E+03,
     4 8.721E+03, 8.398E+03, 8.091E+03, 7.799E+03, 7.520E+03, 7.255E+03,
     5 7.002E+03, 6.760E+03, 6.530E+03, 6.310E+03, 6.100E+03, 5.898E+03,
     6 5.706E+03, 5.522E+03, 5.346E+03, 5.177E+03, 5.015E+03, 4.860E+03,
     7 4.711E+03, 4.569E+03, 4.432E+03, 4.300E+03/
C
       GAMMA = 0.
       GL = 0.
       GU = 0.
       IF (NLOW.EQ.1) THEN
         IF (NUPP.LE.100) GAMMA = ASUMLYMAN(NUPP)
       ELSE
         IF (NLOW.LE.100) GL = ASUM(NLOW)
         IF (NUPP.LE.100) GU = ASUM(NUPP)
         GAMMA = GL + GU
       ENDIF
       WIDTH = GAMMA / 2. 
C
      RETURN
      END


      SUBROUTINE ALI_GRIEM(NLOW,NUPP,WAVE0,H1FRC,WIDTH)
C
C   Computes the Resonance Width (AA) of the Hydrogen Balmer Lines
C   Halfwidth given H1FRC perturbers per cm^3 
C
C   NB: assumes 2p level only! more correct than adding the 3p or 4p etc.
C   since 2p-nd dominates
C   so LINE is not used actually
C
C  GAMMA = lambda^2/(2 pi c) * 1.15 * e^2/(mc) f12 lambda12 N(H)   fuhrmann
C  GAMMA = lambda^2/(2 pi c) * 1.108 * e^2/(mc) f12 lambda12 N(H)   ali and griem
C
      IMPLICIT NONE
      REAL*4 H1FRC
      REAL*8 WAVE0,WIDTH,FW
      INTEGER NLOW,NUPP,LINE
C      
      IF (NLOW.NE.2)  THEN
         WIDTH = 0.d0
         RETURN
      ENDIF  
C
C     WAVE0 in AA, WIDTH in AA  
C
C   the two different leading numbers makes only a minute difference
C     FW=2.6037D-27*WAVE0*WAVE0   !balmer lower level (fuhrmann 2 pi)
C
      FW=2.516D-27*WAVE0*WAVE0   !balmer lower level (ali and griem 1.92 pi)
      FW=FW*H1FRC
C      
      WIDTH=FW/2.d0
C
      RETURN
C
      END


C
C  THE FOLLOWING IS KURUCZ'S SUBROUTINE
C
C  SLIGHTLY ALTERED BY PB & REGNER TRAMPEDACH
C  UNIT 25->1 NEEDED TO WORK WITH MULTI
C
      SUBROUTINE VCS(PR,XNE,T,DEL,II,N,M,HVCSFILE)     
C     CALCULATES VIDAL COOPER AND SMITH PROFILES FOR FIRST FOUR BALMER 
C     LINES.  THEY ARE RETURNED IN ANGSTROM UNITS.                     
C     ASSUMES RAW PROFILES ARE IN ALPHA UNITS.                      
      IMPLICIT REAL*8 (A-H,O-Z)
C  (assumed size: II values; was PR(40), DEL(40))
      DIMENSION PR(*),DEL(*)
C     DIMENSION ALPHA(40),PRALPH(40)                                  
      DIMENSION PRALPH(40)                                   
      DIMENSION SVCS(6,17,40,4),ALPHA0(4)  
      CHARACTER*(*) HVCSFILE                                
      DATA SVCS(1,1,1,1)/0./,ALPHA0/-3.,-3.,-3.,-3./   
      SAVE SVCS
C
      EXP10(X)=EXP(2.30258509299405E0*X)                       
      IF(SVCS(1,1,1,1).NE.0.)GO TO 3
C
C     READ IN VCS ARRAYS                                          
      OPEN(1,FILE=HVCSFILE,FORM='FORMATTED',STATUS='OLD')
      READ(1,10)
      DO 20 LINE=1,4
        READ(1,10)
        READ(1,10)(((SVCS(I,J,K,LINE),J=1,17),I=1,6),K=1,40)
   10   FORMAT (10F8.4)
   20 CONTINUE
      CLOSE(1)
C
    3 LINE=M-N                                  
C     ONLY THE FIRST FOUR BALMER LINES ARE TABULATED
      IF ((N.NE.2).OR.(LINE.GT.4)) THEN
        PR(1)=-1.
        RETURN
      END IF
C     TEMPERATURE AND ELECTRON DENSITY INTERPOLATION        
      AT=LOG10(T)                                                              
      BTEMP=(AT-3.3979411)/.3010300+1.                     
      ITEMP=BTEMP                                       
C     ITEMP=(LOG10(T)-LOG10(2500.))/LOG10(2.)+1.                             
      ITEMP=MAX(MIN(ITEMP,5),1)
      WTTEMP=BTEMP-ITEMP                                   
      ANE=LOG10(XNE)                                                           
C     TO AVOID BAD EXTRAPOLATION AT LOW DENSITY
      ANE=MAX(10.D0,ANE)
      BNE=(ANE-10.)/.5+1.                  
      INE=BNE                      
      INE=MAX(MIN(INE,16),1)
      WTXNE=BNE-INE                             
      DO 21 I=1,40                                                
   21 PRALPH(I)=(1.-WTXNE)*(1.-WTTEMP)*SVCS(ITEMP  ,INE  ,I,LINE)+
     1               (1.-WTXNE)*WTTEMP*SVCS(ITEMP+1,INE  ,I,LINE)+           
     2               WTXNE*(1.-WTTEMP)*SVCS(ITEMP  ,INE+1,I,LINE)+
     3                    WTXNE*WTTEMP*SVCS(ITEMP+1,INE+1,I,LINE)   
C     NOW ALPHA INTERPOLATION                                             
      FO=1.25E-9*XNE**.66666667                                    
      DO 50 I=1,II                                         
      IF (DEL(I).NE.0.) GO TO 25                            
      PR(I)=EXP10(PRALPH(1))/FO                               
      GO TO 50                     
   25 AALP=LOG10(ABS(DEL(I))/FO)                                               
      BALP=(AALP-ALPHA0(LINE))/.2+1.          
      IALP=BALP                                       
      IF (IALP.GT.1) GO TO 30                 
      PR(I)=EXP10(PRALPH(1))/FO           
      GO TO 50                                                  
   30 IF (IALP.GT.39) GO TO 40                
      WT=BALP-IALP                        
      PR(I)=EXP10((1.-WT)*PRALPH(IALP)+WT*PRALPH(IALP+1))/FO    
      GO TO 50                                        
   40 PR(I)=EXP10(PRALPH(40)+2.5*(ALPHA0(LINE)+7.8-AALP))/FO      
   50 CONTINUE
      RETURN                       
      END  
