C  Test driver: HBOP_LEVELS, HBOP_BF and HBOP_BF_DISSOLVED against HBOP
C  (hbop.f).
C  stdin, line 1: NH NHE NE T
C         then one wavelength (A) per line
C  stdout, line 1: sum of the 100 LTE populations / NH, W(1), W(15)
C          then per wavelength: WAVE  CONTIN (HBOP, LTE, no lines)  HBOP_BF
C          HBOP_BF_DISSOLVED with NBFMAX = 15, and with NBFMAX = 4
      PROGRAM DRVHBF
      INTEGER NLO(1), NUP(1)
      REAL NH, NHE, NE, T, DOP, NPOP(1), TOT, CONT, CB, S, D15, D4
      REAL NLTE(100)
      REAL*8 WAVE, WAVEH(1), W(100)
      INTEGER I
      READ(*,*) NH, NHE, NE, T
      CALL HBOP_LEVELS(NH, NHE, NE, T, 100, W, NLTE)
      S = 0.
      DO 5 I = 1, 100
        S = S + NLTE(I)
    5 CONTINUE
      WRITE(*,'(1P,5E25.16)') S/NH, W(1), W(15), 0., 0.
      DOP = 1.E-5
   10 READ(*,*,END=99) WAVE
      CALL HBOP(WAVE, 0, NLO, NUP, WAVEH, NH, NHE, NE, T, DOP,
     *          NPOP, 0, TOT, CONT)
      CALL HBOP_BF(WAVE, NH, NHE, NE, T, W, NLTE, NLTE, CB)
      CALL HBOP_BF_DISSOLVED(WAVE, NH, NHE, NE, T, W, NLTE, NLTE, 15,
     *                       D15)
      CALL HBOP_BF_DISSOLVED(WAVE, NH, NHE, NE, T, W, NLTE, NLTE, 4, D4)
      WRITE(*,'(1P,5E25.16)') WAVE, CONT, CB, D15, D4
      GOTO 10
   99 CONTINUE
      END
