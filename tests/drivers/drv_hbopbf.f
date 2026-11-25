C  Test driver: HBOP_LEVELS and HBOP_BF against HBOP (hbop.f).
C  stdin, line 1: NH NHE NE T
C         then one wavelength (A) per line
C  stdout, line 1: sum of the 100 LTE populations / NH, W(1), W(15)
C          then per wavelength: WAVE  CONTIN (HBOP, LTE, no lines)  HBOP_BF
      PROGRAM DRVHBF
      INTEGER NLO(1), NUP(1)
      REAL NH, NHE, NE, T, DOP, NPOP(1), TOT, CONT, CB, S
      REAL NLTE(100)
      REAL*8 WAVE, WAVEH(1), W(100)
      INTEGER I
      READ(*,*) NH, NHE, NE, T
      CALL HBOP_LEVELS(NH, NHE, NE, T, 100, W, NLTE)
      S = 0.
      DO 5 I = 1, 100
        S = S + NLTE(I)
    5 CONTINUE
      WRITE(*,'(1P,3E25.16)') S/NH, W(1), W(15)
      DOP = 1.E-5
   10 READ(*,*,END=99) WAVE
      CALL HBOP(WAVE, 0, NLO, NUP, WAVEH, NH, NHE, NE, T, DOP,
     *          NPOP, 0, TOT, CONT)
      CALL HBOP_BF(WAVE, NH, NHE, NE, T, W, NLTE, NLTE, CB)
      WRITE(*,'(1P,3E25.16)') WAVE, CONT, CB
      GOTO 10
   99 CONTINUE
      END
