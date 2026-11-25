C  Test driver: HLINOP_SET_MEDIUM (hlinop.f, hbop.f).
C  stdin, line 1: IAIR (1 air, 0 vacuum; negative: do not call)
C         then per line: WAVE NLO NUP WAVEH
C  stdout: HLCUT  HLINOP (T 6000 K, n_e 1e14, n_H 1e16, n_He 1e15)
C          HBOP_BF at WAVE (n_H 1e17, n_He 1e16, n_e 1e14, T 8000 K)
      PROGRAM DRVMED
      INTEGER IAIR, NLO, NUP
      REAL*8 WAVE, WAVEH, W(100)
      REAL HLCUT, HLINOP, CUT, HL, CB, NLTE(100)
      REAL NH, NHE, NE, T
      READ(*,*) IAIR
      IF (IAIR.GE.0) CALL HLINOP_SET_MEDIUM(IAIR)
      NH = 1.E17
      NHE = 1.E16
      NE = 1.E14
      T = 8000.
      CALL HBOP_LEVELS(NH, NHE, NE, T, 100, W, NLTE)
   10 READ(*,*,END=99) WAVE, NLO, NUP, WAVEH
      CUT = HLCUT(WAVE, NLO, NUP, WAVEH)
      HL = HLINOP(WAVE, NLO, NUP, WAVEH, 6000., 1.E14, 1.E16, 1.E15,
     *            1.E-5)
      CALL HBOP_BF(WAVE, NH, NHE, NE, T, W, NLTE, NLTE, CB)
      WRITE(*,'(1P,3E25.16)') CUT, HL, CB
      GOTO 10
   99 CONTINUE
      END
