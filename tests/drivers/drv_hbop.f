C  Test driver: HBOP total and continuous opacity (hbop.f), LTE.
C  stdin, line 1:  NH NHE NE T DOP NLINES
C         then NLINES lines:  NLO NUP WAVEH
C         then one wavelength (A) per line
C  stdout: WAVE TOTAL CONTIN   (cm^-1)
      PROGRAM DRVHBOP
      INTEGER NMAX
      PARAMETER (NMAX = 200)
      INTEGER N, I, NLO(NMAX), NUP(NMAX)
      REAL*8 WAVE, WAVEH(NMAX)
      REAL NH, NHE, NE, T, DOP, NPOP(1), TOTAL, CONTIN
      READ(*,*) NH, NHE, NE, T, DOP, N
      IF (N.GT.NMAX) STOP 'DRVHBOP: too many lines'
      DO 5 I = 1, N
    5 READ(*,*) NLO(I), NUP(I), WAVEH(I)
   10 READ(*,*,END=99) WAVE
      CALL HBOP(WAVE, N, NLO, NUP, WAVEH, NH, NHE, NE, T, DOP,
     *          NPOP, 0, TOTAL, CONTIN)
      WRITE(*,'(1P,3E25.16)') WAVE, TOTAL, CONTIN
      GOTO 10
   99 CONTINUE
      END
