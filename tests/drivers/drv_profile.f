C  Test driver: HPROFL, HLINOP, HLCUT and STARK1 (hlinop.f).
C  stdin, line 1:  N M WAVEH T XNE H1FRC HE1FRC DOP MODE
C         MODE >= 0: CALL HLINOP_SET_STARK(MODE) first (negative: default)
C         then one wavelength (A) per line
C  stdout: WAVE HPROFL HLINOP HLCUT STARK1 (profiles per unit frequency)
      PROGRAM DRVPROF
      INTEGER N, M, MODE
      REAL*8 WAVE, WAVEH
      REAL T, XNE, H1, HE1, DOP, HPROFL, HLINOP, HLCUT, STARK1
      REAL P, Q, W, S
      READ(*,*) N, M, WAVEH, T, XNE, H1, HE1, DOP, MODE
      IF (MODE.GE.0) CALL HLINOP_SET_STARK(MODE)
   10 READ(*,*,END=99) WAVE
      P = HPROFL(N, M, WAVE, WAVEH, T, XNE, H1, HE1, DOP)
      Q = HLINOP(WAVE, N, M, WAVEH, T, XNE, H1, HE1, DOP)
      W = HLCUT(WAVE, N, M, WAVEH)
      S = STARK1(N, M, WAVE, WAVEH, T, XNE)
      WRITE(*,'(1P,5E25.16)') WAVE, P, Q, W, S
      GOTO 10
   99 CONTINUE
      END
