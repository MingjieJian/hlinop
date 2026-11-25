C  Test driver: HPROFL, HLINOP and HLCUT (hlinop.f).
C  stdin, line 1:  N M WAVEH T XNE H1FRC HE1FRC DOP
C         then one wavelength (A) per line
C  stdout: WAVE HPROFL HLINOP HLCUT  (profiles per unit frequency, Hz^-1)
      PROGRAM DRVPROF
      INTEGER N, M
      REAL*8 WAVE, WAVEH
      REAL T, XNE, H1, HE1, DOP, HPROFL, HLINOP, HLCUT, P, Q, W
      READ(*,*) N, M, WAVEH, T, XNE, H1, HE1, DOP
   10 READ(*,*,END=99) WAVE
      P = HPROFL(N, M, WAVE, WAVEH, T, XNE, H1, HE1, DOP)
      Q = HLINOP(WAVE, N, M, WAVEH, T, XNE, H1, HE1, DOP)
      W = HLCUT(WAVE, N, M, WAVEH)
      WRITE(*,'(1P,4E25.16)') WAVE, P, Q, W
      GOTO 10
   99 CONTINUE
      END
