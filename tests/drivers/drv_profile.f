C  Test driver: HPROFL and HLINOP profiles (hlinop.f).
C  stdin, line 1:  N M WAVEH T XNE H1FRC HE1FRC DOP
C         then one wavelength (A) per line
C  stdout: WAVE HPROFL HLINOP   (both per unit frequency, Hz^-1)
      PROGRAM DRVPROF
      INTEGER N, M
      REAL*8 WAVE, WAVEH
      REAL T, XNE, H1, HE1, DOP, HPROFL, HLINOP, P, Q
      READ(*,*) N, M, WAVEH, T, XNE, H1, HE1, DOP
   10 READ(*,*,END=99) WAVE
      P = HPROFL(N, M, WAVE, WAVEH, T, XNE, H1, HE1, DOP)
      Q = HLINOP(WAVE, N, M, WAVEH, T, XNE, H1, HE1, DOP)
      WRITE(*,'(1P,3E25.16)') WAVE, P, Q
      GOTO 10
   99 CONTINUE
      END
