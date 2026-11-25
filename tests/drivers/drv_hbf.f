C  Test driver: bound-free cross section HBF (hbop.f).
C  stdin, one case per line:  WAVE NLO
C  stdout: WAVE NLO HBF   (cm^2)
      PROGRAM DRVHBF
      INTEGER NLO
      REAL*8 WAVE
      REAL HBF
   10 READ(*,*,END=99) WAVE, NLO
      WRITE(*,'(1P,E25.16,I4,E25.16)') WAVE, NLO, HBF(WAVE, NLO)
      GOTO 10
   99 CONTINUE
      END
