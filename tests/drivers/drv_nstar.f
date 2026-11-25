C  Test driver: HBNSTAR (hbop.f) and the air/vacuum conversions (hlinop.f).
C  stdin, one case per line:  WAVE NLO
C  stdout: WAVE AIRVAC(WAVE) VACAIR(WAVE) HBNSTAR(WAVE,NLO) HBNSTAR(AIRVAC(WAVE),NLO)
      PROGRAM DRVNSTAR
      INTEGER NLO
      REAL*8 WAVE, AIRVAC, VACAIR, HBNSTAR
   10 READ(*,*,END=99) WAVE, NLO
      WRITE(*,'(1P,5E25.16)') WAVE, AIRVAC(WAVE), VACAIR(WAVE),
     *    HBNSTAR(WAVE, NLO), HBNSTAR(AIRVAC(WAVE), NLO)
      GOTO 10
   99 CONTINUE
      END
