C  Test driver: HLINPROF and HLINOP for a sequence of states (hlinop.f,
C  hlinprof.f), in the order given, in one process (cache behaviour).
C  stdin, one call per line:
C         NLOW NUP WAVE0 TEMP XNE H1FRC HE1FRC DOP WAVE
C         (NLOW = 0: CALL HLINOP_CLEAR_CACHE, no output line)
C  stdout: PROF (HLINPROF, per A)  HL (HLINOP, per Hz)
      PROGRAM DRVSTA
      INTEGER NLOW, NUP
      REAL*8 WAVE, WAVE0
      REAL T, XNE, H1, HE1, DOP, PROF, HL, HLINOP
   10 READ(*,*,END=99) NLOW, NUP, WAVE0, T, XNE, H1, HE1, DOP, WAVE
      IF (NLOW.EQ.0) THEN
        CALL HLINOP_CLEAR_CACHE
        GOTO 10
      ENDIF
      CALL HLINPROF(WAVE, WAVE0, T, XNE, NLOW, NUP, H1, HE1, DOP, PROF)
      HL = HLINOP(WAVE, NLOW, NUP, WAVE0, T, XNE, H1, HE1, DOP)
      WRITE(*,'(1P,2E25.16)') PROF, HL
      GOTO 10
   99 CONTINUE
      END
