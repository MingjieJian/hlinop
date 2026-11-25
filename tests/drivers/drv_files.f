C  Test driver: HLINPROF_SET_FILES, HLINPROF_SOURCE and 
C  HLINPROF_SET_MESSAGES (hlinprof.f).
C  stdin, line 1: path of the Stehle table, or - (keep DATA.FILES)
C         line 2: message mode (1 or 0; negative: do not call)
C         then one call per line:
C         NLOW NUP WAVE0 TEMP XNE H1FRC HE1FRC DOP WAVE
C  stdout: PROF (per A)  SOURCE
      PROGRAM DRVFIL
      INTEGER NLOW, NUP, HLINPROF_SOURCE, MSG
      REAL*8 WAVE, WAVE0
      REAL T, XNE, H1, HE1, DOP, PROF
      CHARACTER*592 PATH
      READ(*,'(A)') PATH
      READ(*,*) MSG
      IF (PATH.NE.'-') CALL HLINPROF_SET_FILES(PATH, ' ', ' ')
      IF (MSG.GE.0) CALL HLINPROF_SET_MESSAGES(MSG)
   10 READ(*,*,END=99) NLOW, NUP, WAVE0, T, XNE, H1, HE1, DOP, WAVE
      CALL HLINPROF(WAVE, WAVE0, T, XNE, NLOW, NUP, H1, HE1, DOP, PROF)
      WRITE(*,'(1P,E25.16,I4)') PROF, HLINPROF_SOURCE()
      GOTO 10
   99 CONTINUE
      END
