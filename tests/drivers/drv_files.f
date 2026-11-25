C  Test driver: HLINPROF_SET_FILES and HLINPROF_SOURCE (hlinprof.f).
C  stdin, line 1: path of the Stehle table, or - (keep DATA.FILES)
C         then one call per line:
C         NLOW NUP WAVE0 TEMP XNE H1FRC HE1FRC DOP WAVE
C  stdout: PROF (per A)  SOURCE
      PROGRAM DRVFIL
      INTEGER NLOW, NUP, HLINPROF_SOURCE
      REAL*8 WAVE, WAVE0
      REAL T, XNE, H1, HE1, DOP, PROF
      CHARACTER*592 PATH
      READ(*,'(A)') PATH
      IF (PATH.NE.'-') CALL HLINPROF_SET_FILES(PATH, ' ', ' ')
   10 READ(*,*,END=99) NLOW, NUP, WAVE0, T, XNE, H1, HE1, DOP, WAVE
      CALL HLINPROF(WAVE, WAVE0, T, XNE, NLOW, NUP, H1, HE1, DOP, PROF)
      WRITE(*,'(1P,E25.16,I4)') PROF, HLINPROF_SOURCE()
      GOTO 10
   99 CONTINUE
      END
