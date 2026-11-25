C  Test driver: HLINPROF (hlinprof.f).
C  stdin, line 1:  NLOW NUP WAVE0 TEMP XNE H1FRC HE1FRC DOP
C         then one wavelength (A) per line
C  stdout: WAVE PROF   (area normalised per Angstrom)
      PROGRAM DRVHLP
      INTEGER NLOW, NUP
      REAL*8 WAVE, WAVE0
      REAL T, XNE, H1, HE1, DOP, PROF
      READ(*,*) NLOW, NUP, WAVE0, T, XNE, H1, HE1, DOP
   10 READ(*,*,END=99) WAVE
      CALL HLINPROF(WAVE, WAVE0, T, XNE, NLOW, NUP, H1, HE1, DOP, PROF)
      WRITE(*,'(1P,2E25.16)') WAVE, PROF
      GOTO 10
   99 CONTINUE
      END
