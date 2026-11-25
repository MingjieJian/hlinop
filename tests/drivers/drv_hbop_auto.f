C  Test driver: HBOP_AUTO (hbop_auto.f).
C  stdin, line 1:  NMAX NH NHE NE T DOP; then one wavelength (A, air) per line
C  stdout: WAVE TOTAL CONTIN (cm^-1)
      PROGRAM DRVAUTO
      INTEGER NMAX
      REAL*8 WAVE
      REAL NH, NHE, NE, T, DOP, TOTAL, CONTIN
      READ(*,*) NMAX, NH, NHE, NE, T, DOP
   10 READ(*,*,END=99) WAVE
      CALL HBOP_AUTO(WAVE, NMAX, NH, NHE, NE, T, DOP, TOTAL, CONTIN)
      WRITE(*,'(1P,3E25.16)') WAVE, TOTAL, CONTIN
      GOTO 10
   99 CONTINUE
      END
