C  Test driver: occupation probability WCALC (hbop.f).
C  stdin, one case per line:  NH NE NHE NS T SCALE
C         SCALE >= 0: CALL HBOP_SET_EXCLUDED_VOLUME(SCALE) first;
C         SCALE < 0: leave the setting as it is (default on start-up)
C  stdout: WCALC
      PROGRAM DRVWCALC
      REAL NH, NE, NHE, NS, T, SCALE
      REAL*8 WCALC
   10 READ(*,*,END=99) NH, NE, NHE, NS, T, SCALE
      IF (SCALE.GE.0.) CALL HBOP_SET_EXCLUDED_VOLUME(SCALE)
      WRITE(*,'(1P,E25.16)') WCALC(NH, NE, NHE, NS, T)
      GOTO 10
   99 CONTINUE
      END
