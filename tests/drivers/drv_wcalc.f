C  Test driver: occupation probability WCALC (hbop.f).
C  stdin, one case per line:  NH NE NHE NS T
C  stdout: WCALC
      PROGRAM DRVWCALC
      REAL NH, NE, NHE, NS, T
      REAL*8 WCALC
   10 READ(*,*,END=99) NH, NE, NHE, NS, T
      WRITE(*,'(1P,E25.16)') WCALC(NH, NE, NHE, NS, T)
      GOTO 10
   99 CONTINUE
      END
