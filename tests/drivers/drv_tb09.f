C  Test driver: HLINOP_SET_TB09 and the Stehle/TB09 blending in STEHL2
C  (hlinprof.f).
C  stdin, line 1: path of the Stehle table (argument HFILE of STEHL2)
C         then commands, one per line:
C         1               next line: MODE DLO DHI; next line: TB09 path
C                         -> CALL HLINOP_SET_TB09
C         2 NUP WAVE0 T XNE DLAM IBOTH
C                         -> STEHL2 and STEHLS at WAVE0 + DLAM (Balmer)
C            stdout: HLIN HLINM (STEHL2)  HLIN HLINM (STEHLS)  ANESL
C         3 NUP WAVE0 T XNE DOP DLAM
C                         -> HLINPROF (n_H = n_He = 0) at WAVE0 + DLAM
C            stdout: PROF
      PROGRAM DRVTB9
      INTEGER ICMD, MODE, NUP, IBOTH
      REAL*8 WAVE0, DLAM
      REAL DLO, DHI, T, XNE, DOP, H, HM, H0, HM0, PROF
      CHARACTER*592 SFILE, TFILE
      CHARACTER*256 LINE
      COMMON /HLSLIM/ ANESL
      READ(*,'(A)') SFILE
   10 READ(*,'(A)',END=99) LINE
      IF (LEN_TRIM(LINE).EQ.0) GOTO 10
      READ(LINE,*) ICMD
      IF (ICMD.EQ.1) THEN
        READ(*,*) MODE, DLO, DHI
        READ(*,'(A)') TFILE
        CALL HLINOP_SET_TB09(TFILE(1:LEN_TRIM(TFILE)), MODE, DLO, DHI)
      ELSE IF (ICMD.EQ.2) THEN
        READ(LINE,*) ICMD, NUP, WAVE0, T, XNE, DLAM, IBOTH
        HM = 0.
        HM0 = 0.
        CALL STEHL2(SFILE(1:LEN_TRIM(SFILE)), WAVE0 + DLAM, WAVE0, T,
     *              XNE, 2, NUP, H, HM, IBOTH)
        CALL STEHLS(SFILE(1:LEN_TRIM(SFILE)), WAVE0 + DLAM, WAVE0, T,
     *              XNE, 2, NUP, H0, HM0, IBOTH)
        WRITE(*,'(1P,5E25.16)') H, HM, H0, HM0, ANESL
      ELSE IF (ICMD.EQ.3) THEN
        READ(LINE,*) ICMD, NUP, WAVE0, T, XNE, DOP, DLAM
        CALL HLINPROF(WAVE0 + DLAM, WAVE0, T, XNE, 2, NUP, 0., 0., DOP,
     *                PROF)
        WRITE(*,'(1P,E25.16)') PROF
      END IF
      GOTO 10
   99 CONTINUE
      END
