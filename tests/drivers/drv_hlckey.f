C  Test driver: the cache index HLCKEY (hlinop.f) with NSLOT = 16 hash
C  slots, NENT = 12 entries and one-number keys.
C  stdin: one key per line; a negative key calls HLINOP_CLEAR_CACHE
C         (no output line)
C  stdout: IE HIT NUSED
      PROGRAM DRVHCK
      INTEGER NSLOT, NENT
      PARAMETER (NSLOT=16, NENT=12)
      INTEGER NUSED, IGEN, IE, IHIT, IDX(NSLOT)
      REAL*8 KEY(1), KEYS(1,NENT)
      LOGICAL HIT
      DATA IDX/NSLOT*0/, NUSED/0/, IGEN/0/
   10 READ(*,*,END=99) KEY(1)
      IF (KEY(1).LT.0.D0) THEN
        CALL HLINOP_CLEAR_CACHE
        GOTO 10
      ENDIF
      CALL HLCKEY(KEY, 1, KEYS, NSLOT, IDX, NENT, NUSED, IGEN, IE, HIT)
      IHIT = 0
      IF (HIT) IHIT = 1
      WRITE(*,'(3I6)') IE, IHIT, NUSED
      GOTO 10
   99 CONTINUE
      END
