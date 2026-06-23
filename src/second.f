      SUBROUTINE SECONDS(TSEC)
      REAL*8 TSEC
C
C...Returns elapsed wall-clock time in seconds
C   Replacement for non-standard SECNDS intrinsic
C
      INTEGER COUNT, COUNT_RATE
      CALL SYSTEM_CLOCK(COUNT, COUNT_RATE)
      TSEC = DBLE(COUNT) / DBLE(COUNT_RATE)
      END