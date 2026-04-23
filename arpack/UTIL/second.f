      REAL FUNCTION SECOND()
      IMPLICIT NONE
      REAL :: T

      CALL CPU_TIME(T)

      IF (T < 0.0) THEN
         SECOND = 0.0          ! fallback if cpu_time fails
      ELSE
         SECOND = T
      END IF

      RETURN
      END
