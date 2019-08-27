        MODULE reading_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading(my_d, my_array, my_file)

        INTEGER :: my_d
        REAL(KIND=dp), DIMENSION(0:f,m,a+1,my_d) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja
        INTEGER :: a_max_temp

        WRITE(6,'(2a)') '# SUBROUTINE = reading; filename = ', &
          TRIM(my_file)
        WRITE(6,*) 'f, m, a, my_d = ', f, m, a, my_d
        CALL FLUSH(6)

        my_array = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        DO jf = 0, f
          READ(3,*) a_max_temp
          READ(3,*)
!         WRITE(7,*) a_max_temp
!         WRITE(7,*)
          CALL FLUSH(7)
          DO jm = 1, m
            DO ja = 1, a
              READ(3,*) &
                symbol(jm,ja), my_array(jf,jm,ja,:), molecule(jm,ja)
!             WRITE(7,100) &
!               symbol(jm,ja), my_array(jf,jm,ja,:), molecule(jm,ja)
!             CALL FLUSH(7)
            END DO
            symbol(jm,4) = 'M'
            molecule(jm,4) = 'H2O'
          END DO
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(2x,a,1x,3(f20.10),1x,a)

        END SUBROUTINE reading

        END MODULE reading_mod
