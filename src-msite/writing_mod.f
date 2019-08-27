        MODULE writing_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE writing(my_d, my_array, my_file)

        INTEGER :: my_d
        REAL(KIND=dp), DIMENSION(0:f,m,a+1,my_d) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja

        WRITE(6,'(2a)') '# SUBROUTINE = writing; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        OPEN(3,file=TRIM(my_file),status='unknown')

        DO jf = 0, f
!         WRITE(6,*) 'jf = ', jf
!         CALL FLUSH(6)
          WRITE(3,*) m*(a+1)
          WRITE(3,*)
          DO jm = 1, m
            DO ja = 1, a+1
              WRITE(3,100) &
                symbol(jm,ja), my_array(jf,jm,ja,:), molecule(jm,ja)
            END DO
          END DO
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(2x,a,1x,3(f20.10),1x,a)

        END SUBROUTINE writing

        END MODULE writing_mod
