        MODULE reading_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading(my_file)

        CHARACTER(len=80) my_file
 
        INTEGER :: jf, jm, ja
        INTEGER :: a_max_temp
        CHARACTER(len=80) :: string80

        WRITE(6,'(2a)') '# SUBROUTINE = reading; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        coor = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        DO jf = 0, f
!         WRITE(6,*) 'jf = ', jf
!         CALL FLUSH(6)
          READ(3,*) a_max_temp
          READ(3,*) string80
          DO jm = 1, m
            DO ja = 1, a
              READ(3,*) &
                symbol(jm,ja), coor(jf,jm,ja,:)
            END DO
          END DO
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

120     FORMAT(1x, &
          a3,1x,i8, &
          a8,1x,f12.3, &
          a5,1x,f20.10 &
          )
220     FORMAT(2x,a,3(f20.10))

        END SUBROUTINE reading

        END MODULE reading_mod
