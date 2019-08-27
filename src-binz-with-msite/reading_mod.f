        MODULE reading_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading(my_d, my_array, &
          my_file)

        INTEGER :: my_d
        REAL(KIND=dp), DIMENSION(0:f) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja
        INTEGER :: a_max_temp
        REAL(KIND=dp) :: temp
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: coor
        CHARACTER(len=80) :: string80
        CHARACTER(len=1) symbol

        WRITE(6,'(2a)') '# SUBROUTINE = reading; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        my_array = zero

        ALLOCATE ( coor(d) )
        coor = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        DO jf = 0, f
!         WRITE(6,*) 'jf = ', jf
!         CALL FLUSH(6)
          READ(3,*) a_max_temp
!         WRITE(6,*) a_max_temp
          CALL FLUSH(6)
          READ(3,*) string80
!         WRITE(6,*) string80
!         CALL FLUSH(6)
          DO jm = 1, m_w
            DO ja = 1, a_w
              READ(3,*) &
                symbol, coor(:)
            END DO
          END DO
!         DO jm = 1, m_s ! m_s is always 1
            DO ja = 1, a_s
              READ(3,*) &
                symbol, coor(:)
                IF (ja == which) &
                  my_array(jf) = coor(my_d) - middle
            END DO
!         END DO
        END DO

        CLOSE(3)

        WRITE(6,*) 'my_array(0) before = ', my_array(0) ! handle the h2o z=-1 instead of z=1 centered problem
        IF ( ABS(-my_array(0)-z_s) <  ABS(my_array(0)-z_s) ) &
          my_array = -my_array
        WRITE(6,*) 'my_array(0) after = ', my_array(0)

        DEALLOCATE ( coor )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        END SUBROUTINE reading

        END MODULE reading_mod
