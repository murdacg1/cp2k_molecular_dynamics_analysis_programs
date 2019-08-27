        MODULE readit_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE readit(my_t, my_titles, my_m, my_n, my_array, my_file)

        INTEGER :: my_t, my_m, my_n
        CHARACTER(len=80), DIMENSION(t) ::  my_titles
        REAL(KIND=dp), DIMENSION(m, n) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: jt, jm

        WRITE(6,'(2a)') '# SUBROUTINE = readit; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        my_titles = ' '
        my_array = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

!       WRITE(6,*)
!       WRITE(6,*)

        DO jt = 1, my_t
          READ(3,*) my_titles(jt)
!         WRITE(6,*) 'my_titles(jt) = ', TRIM(my_titles(jt))
!         CALL FLUSH(6)
        END DO

        DO jm = 1, my_m
          READ(3,*) my_array(jm,:)
!         WRITE(6,*) 'jm, my_array(jm,1:my_n) = ', jm, my_array(jm,1:my_n)
!         CALL FLUSH(6)
        END DO
        CLOSE(3)

        END SUBROUTINE readit

        END MODULE readit_mod
