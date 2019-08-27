        MODULE nddo_restraints_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE nddo_restraints(my_d, my_array1, my_array2, &
          my_factor, my_symbol, &
          my_file)

        INTEGER :: my_d
        REAL(KIND=dp), DIMENSION(0:f, m, a, my_d) :: my_array1
        REAL(KIND=dp), DIMENSION(0:f, m, my_d) :: my_array2
        REAL(KIND=dp) :: my_factor
        CHARACTER(len=1) my_symbol
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, jd, ja

        WRITE(6,'(2a)') '# SUBROUTINE = nddo_restraints; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        my_array2 = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(2a)') &
          '#  1f  2tim  3m  ', &
          '4ptot  5f'

        DO jf = 0, f
          DO jm = 1, m

            DO jd = 1, my_d
              DO ja = 1, a
                IF ( symbol(ja) == my_symbol ) &
                  my_array2(jf,jm,jd) = &
                  my_array2(jf,jm,jd) + &
                  my_array1(jf,jm,ja,jd)*my_factor
              END DO
            END DO

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (my_array2(jf,jm,jd),jd=1,my_d)

          END DO
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,e20.10),1(1x,i4),7(1x,e20.10))

        END SUBROUTINE nddo_restraints

        END MODULE nddo_restraints_mod
