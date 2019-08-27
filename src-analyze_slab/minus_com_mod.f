        MODULE minus_com_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE minus_com(my_array, my_com, &
          my_file)

        REAL(KIND=dp), DIMENSION(0:f,m,a,d) :: my_array
        REAL(KIND=dp), DIMENSION(0:f,d+1) :: my_com
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja, jd
        REAL(KIND=dp) :: my_mass_total
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: my_mass

        WRITE(6,'(2a)') '# SUBROUTINE = minus_com; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        my_mass_total = REAL(m,dp)*SUM(mass)
        ALLOCATE( my_mass(a) )
        my_mass = zero
        DO ja = 1, a
          my_mass(ja) = mass(ja)/my_mass_total
        END DO

        WRITE(6,*) 'mass = ', mass(:)
        WRITE(6,*) 'my_mass = ', my_mass(:)

        OPEN(3,file=TRIM(my_file),status='unknown')

        my_com = zero
        DO jf = 0, f
          DO jd = 1, d
            DO jm = 1, m
              DO ja = 1, a
                my_com(jf,jd) = &
                  my_com(jf,jd) + my_mass(ja)*my_array(jf,jm,ja,jd)
              END DO
            END DO
            my_com(jf,d+1) = &
              my_com(jf,d+1) + my_com(jf,jd)**2
          END DO
          my_com(jf,d+1) = &
            DSQRT( my_com(jf,d+1) )
        END DO

        WRITE(3,'(a)') '# jf, tim(jf), my_com(jf,:)'
        DO jf = 0, f
          WRITE(3,500) jf, tim(jf), my_com(jf,:)
        END DO

        DO jf = 0, f
          DO jm = 1, m
            DO ja = 1, a
              DO jd = 1, d
                my_array(jf,jm,ja,jd) = &
                  my_array(jf,jm,ja,jd) - my_com(jf,jd)
              END DO
            END DO
          END DO
        END DO

        WRITE(3,*)
        WRITE(3,*)

        my_com = zero
        DO jf = 0, f
          DO jd = 1, d
            DO jm = 1, m
              DO ja = 1, a
                my_com(jf,jd) = &
                  my_com(jf,jd) + my_mass(ja)*my_array(jf,jm,ja,jd)
              END DO
            END DO
            my_com(jf,d+1) = &
              my_com(jf,d+1) + my_com(jf,jd)**2
          END DO
          my_com(jf,d+1) = &
            DSQRT( my_com(jf,d+1) )
        END DO

        WRITE(3,'(a)') '# jf, tim(jf), my_com(jf,:)'
        DO jf = 0, f
          WRITE(3,500) jf, tim(jf), my_com(jf,:)
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( my_mass )

500     FORMAT(1(1x,i10),1(1x,f20.10),4(1x,f20.10))

        END SUBROUTINE minus_com

        END MODULE minus_com_mod
