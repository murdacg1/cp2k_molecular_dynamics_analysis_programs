        MODULE reading_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading(my_d, total, wat, sol, &
          my_file)

        INTEGER :: my_d
        REAL(KIND=dp), DIMENSION(0:f,m,a,my_d) :: total
        REAL(KIND=dp), DIMENSION(0:f,m_wat,a_wat,my_d) :: wat
        REAL(KIND=dp), DIMENSION(0:f,m_sol,a_sol,my_d) :: sol
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja, jd
        INTEGER :: a_temp
        INTEGER :: jf_temp
        REAL(KIND=dp) :: tim_temp, pot_temp
        REAL(KIND=dp) :: mass_temp
        REAL(KIND=dp) :: shift_mag
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: middle_old, middle
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: shift
        CHARACTER(len=3) char3_temp
        CHARACTER(len=8) char8_temp
        CHARACTER(len=5) char5_temp
        CHARACTER(len=2) symbol_temp

        WRITE(6,'(2a)') '# SUBROUTINE = reading; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        tim = zero
        pot = zero
        symbol = '  '
        symbol_wat = '  '
        symbol_sol = '  '
        mass = zero
        mass_wat = zero
        mass_sol = zero
        total = zero
        wat = zero
        sol = zero

!       WRITE(6,*)
!       WRITE(6,*)
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jf = 0, f
          READ(3,*) a_temp
          READ(3,200) &
            char3_temp, jf_temp, &
            char8_temp, tim_temp, &
            char5_temp, pot_temp
          tim(jf) = tim_temp*fs2ps
          pot(jf) = pot_temp
!         WRITE(6,200) &
!           char3_temp, jf_temp, &
!           char8_temp, tim(jf)/fs2ps, &
!           char5_temp, pot(jf)
!         CALL FLUSH(6)
! total
          DO jm = 1, m
            DO ja = 1, a
              READ(3,300)  symbol_temp, (total(jf,jm,ja,jd),jd=1,my_d)
              SELECT CASE (symbol_temp)
              CASE ('Si')
                mass_temp = 27.9769271_dp
              CASE (' O')
                mass_temp = 15.99491463_dp
              CASE (' H')
                mass_temp = 1.0078250321_dp
              CASE ('Cl')
                mass_temp = 34.968852721_dp
              END SELECT
              symbol(ja) = symbol_temp
              mass(ja) = mass_temp
!             WRITE(6,300) symbol(ja), (total(jf,jm,ja,jd),jd=1,my_d)
!             CALL FLUSH(6)
            END DO
          END DO
        END DO
        CLOSE(3)
! total

!       WRITE(6,*)
!       WRITE(6,*)
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jf = 0, f
          READ(3,*) a_temp
          READ(3,200) &
            char3_temp, jf_temp, &
            char8_temp, tim_temp, &
            char5_temp, pot_temp
          tim(jf) = tim_temp*fs2ps
          pot(jf) = pot_temp
!         WRITE(6,200) &
!           char3_temp, jf_temp, &
!           char8_temp, tim(jf)/fs2ps, &
!           char5_temp, pot(jf)
!         CALL FLUSH(6)
! wat
!         WRITE(6,*)
!         WRITE(6,*) 'wat'
          DO jm = 1, m_wat
            DO ja = 1, a_wat
              READ(3,300)  symbol_temp, (wat(jf,jm,ja,jd),jd=1,my_d)
              SELECT CASE (symbol_temp)
              CASE ('Si')
                mass_temp = 27.9769271_dp
              CASE (' O')
                mass_temp = 15.99491463_dp
              CASE (' H')
                mass_temp = 1.0078250321_dp
              CASE ('Cl')
                mass_temp = 34.968852721_dp
              END SELECT
              symbol_wat(ja) = symbol_temp
              mass_wat(ja) = mass_temp
!             WRITE(6,300) symbol_wat(ja), (wat(jf,jm,ja,jd),jd=1,my_d)
!             CALL FLUSH(6)
            END DO
          END DO
! wat
! sol
!         WRITE(6,200)
!         WRITE(6,200) 'sol'
          DO jm = 1, m_sol
            DO ja = 1, a_sol
              READ(3,300)  symbol_temp, (sol(jf,jm,ja,jd),jd=1,my_d)
              SELECT CASE (symbol_temp)
              CASE ('Si')
                mass_temp = 27.9769271_dp
              CASE (' O')
                mass_temp = 15.99491463_dp
              CASE (' H')
                mass_temp = 1.0078250321_dp
              CASE ('Cl')
                mass_temp = 34.968852721_dp
              END SELECT
              symbol_sol(ja) = symbol_temp
              mass_sol(ja) = mass_temp
!             WRITE(6,300) symbol_sol(ja), (sol(jf,jm,ja,jd),jd=1,my_d)
!             CALL FLUSH(6)
            END DO
          END DO
! sol
        END DO
        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

200     FORMAT(1x, &
          a3,1x,i8, &
          a8,1x,f12.3, &
          a5,1x,f20.10 &
          )
300     FORMAT(1x,a,1x,3(f20.10))
500     FORMAT(1(1x,i10),1(1x,f20.10),4(1x,f20.10))

        END SUBROUTINE reading

        END MODULE reading_mod
