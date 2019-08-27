        MODULE reading_and_calc_ke_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading_and_calc_ke(total, wat, sol, &
          kin, kin_slab, kin_wat, kin_sol, kin_com_sol, kin_vibrot_sol, uu_com_sol,  &
          my_file)

        REAL(KIND=dp), DIMENSION(0:f,a,d3) :: total
        REAL(KIND=dp), DIMENSION(0:f,m_wat,a_wat,d3) :: wat
        REAL(KIND=dp), DIMENSION(0:f,m_sol,a_sol,d3) :: sol
        CHARACTER(len=80) my_file

        REAL(KIND=dp), DIMENSION(0:f,1,d1) :: kin
        REAL(KIND=dp), DIMENSION(0:f,1,d1) :: kin_slab
        REAL(KIND=dp), DIMENSION(0:f,m_wat,d1) :: kin_wat
        REAL(KIND=dp), DIMENSION(0:f,m_sol,d1) :: kin_sol
        REAL(KIND=dp), DIMENSION(0:f,m_sol,d4) :: kin_com_sol
        REAL(KIND=dp), DIMENSION(0:f,m_sol,d1) :: kin_vibrot_sol
        REAL(KIND=dp), DIMENSION(0:f,m_sol,d4) :: uu_com_sol

        INTEGER :: jf, ja, jd
        INTEGER :: jm
        INTEGER :: a_temp
        INTEGER :: jf_temp
        REAL(KIND=dp) :: tim_temp, pot_temp
        REAL(KIND=dp) :: mass_temp, uu_com_sol_sq, uu_com_sol_comp_sq
        CHARACTER(len=3) char3_temp
        CHARACTER(len=8) char8_temp
        CHARACTER(len=5) char5_temp
        CHARACTER(len=2) symbol_temp
        REAL(KIND=dp), DIMENSION(d3) :: temp_array = zero

        REAL(KIND=dp) :: mass_sol_tot
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: v_com_sol

        WRITE(6,'(2a)') &
          '# SUBROUTINE = reading_and_calc_ke; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        tim = zero
        symbol = '  '
        symbol_wat = '  '
        symbol_sol = '  '
        mass = zero
        mass_wat = zero
        mass_sol = zero
        total = zero
        wat = zero
        sol = zero

        kin = zero
        kin_slab = zero
        kin_wat = zero
        kin_sol = zero
        kin_com_sol = zero
        kin_vibrot_sol = zero
        uu_com_sol = zero

        ALLOCATE( v_com_sol(0:f,m_sol,d3) )
        v_com_sol = zero

!!! first reading to fill some arrays for symbol and masses, no calculations here !!!
!
!       WRITE(6,*)
!       WRITE(6,*)
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jf = 0, 0
          READ(3,*) a_temp
          READ(3,200) &
            char3_temp, jf_temp, &
            char8_temp, tim_temp, &
            char5_temp, pot_temp
          tim(jf) = tim_temp*fs2ps
!         WRITE(6,200) &
!           char3_temp, jf_temp, &
!           char8_temp, tim(jf), &
!           char5_temp, pot_temp
!         CALL FLUSH(6)
! total
          DO ja = 1, a
            READ(3,300)  symbol_temp, (temp_array(jd),jd=1,d3)
              SELECT CASE (symbol_temp)
              CASE (' O')
                mass_temp = mass_O
              CASE (' H')
                mass_temp = mass_H
              CASE (' C')
                mass_temp = mass_C
              CASE ('Cl')
                mass_temp = mas_Cl
              CASE (' N')
                mass_temp = mass_N
              CASE (' M')
                mass_temp = mass_M
              CASE ('M2')
                mass_temp = mass_M
              CASE ('N2')
                mass_temp = mass_N
              CASE (' S')
                mass_temp = mass_S
              CASE ('Si')
                mass_temp = mass_Si
              END SELECT
            symbol(ja) = symbol_temp
            mass(ja) = mass_temp * m_u2m_amu
!           WRITE(6,300) symbol(ja), (temp_array(jd),jd=1,d3)
!           CALL FLUSH(6)
          END DO
        END DO
        CLOSE(3)
! total
!
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
!         WRITE(6,200) &
!           char3_temp, jf_temp, &
!           char8_temp, tim(jf), &
!           char5_temp, pot_temp
!         CALL FLUSH(6)
! wat
!         WRITE(6,*)
!         WRITE(6,*) 'wat'
          DO jm = 1, m_wat
            DO ja = 1, a_wat
              READ(3,300)  symbol_temp, (temp_array(jd),jd=1,d3)
              SELECT CASE (symbol_temp)
              CASE (' O')
                mass_temp = mass_O
              CASE (' H')
                mass_temp = mass_H
              CASE (' C')
                mass_temp = mass_C
              CASE ('Cl')
                mass_temp = mas_Cl
              CASE (' N')
                mass_temp = mass_N
              CASE (' M')
                mass_temp = mass_M
              CASE ('M2')
                mass_temp = mass_M
              CASE ('N2')
                mass_temp = mass_N
              CASE (' S')
                mass_temp = mass_S
              CASE ('Si')
                mass_temp = mass_Si
              END SELECT
              symbol_wat(ja) = symbol_temp
              mass_wat(ja) = mass_temp * m_u2m_amu
!             WRITE(6,300) symbol_wat(ja), (temp_array(jd),jd=1,d3)
!             CALL FLUSH(6)
            END DO
          END DO
! wat
!
! sol
!         WRITE(6,200)
!         WRITE(6,200) 'sol'
          DO jm = 1, m_sol
            DO ja = 1, a_sol
              READ(3,300)  symbol_temp, (temp_array(jd),jd=1,d3)
              SELECT CASE (symbol_temp)
              CASE (' O')
                mass_temp = mass_O
              CASE (' H')
                mass_temp = mass_H
              CASE (' C')
                mass_temp = mass_C
              CASE ('Cl')
                mass_temp = mas_Cl
              CASE (' N')
                mass_temp = mass_N
              CASE (' M')
                mass_temp = mass_M
              CASE ('M2')
                mass_temp = mass_M
              CASE ('N2')
                mass_temp = mass_N
              CASE (' S')
                mass_temp = mass_S
              CASE ('Si')
                mass_temp = mass_Si
              END SELECT
              symbol_sol(ja) = symbol_temp
              mass_sol(ja) = mass_temp * m_u2m_amu
!             WRITE(6,300) symbol_sol(ja), (temp_array(jd),jd=1,d3)
!             CALL FLUSH(6)
            END DO
          END DO
! sol
        END DO
        CLOSE(3)
        mass_sol_tot = SUM( mass_sol )
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

!!! actual reading and calculations !!!
!
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jf = 0, f
          READ(3,*) a_temp
          READ(3,200) &
            char3_temp, jf_temp, &
            char8_temp, tim_temp, &
            char5_temp, pot_temp
          tim(jf) = tim_temp*fs2ps
! total
          DO ja = 1, a
            READ(3,300)  symbol_temp, (total(jf,ja,jd),jd=1,d3)
         END DO
        END DO
        CLOSE(3)
! total
!
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jf = 0, f
          READ(3,*) a_temp
          READ(3,200) &
            char3_temp, jf_temp, &
            char8_temp, tim_temp, &
            char5_temp, pot_temp
          tim(jf) = tim_temp*fs2ps
! wat
          DO jm = 1, m_wat
            DO ja = 1, a_wat
              READ(3,300)  symbol_temp, (wat(jf,jm,ja,jd),jd=1,d3)
              kin_wat(jf,jm,d1) = kin_wat(jf,jm,d1) + &
                0.5_dp * mass_wat(ja) * DOT_PRODUCT( wat(jf,jm,ja,:), wat(jf,jm,ja,:) )
            END DO
          END DO
! wat
!
! sol
          DO jm = 1, m_sol
            DO ja = 1, a_sol
              READ(3,300)  symbol_temp, (sol(jf,jm,ja,jd),jd=1,d3)
              kin_sol(jf,jm,d1) = kin_sol(jf,jm,d1) + &
                0.5_dp * mass_sol(ja) * DOT_PRODUCT( sol(jf,jm,ja,:), sol(jf,jm,ja,:) )
              v_com_sol(jf,jm,:) = v_com_sol(jf,jm,:) + mass_sol(ja) * sol(jf,jm,ja,:)
            END DO
            v_com_sol(jf,jm,:) = v_com_sol(jf,jm,:) / mass_sol_tot
            uu_com_sol_sq = DOT_PRODUCT( v_com_sol(jf,jm,:), v_com_sol(jf,jm,:) )
            kin_com_sol(jf,jm,d4) = 0.5_dp * mass_sol_tot * uu_com_sol_sq
            kin_vibrot_sol(jf,jm,d1) = kin_sol(jf,jm,d1) - kin_com_sol(jf,jm,d4)
            uu_com_sol(jf,jm,d4) = SQRT( uu_com_sol_sq )
            DO jd = 1, d3
              uu_com_sol_comp_sq = ( v_com_sol(jf,jm,jd) )**2
              kin_com_sol(jf,jm,jd) = 0.5_dp * mass_sol_tot * uu_com_sol_comp_sq
              uu_com_sol(jf,jm,jd) = SQRT( uu_com_sol_comp_sq )
            END DO
          END DO
! sol
!
! slab and whole cell 
          kin_slab(jf,1,d1) = SUM( kin_wat(jf,1:m_wat,d1) )
          kin(jf,1,d1) = kin_slab(jf,1,d1) + SUM( kin_sol(jf,1:m_sol,d1) )
! slab and whole cell 
!
        END DO
        CLOSE(3)

        DEALLOCATE( v_com_sol )

        u_com_sol = u_com_sol * au_velocity_to_ang_per_ps

        my_file='ke_wat.dat'
        OPEN(11,file=TRIM(my_file),status='unknown')
        WRITE(11,'(4a)') &
          '#          1f         2tim   3m', &
          '               4ke_wat'
        DO jf = 0, f
          DO jm = 1, m_wat
            WRITE(11,101) jf, tim(jf), jm, &
              kin_wat(jf,jm,d1)
          END DO
        END DO
        CALL FLUSH(11)
        CLOSE(11)

        my_file='ke_sol.dat'
        OPEN(12,file=TRIM(my_file),status='unknown')
        WRITE(12,'(7a)') &
          '#          1f         2tim   3m', &
          '           4ke_com_sol', &
          '        5ke_vibrot_sol', &
          '               6ke_sol', &
          '         7ke_com_sol_x', &
          '         8ke_com_sol_y', &
          '         9ke_com_sol_z'
        my_file='u_com_sol.dat'
        OPEN(13,file=TRIM(my_file),status='unknown')
        WRITE(13,'(5a)') &
          '#          1f         2tim   3m', &
          '            4u_com_sol', &
          '          5u_com_sol_x', &
          '          6u_com_sol_y', &
          '          7u_com_sol_z'
        DO jf = 0, f
          DO jm = 1, m_sol
            WRITE(12,106) jf, tim(jf), jm, &
              kin_com_sol(jf,jm,d4), &
              kin_vibrot_sol(jf,jm,d1), &
              kin_sol(jf,jm,d1), &
              kin_com_sol(jf,jm,1:d3)
            WRITE(13,104) jf, tim(jf), jm, &
              u_com_sol(jf,jm,d4), &
              u_com_sol(jf,jm,1:d3)
          END DO
        END DO
        CALL FLUSH(12)
        CLOSE(12)
        CALL FLUSH(13)
        CLOSE(13)

        my_file='ke.dat'
        OPEN(13,file=TRIM(my_file),status='unknown')
        WRITE(13,'(4a)') &
          '#          1f         2tim   3m', &
          '              4ke_slab', &
          '               5ke_sol', &
          '              6ke_cell'
        DO jf = 0, f
          DO jm = 1, 1
            WRITE(13,103) jf, tim(jf), jm, &
              kin_slab(jf,1,d1), &
              SUM( kin_sol(jf,1:m_sol,d1) ), &
              kin(jf,1,d1)
          END DO
        END DO
        CALL FLUSH(13)
        CLOSE(13)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

101     FORMAT(1(1x,i12),1(1x,f12.6),1(1x,i4),1(1x,f21.12))
103     FORMAT(1(1x,i12),1(1x,f12.6),1(1x,i4),3(1x,f21.12))
104     FORMAT(1(1x,i12),1(1x,f12.6),1(1x,i4),4(1x,f21.12))
106     FORMAT(1(1x,i12),1(1x,f12.6),1(1x,i4),6(1x,f21.12))
200     FORMAT(1x, &
          a3,1x,i8, &
          a8,1x,f12.3, &
          a5,1x,f20.10 &
          )
300     FORMAT(1x,a,1x,3(f20.10))

        END SUBROUTINE reading_and_calc_ke

        END MODULE reading_and_calc_ke_mod
