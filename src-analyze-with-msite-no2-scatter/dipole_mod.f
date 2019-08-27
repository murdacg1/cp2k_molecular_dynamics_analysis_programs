        MODULE dipole_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE dipole(my_d_chrg, my_d_dipl, my_m_dipl_cell, &
          my_file)

        INTEGER :: my_d_chrg, my_d_dipl, my_m_dipl_cell
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja, jd, f_counter
        REAL(KIND=dp) :: dist, dist_cell, &
          volume, ms_dipl_cell, rms_dipl_cell, &
          dielectric_const
        REAL(KIND=dp) :: X, Y, Z
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: dipl_cell_temp
        REAL(KIND=dp) :: pi

        WRITE(6,'(2a)') '# SUBROUTINE = dipole; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        X = zero
        Y = zero
        Z = zero

        ALLOCATE( dipl_cell_temp(d) )

        dipl_cell_temp = zero

        dipl_mol = zero
        dipl_sys = zero
        dipl_cell = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(2a)') &
          '#  1f  2tim  3m  ', &
          '4dipx  5dipy  6dipz  7dip'

        DO jf = f0+1, f

          DO jm = 1, m

            DO jd = 1, d

              DO ja = 1, a

                dist = &
                  (coor(jf, jm, ja, jd) - &
                  coor(jf, jm, 1, jd))/bohr2ang

                dist_cell = &
                  (coor(jf, jm, ja, jd))/bohr2ang

                dipl_mol(jf, jm, jd) = &
                  dipl_mol(jf, jm, jd) + &
                  dipl(jf, jm, ja, jd) + &
                  dist*chrg(jf, jm, ja, my_d_chrg)

                dipl_cell(jf, my_m_dipl_cell, jd) = &
                  dipl_cell(jf, my_m_dipl_cell, jd) + &
                  (dipl(jf, jm, ja, jd) + &
                  dist_cell*chrg(jf, jm, ja, my_d_chrg))*au2debye

                IF (jf == 0) THEN
                  dipl_cell_temp(jd) = &
                  dipl_cell_temp(jd) + &
                  (dipl(jf, jm, ja, jd) + &
                  dist_cell*chrg(jf, jm, ja, my_d_chrg))*au2debye
                END IF

              END DO

              dipl_mol(jf, jm, jd) = dipl_mol(jf, jm, jd)*au2debye
              dipl_mol(jf, jm, my_d_dipl) = &
                dipl_mol(jf, jm, my_d_dipl) + &
                dipl_mol(jf, jm, jd)**2

              dipl_sys(jf, my_m_dipl_cell, jd) = &
                dipl_sys(jf, my_m_dipl_cell, jd) + &
                dipl_mol(jf, jm, jd)

            END DO

            dipl_mol(jf, jm, my_d_dipl) = &
              SQRT(dipl_mol(jf, jm, my_d_dipl))

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (dipl_mol(jf, jm, jd),jd=1,my_d_dipl)

            IF (jf == 0) THEN

              WRITE(9,*) &
                jf, tim(jf), jm, &
                (dipl_mol(jf, jm, jd),jd=1,my_d_dipl)

                X = X + dipl_mol(jf, jm, 1)
                Y = Y + dipl_mol(jf, jm, 2)
                Z = Z + dipl_mol(jf, jm, 3)

            END IF

          END DO

        END DO

        f_counter = 0
        ms_dipl_cell = zero
        DO jf = f0+1, f
          f_counter = f_counter + 1
          DO jd = 1, d
            ms_dipl_cell = &
              ms_dipl_cell + dipl_cell(jf, my_m_dipl_cell, jd)**2
          END DO
        END DO
        ms_dipl_cell = ms_dipl_cell/REAL(f_counter,dp)
        rms_dipl_cell = SQRT(ms_dipl_cell)

! put all quantities in atomic units for calculating dielectric_const

        pi = ACOS(-1.0_dp)
        ms_dipl_cell = ms_dipl_cell/(au2debye**2)
        temperature = temperature/hart2kelvin
        volume = 1.0_dp
        DO jd = 1, d
          volume = volume*cell(jd)/bohr2ang
        END DO
        dielectric_const = &
          1.0_dp + (4.0_dp*pi/3.0_dp) * &
          ms_dipl_cell / (volume*temperature)

        WRITE(6,150) '# temperature [Hartree] = ', temperature
        WRITE(6,150) '# volume [bohr^3] = ', volume
        WRITE(6,150) '# ms_dipl_cell [au^2] = ', ms_dipl_cell
        WRITE(6,150) '# rms_dipl_cell [Debye^2] = ', rms_dipl_cell
        WRITE(6,150) '# dielectric_const [unitless] = ', &
          dielectric_const

        WRITE(9,*) 'X, Y, Z = ', X, Y, Z

        WRITE(10,*) 'dipl_cell_temp = ', (dipl_cell_temp(jd),jd=1,d)

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( dipl_cell_temp )

100     FORMAT(1(1x,i10),1(1x,e20.10),1(1x,i4),4(1x,e20.10))
150     FORMAT(1(1x,a),1(1x,e20.10))

        END SUBROUTINE dipole

        END MODULE dipole_mod
