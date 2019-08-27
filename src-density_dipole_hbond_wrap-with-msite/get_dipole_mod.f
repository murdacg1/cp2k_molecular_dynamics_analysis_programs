        MODULE get_dipole_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE get_dipole(name3, my_d_chrg, my_file)

        CHARACTER(len=6) name3
        INTEGER :: my_d_chrg
        CHARACTER(len=80) my_file

        INTEGER :: f0, jf, jm, ja, jd, jb
        REAL(KIND=dp) :: z_temp, contrib, &
          m_check, dist

        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          dipl_mol, number_dipl_mol, sd, &
          contrib_vec
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: &
          dipl_mol_frame_molec

        WRITE(6,'(2a)') '# SUBROUTINE = get_dipole; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        ALLOCATE( dipl_mol(0:b) )
        ALLOCATE( number_dipl_mol(0:b) )
        ALLOCATE( sd(0:b) )
        ALLOCATE( dipl_mol_frame_molec(0:f,1:m,0:b) )
        dipl_mol = zero
        number_dipl_mol = zero
        sd = zero
        dipl_mol_frame_molec = zero

        ALLOCATE( contrib_vec(d) )

        f0 = f - f_avg
        DO jf = f0+1, f

          DO jm = 1, m

            contrib = zero

            DO jd = 1, d

              contrib_vec(jd) = zero 

              DO ja = 1, a

                dist = &
                  (coor(jf, jm, ja, jd) - &
                   coor(jf, jm,  1, jd)) / bohr2ang

                contrib_vec(jd) = &
                  contrib_vec(jd) + &
                  dipl(jf, jm, ja, jd) + &
                  dist*chrg(jf, jm, ja, my_d_chrg)

              END DO

              contrib_vec(jd) = contrib_vec(jd)*au2debye
              contrib = contrib + contrib_vec(jd)**2

            END DO

            contrib = SQRT(contrib)

            z_temp = coor(jf,jm,1,d)

            IF ( z_temp <= cell(d)/2.0_dp ) THEN ! use both surfaces
              z_temp = -z_temp + cell(d)/2.0_dp
            ELSE
              z_temp =  z_temp - cell(d)/2.0_dp
            END IF

            IF ( bulk_versus_surface == 0 ) THEN
              jb = INT(z_temp/dz) ! NINT
              IF ( ABS(dz - cell(d) ) < eps ) jb = 0 ! bulk only calc
            ELSE
              IF ( z_temp < (z_gds-2.0_dp*delta) ) THEN
                jb = 0 ! bulk region
              ELSE
                jb = 1 ! surface region
              END IF
            END IF

            IF  ( (jb >= 0) .AND. (jb <= b) ) THEN

              dipl_mol_frame_molec(jf,jm,jb) = contrib
              number_dipl_mol(jb) = number_dipl_mol(jb) + 1.0_dp

            ELSE

              WRITE(17,*) 'jf, jm, jb, b = ', jf, jm, jb, b
              WRITE(17,*) 'coor(jf,jm,1,d), z_temp = ', &
                coor(jf,jm,1,d), z_temp
              WRITE(17,*)
              CALL FLUSH(17)

            END IF

          END DO

        END DO

        DO jb = 0, b
          DO jf = f0+1, f
            DO jm = 1, m
              dipl_mol(jb) = &
                dipl_mol(jb) + dipl_mol_frame_molec(jf,jm,jb)
            END DO
          END DO
          IF ( number_dipl_mol(jb) > zero ) &
            dipl_mol(jb) = dipl_mol(jb) / number_dipl_mol(jb)
        END DO
        m_check = SUM(number_dipl_mol) / REAL(f_avg, dp)
        IF ( ABS(m_check - REAL(m,dp) ) > eps ) THEN
          WRITE(17,*) 'm, m_check = ', m, m_check
          FLUSH(17)
        END IF

        DO jb = 0, b
          DO jf = f0+1, f
            DO jm = 1, m
              IF ( dipl_mol_frame_molec(jf,jm,jb) > zero ) &
                sd(jb) = sd(jb) + &
                  ( dipl_mol(jb) - dipl_mol_frame_molec(jf,jm,jb) )**2
            END DO
          END DO
          IF ( number_dipl_mol(jb) > zero ) &
            sd(jb) = dsqrt( sd(jb) / number_dipl_mol(jb) )
        END DO

        WRITE(6,'(2a)') '# name3 = ', name3
        WRITE(6,100) '# m = ', m
        WRITE(6,110) '# m_check = ', m_check
        WRITE(6,100) '# a = ', a
        WRITE(6,110) '# contrib = ', contrib
        WRITE(6,300) '# f_avg, f, f0+1 = ', &
          f_avg, f, f0+1
        WRITE(6,100) '# b = ', b
        WRITE(6,110) '# dz = ', dz

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,420) z(jb), dipl_mol(jb), sd(jb), &
            number_dipl_mol(jb) / REAL(f_avg, dp)
          CALL FLUSH(3)
        END DO
        CLOSE(3)

        DEALLOCATE( dipl_mol )
        DEALLOCATE( number_dipl_mol )
        DEALLOCATE( sd )
        DEALLOCATE( dipl_mol_frame_molec )

        DEALLOCATE( contrib_vec )

100     FORMAT(1x,a40,1(1x,i10))
110     FORMAT(1x,a40,1(1x,f20.10))
200     FORMAT(1x,a40,2(1x,i10))
210     FORMAT(1x,a40,2(1x,f20.10))
300     FORMAT(1x,a40,3(1x,i10))
420     FORMAT(4(1x,e20.10))

        END SUBROUTINE get_dipole

        END MODULE get_dipole_mod
