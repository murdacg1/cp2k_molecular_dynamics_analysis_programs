        MODULE dist_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE dist

        INTEGER :: jn, jf, jm, ja1, ja2
        CHARACTER(len=80) :: filename12
        REAL(KIND=dp) :: dx, dy, dz, r

        WRITE(6,'(a)') '# SUBROUTINE = dist'
        CALL FLUSH(6)

        jm = 1
        WRITE(6,'(a)') 'jn, jf, dx1, dx, ABS(dx-dx1), cell(1):'
        DO jn = 1, n_dist
          ja1 = atom1_array(jn) + 1
          ja2 = atom2_array(jn) + 1
          filename12 = filename12_array(jn)
          OPEN(3,file=TRIM(filename12),status='unknown')
          DO jf = 0, f
            dx = coor(jf,jm,ja2,1) - coor(jf,jm,ja1,1)
            dy = coor(jf,jm,ja2,2) - coor(jf,jm,ja1,2)
            dz = coor(jf,jm,ja2,3) - coor(jf,jm,ja1,3)
            dx = dx - cell(1)*ANINT(dx/cell(1))
            dy = dy - cell(2)*ANINT(dy/cell(2))
            dz = dz - cell(3)*ANINT(dz/cell(3))
            r = SQRT(dx**2 + dy**2 + dz**2)
            WRITE(3,100) jf, r
          END DO
          CLOSE(3)
        END DO

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,f12.6))

        END SUBROUTINE dist

        END MODULE dist_mod
