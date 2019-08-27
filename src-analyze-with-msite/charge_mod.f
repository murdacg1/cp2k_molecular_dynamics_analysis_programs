        MODULE charge_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE charge(my_d_chrg, &
          my_file)

        INTEGER :: my_d_chrg
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, ja, jd

        WRITE(6,'(2a)') '# SUBROUTINE = charge; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        chrg_mol = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(2a)') &
          '#  1f  2tim  3m  ', &
          '4chrg_mol'

        DO jf = f0+1, f
          DO jm = 1, m

            DO jd = 1, my_d_chrg
              DO ja = 1, a
                chrg_mol(jf,jm,jd) = &
                  chrg_mol(jf,jm,jd) +  &
                  chrg(jf,jm,ja,jd)
              END DO
            END DO

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (chrg_mol(jf,jm,jd),jd=1,my_d_chrg)

          END DO
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,e20.10),1(1x,i4),7(1x,e20.10))

        END SUBROUTINE charge

        END MODULE charge_mod
