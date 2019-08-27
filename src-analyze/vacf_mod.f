        MODULE vacf_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE vacf(my_file)

        CHARACTER(len=80) my_file

        INTEGER :: f0, jf, jm, jd
        REAL(kind=dp) :: conversion_factor

        WRITE(6,'(2a)') '# SUBROUTINE = vacf; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        conversion_factor = bohr2ang/(autime2fs*fs2ps)
        conversion_factor = conversion_factor**2

        f0 = f - f_avg

        OPEN(3,file=TRIM(my_file),status='unknown')

        vacf_mol = zero

        WRITE(3,'(a)') '#  1f  2tim  3vacf_mol  3vacf_mol_normed'

        DO jf = f0, f

          DO jm = 1, m

            DO jd = 1, d
              vacf_mol(jf) = &
                vacf_mol(jf) + &
                velo_mol(f0,jm,jd) * velo_mol(jf,jm,jd) * conversion_factor
            END DO

          END DO

          vacf_mol(jf) = vacf_mol(jf)/REAL(m,dp)

          WRITE(3,100) jf, tim(jf), vacf_mol(jf), vacf_mol(jf)/vacf_mol(f0)

        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,e20.10),2(1x,e20.10))

        END SUBROUTINE vacf

        END MODULE vacf_mod
