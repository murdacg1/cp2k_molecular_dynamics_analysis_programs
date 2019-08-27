        MODULE d_vacf_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE d_vacf

        INTEGER :: jf
        REAL(KIND=dp) :: d_vacf_mol, dt

        WRITE(6,'(a)') '# SUBROUTINE = d_vacf'
        CALL FLUSH(6)

        dt = tim(1) - tim(0)
        d_vacf_mol = zero

        DO jf = f0, f
          d_vacf_mol = d_vacf_mol + vacf_mol(jf)
        END DO
        d_vacf_mol = d_vacf_mol*dt/3.0_dp

        WRITE(6,100) '# d_vacf_mol = ', d_vacf_mol

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,a),1(1x,e20.10))

        END SUBROUTINE d_vacf

        END MODULE d_vacf_mod
