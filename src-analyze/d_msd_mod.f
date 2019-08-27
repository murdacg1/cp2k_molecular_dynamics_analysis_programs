        MODULE d_msd_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE d_msd

        INTEGER :: f0
        REAL(KIND=dp) Delta_t, Delta_msd_mol, d_msd_mol

        WRITE(6,'(a)') '# SUBROUTINE = d_msd'
        CALL FLUSH(6)

        f0 = f - f_avg

        Delta_t = tim(f) - tim(f0)
        Delta_msd_mol = msd_mol(f) - msd_mol(f0)

        d_msd_mol = (Delta_msd_mol/Delta_t)/6.0_dp

        WRITE(6,100) '# d_msd_mol = ', d_msd_mol

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,a),1(1x,e20.10))

        END SUBROUTINE d_msd

        END MODULE d_msd_mod
