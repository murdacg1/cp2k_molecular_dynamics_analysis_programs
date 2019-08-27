        MODULE msd_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE msd(my_file)

        CHARACTER(len=80) my_file

        INTEGER :: f0, jf, jm, jd
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: dR

        WRITE(6,'(2a)') '# SUBROUTINE = msd; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        msd_mol = zero

        ALLOCATE( dR(d) )

        dR = zero

        f0 = f - f_avg

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(a)') '#  1f  2tim  3msd_mol'

        DO jf = f0+1, f

          DO jm = 1, m

            DO jd = 1, d
              dR(jd) = coor_mol(jf,jm,jd) - coor_mol(f0,jm,jd)
!             dR(jd) = dR(jd) - cell(jd)*ANINT(dR(jd)/cell(jd))
              msd_mol(jf) = msd_mol(jf) +  dR(jd)**2
            END DO

          END DO

          msd_mol(jf) = msd_mol(jf)/REAL(m,dp)

          WRITE(3,100) jf, tim(jf), msd_mol(jf)

        END DO

        CLOSE(3)

        DEALLOCATE( dR )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,e20.10),1(1x,e20.10))

        END SUBROUTINE msd

        END MODULE msd_mod
