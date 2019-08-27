        MODULE get_msite_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE get_msite

        INTEGER :: jf, jm
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          r1vec, r2vec, rbvec, rbhat
        REAL(KIND=dp) :: rb

        WRITE(6,'(2a)') '# SUBROUTINE = get_msite'
        CALL FLUSH(6)

        ALLOCATE( r1vec(d) )
        ALLOCATE( r2vec(d) )
        ALLOCATE( rbvec(d) )
        ALLOCATE( rbhat(d) )
        r1vec = zero
        r2vec = zero
        rbvec = zero
        rbhat = zero

        DO jf = 0, f

          DO jm = 1, m

            r1vec = coor(jf,jm,2,:) - coor(jf,jm,1,:) ! H_1-O_1
            r2vec = coor(jf,jm,3,:) - coor(jf,jm,1,:) ! H_2-O_1
            rbvec = r1vec + r2vec
            rb = SQRT( DOT_PRODUCT(rbvec,rbvec) )
            rbhat = rbvec
            coor(jf,jm,4,:) = coor(jf,jm,1,:) + dist*rbhat(:)

          END DO

        END DO

        DEALLOCATE( r1vec )
        DEALLOCATE( r2vec )
        DEALLOCATE( rbvec )
        DEALLOCATE( rbhat )

        END SUBROUTINE get_msite

        END MODULE get_msite_mod
