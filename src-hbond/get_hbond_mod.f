        MODULE get_hbond_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE get_hbond(name1, my_file)

        CHARACTER(len=10) name1
        CHARACTER(len=80) my_file

        INTEGER :: f0, jf, jm1, jm2, jm, ja, &
          jaa, jdd, &
          first
        REAL(KIND=dp) :: increment, &
          rcut, qcut, qcut_sup, & ! r, theta cutoffs
          Z
        REAL(KIND=dp) :: &
          rOO, r1a, r1b, r1, r2, q, q_sup
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          rOOvec, r1avec, r1bvec, r1vec, r2vec, &
          rOOhat, r1ahat, r1bhat, r1hat, r2hat
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: acceptor, donor
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: counts, pcounts ! percentage counts
        REAL(KIND=dp) :: pi

        WRITE(6,'(2a)') '# SUBROUTINE = get_hbond; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        pi = ACOS(-1.0_dp)

        SELECT CASE (TRIM(name1))
        CASE ('160')
          qcut = 160.0_dp
          rcut = 2.27_dp
        CASE ('150')
          qcut = 150.0_dp
          rcut = 2.27_dp
        CASE ('140')
          qcut = 140.0_dp
          rcut = 2.27_dp
        CASE ('wernet')
          qcut = 140.0_dp
          rcut = zero ! for now
        END SELECT
        qcut_sup =  180.0_dp - qcut
        increment = 1.0_dp

        ALLOCATE( rOOvec(d) )
        ALLOCATE( r1avec(d) )
        ALLOCATE( r1bvec(d) )
        ALLOCATE( r1vec(d) )
        ALLOCATE( r2vec(d) )
        ALLOCATE( rOOhat(d) )
        ALLOCATE( r1ahat(d) )
        ALLOCATE( r1bhat(d) )
        ALLOCATE( r1hat(d) )
        ALLOCATE( r2hat(d) )
        rOOvec = zero
        r1avec = zero
        r1bvec = zero
        r1vec = zero
        r2vec = zero
        rOOhat = zero
        r1ahat = zero
        r1bhat = zero
        r1hat = zero
        r2hat = zero

        ALLOCATE( acceptor(m) )
        ALLOCATE( donor(m) )
        acceptor = zero
        donor = zero

        ALLOCATE( counts(0:aa,0:dd) )
        ALLOCATE( pcounts(0:aa,0:dd) )
        counts = zero
        pcounts = zero

        first = 1
        f0 = f - f_avg
        DO jf = f0+1, f

          DO jm1 = 1, m ! 2, 2 ! 1, m ! in the role of acceptor: O_1, H_1, H_2
            DO jm2 = 1, m ! 19, 19 ! 1, m ! in the role of donor: O_2, H_3, H_4

              IF (jm1 == jm2) CYCLE

              rOOvec = coor(jf,jm2,1,:) - coor(jf,jm1,1,:) ! O_2-O_1
              rOOvec(:) = rOOvec(:) - cell(:)*ANINT(rOOvec(:)/cell(:)) ! do not forget to appy PBCs
              rOO = SQRT( DOT_PRODUCT(rOOvec,rOOvec) )
              rOOhat = rOOvec/rOO

              r1avec = coor(jf,jm2,2,:) - coor(jf,jm1,1,:) ! H_3-O_1
              r1avec(:) = r1avec(:) - cell(:)*ANINT(r1avec(:)/cell(:)) ! do not forget to appy PBCs
              r1a = SQRT( DOT_PRODUCT(r1avec,r1avec) )
              r1ahat = r1avec/r1a

              r1bvec = coor(jf,jm2,3,:) - coor(jf,jm1,1,:) ! H_4-O_1
              r1bvec(:) = r1bvec(:) - cell(:)*ANINT(r1bvec(:)/cell(:)) ! do not forget to appy PBCs
              r1b = SQRT( DOT_PRODUCT(r1bvec,r1bvec) )
              r1bhat = r1bvec/r1b

              IF ( r1a < r1b) THEN

                r1vec = r1avec
                r1 = r1a
                r1hat = r1ahat

                r2vec = coor(jf,jm2,2,:) - coor(jf,jm2,1,:) ! O_2-H_3
                r2vec(:) = r2vec(:) - cell(:)*ANINT(r2vec(:)/cell(:)) ! do not forget to appy PBCs
                r2 = SQRT( DOT_PRODUCT(r2vec,r2vec) )
                r2hat = r2vec/r2

              ELSE

                r1vec = r1bvec
                r1 = r1b
                r1hat = r1bhat

                r2vec = coor(jf,jm2,3,:) - coor(jf,jm2,1,:) ! O_2-H_4
                r2vec(:) = r2vec(:) - cell(:)*ANINT(r2vec(:)/cell(:)) ! do not forget to appy PBCs
                r2 = SQRT( DOT_PRODUCT(r2vec,r2vec) )
                r2hat = r2vec/r2

              END IF

              IF ( TRIM(name1) == 'wernet' ) THEN
                r1vec = rOOvec
                r1 = rOO
                r1hat = rOOhat
              END IF

              q = ACOS( DOT_PRODUCT(r1hat,r2hat) ) * 180.0_dp / pi
              q_sup = 180.0_dp - q

              IF (TRIM(name1) == 'wernet') THEN
                rcut = zero
                IF (ABS(q_sup) <= qcut_sup) &
                  rcut = -0.00044_dp*q_sup**2 + 3.3_dp
              END IF

!             IF (first == 1) THEN
!               WRITE(6,700) &
!               'jf, jm1, jm2, ', &
!               'rcut, qcut, qcut_sup, ', &
!               'rOO, r1a, r1b, r1, r2, q, q_sup  = '
!               first = 0
!             END IF
!             WRITE(6,710) &
!                jf, jm1, jm2, &
!                rcut, qcut, qcut_sup, &
!                rOO, r1a, r1b, r1, r2, q, q_sup

!             IF ( (jf == 1) .AND. (jm1 == 2 ) .AND. (jm2 == 19) ) THEN
!               OPEN(9,file='coor_f1_m2_m19.xyz',status='unknown')
!               WRITE(9,*) 2*a
!               WRITE(9,*) 'coor_f1_m2_m19'
!               DO ja = 1, a
!                 WRITE(9,'(a,1x,3(f20.10))') &
!                   symbol(jm1,ja), coor(jf,jm1,ja,:)
!               END DO
!               DO ja = 1, a
!                 WRITE(9,'(a,1x,3(f20.10))') &
!                   symbol(jm2,ja), coor(jf,jm2,ja,:)
!               END DO
!               CALL FLUSH(9)
!               CLOSE(9)
!             END IF

              IF ( (r1 > rcut) .OR. (ABS(q) < qcut) ) CYCLE
              
              acceptor(jm1) = acceptor(jm1) + increment
              donor(jm2) = donor(jm2) + increment

            END DO
          END DO
        END DO

        acceptor = acceptor / REAL(f_avg,dp)
        donor = donor / REAL(f_avg,dp)

        DO jm = 1, m
          jaa = NINT( acceptor(jm) )
          jdd = NINT( donor(jm) )
          IF ( ( (jaa >= 0) .AND. (jaa <= aa) ) .AND. &
               ( (jdd >= 0) .AND. (jdd <= dd) ) ) &
          counts(jaa,jdd) = counts(jaa,jdd) + 1.0_dp
        END DO
        pcounts = 100.0_dp * counts / REAL(m,dp)

        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,'(a)') '# jm, acceptor(jm), donor(jm)'
        DO jm = 1, m
          WRITE(3,330) jm, acceptor(jm), donor(jm)
        END DO
        WRITE(3,*)
        WRITE(3,*)
        WRITE(3,'(a)') '# (counts(jaa,jdd),jdd=0,dd)'
        DO jaa = 0, aa
          WRITE(3,340) (counts(jaa,jdd),jdd=0,dd)
        END DO
        WRITE(3,*)
        WRITE(3,*)
        WRITE(3,'(a)') '# (pcounts(jaa,jdd),jdd=0,dd)'
         DO jaa = 0, aa
          WRITE(3,340) (pcounts(jaa,jdd),jdd=0,dd)
        END DO
        WRITE(3,*)
        WRITE(3,*)
        WRITE(3,220) &
          '# SUM(counts), SUM(pcounts) = ', &
             SUM(counts), SUM(pcounts)
        CALL FLUSH(3)
        CLOSE(3)

        Z = zero

        jaa = 0
        jdd = 0
        my_file='hbond_' // TRIM(name1) //  '.00'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        jaa = 1
        jdd = 0
        my_file='hbond_' // TRIM(name1) //  '.A'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        jaa = 0
        jdd = 1
        my_file='hbond_' // TRIM(name1) //  '.D'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        jaa = 1
        jdd = 1
        my_file='hbond_' // TRIM(name1) //  '.DA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        jaa = 2
        jdd = 1
        my_file='hbond_' // TRIM(name1) //  '.DAA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        jaa = 1
        jdd = 2
        my_file='hbond_' // TRIM(name1) //  '.ADD'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        jaa = 2
        jdd = 2
        my_file='hbond_' // TRIM(name1) //  '.DDAA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        WRITE(3,800) Z, pcounts(jaa,jdd)    
        CLOSE(3)

        WRITE(6,'(2a)') '# name1 = ', TRIM(name1)
        WRITE(6,210) '# rcut, qcut = ', rcut, qcut
        WRITE(6,300) '# f_avg, f, f0+1 = ', &
          f_avg, f, f0+1
        WRITE(6,100) '# m = ', m
        WRITE(6,220) &
          '# SUM(acceptor), SUM(donor) = ', &
             SUM(acceptor), SUM(donor)
        WRITE(6,220) &
          '# SUM(acceptor)/REAL(m,dp), SUM(donor)/REAL(m,dp) = ', &
             SUM(acceptor)/REAL(m,dp), SUM(donor)/REAL(m,dp)
        WRITE(6,220) &
          '# SUM(counts), SUM(pcounts) = ', &
             SUM(counts), SUM(pcounts)
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( rOOvec )
        DEALLOCATE( r1avec )
        DEALLOCATE( r1bvec )
        DEALLOCATE( r1vec )
        DEALLOCATE( r2vec )
        DEALLOCATE( rOOhat )
        DEALLOCATE( r1ahat )
        DEALLOCATE( r1bhat )
        DEALLOCATE( r1hat )
        DEALLOCATE( r2hat )

        DEALLOCATE( acceptor )
        DEALLOCATE( donor )

        DEALLOCATE( counts )
        DEALLOCATE( pcounts )

100     FORMAT(1x,a40,1(1x,i10))
210     FORMAT(1x,a40,2(1x,f10.4))
220     FORMAT(1x,a80,2(1x,f10.4))
300     FORMAT(1x,a40,3(1x,i10))
330     FORMAT((1x,i10),2(1x,f10.4))
340     FORMAT(5(1x,f10.4))
700     FORMAT(1x,a15,a23,a35)
710     FORMAT(3(1x,i5),10(1x,f10.4))
800     FORMAT(2(1x,f10.4))

        END SUBROUTINE get_hbond

        END MODULE get_hbond_mod
