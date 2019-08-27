        MODULE get_hbond_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE get_hbond(name2, pcounts, &
          num_hbonds_22, num_hbonds_aadd, my_file)

        CHARACTER(len=10) name2
        REAL(KIND=dp), DIMENSION(0:aa,0:dd,0:b) :: pcounts
        REAL(KIND=dp), DIMENSION(0:b) :: &
          num_hbonds_22, num_hbonds_aadd
        CHARACTER(len=80) my_file

        INTEGER :: f0, jf, jm1, jm2, jm, ja
        INTEGER :: jaa, jdd
        INTEGER :: jb
        REAL(KIND=dp) :: rcut, qcut, qcut_sup
        REAL(KIND=dp) :: rOO, r1a, r1b, r1, r2, q, q_sup
        REAL(KIND=dp) :: z_temp, m_check
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          rOOvec, r1avec, r1bvec, r1vec, r2vec, &
          rOOhat, r1ahat, r1bhat, r1hat, r2hat
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: acceptor, donor
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: counts
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          sum_counts, sum_pcounts
        REAL(KIND=dp) :: increment, contrib, pi

        WRITE(6,'(2a)') '# SUBROUTINE = get_hbond; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        increment = 1.0_dp
        contrib = 1.0_dp
        pi = ACOS(-1.0_dp)

        SELECT CASE (TRIM(name2))
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

        ALLOCATE( sum_counts(0:b) )
        sum_counts = zero
        ALLOCATE( sum_pcounts(0:b) )
        sum_pcounts = zero

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

        ALLOCATE( counts(0:aa,0:dd,0:b) )
        counts = zero
        pcounts = zero

        f0 = f - f_avg
        DO jf = f0+1, f

          acceptor = zero
          donor = zero

          DO jm1 = 1, m ! in the role of acceptor: O_1, H_1, H_2

            DO jm2 = 1, m ! in the role of donor: O_2, H_3, H_4

              IF ( jm1 == jm2 ) CYCLE

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

              IF ( r1a < r1b ) THEN

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

              IF ( TRIM(name2) == 'wernet' ) THEN
                r1vec = rOOvec
                r1 = rOO
                r1hat = rOOhat
              END IF

              q = ACOS( DOT_PRODUCT(r1hat,r2hat) ) * 180.0_dp / pi
              q_sup = 180.0_dp - q

              IF ( TRIM(name2) == 'wernet' ) THEN
                rcut = zero
                IF ( ABS(q_sup) <= qcut_sup ) &
                  rcut = -0.00044_dp*q_sup**2 + 3.3_dp
              END IF

              IF ( (r1 > rcut) .OR. (ABS(q) < qcut) ) CYCLE

              acceptor(jm1) = acceptor(jm1) + increment
              donor(jm2) = donor(jm2) + increment

            END DO

          END DO

          DO jm = 1, m

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

              jaa = INT( acceptor(jm) ) ! NINT
              jdd = INT( donor(jm) ) ! NINT

              IF ( ( (jaa >= 0) .AND. (jaa <= aa) ) .AND. &
                   ( (jdd >= 0) .AND. (jdd <= dd) ) ) THEN
                counts(jaa,jdd,jb) = counts(jaa,jdd,jb) + contrib

              ELSE
                WRITE(9,*) 'name2 = ', name2 
                WRITE(9,*) 'jf, jm, jb, b = ', jf, jm, jb, b
                WRITE(9,*) 'jaa, jdd, aa, dd = ', jaa, jdd, aa, dd
                WRITE(9,*) 'counts(jaa,jdd,jb) = ', counts(jaa,jdd,jb)
                WRITE(9,*)
                CALL FLUSH(9)

              END IF

            ELSE

              WRITE(8,*) 'name2 = ', name2 
              WRITE(8,*) 'jf, jm, jb, b = ', jf, jm, jb, b
              WRITE(8,*) 'coor(jf,jm,1,d), z_temp = ', &
                coor(jf,jm,1,d), z_temp
              WRITE(8,*)
              CALL FLUSH(8)

            END IF

          END DO

        END DO

        counts = counts / REAL(f_avg,dp) ! since were summing over f_avg frames

        DO jb = 0, b
          sum_counts(jb) = zero
          DO jaa = 0, aa
            DO jdd = 0, dd
              sum_counts(jb) = sum_counts(jb) + counts(jaa,jdd,jb)
            END DO
          END DO
          DO jaa = 0, aa
            DO jdd = 0, dd
              pcounts(jaa,jdd,jb) = zero
              IF ( ABS(sum_counts(jb)) > eps ) &
              pcounts(jaa,jdd,jb) = &
                100.0_dp * counts(jaa,jdd,jb) / sum_counts(jb)
            END DO
          END DO
          sum_pcounts(jb) = zero
          num_hbonds_22(jb) = zero
          num_hbonds_aadd(jb) = zero
          DO jaa = 0, aa
            DO jdd = 0, dd
              sum_pcounts(jb) = sum_pcounts(jb) + pcounts(jaa,jdd,jb)
              IF ( (jaa <= 2) .AND. (jdd <= 2) ) &
                num_hbonds_22(jb) = num_hbonds_22(jb) + &
                  REAL(jaa+jdd, dp)*pcounts(jaa,jdd,jb)/100.0_dp
              num_hbonds_aadd(jb) = num_hbonds_aadd(jb) + &
                REAL(jaa+jdd, dp)*pcounts(jaa,jdd,jb)/100.0_dp
           END DO
          END DO
        END DO

        OPEN(3,file=TRIM(my_file),status='unknown')

        DO jb = 0, b

          WRITE(3,200) &
            'jb, z(jb), sum_counts(jb), sum_pcounts(jb) = ', &
             jb, z(jb), sum_counts(jb), sum_pcounts(jb)
          WRITE(3,*)

          WRITE(3,'(a)') '# (counts(jaa,jdd,jb),jdd=0,dd)'
          DO jaa = 0, aa
            WRITE(3,340) (counts(jaa,jdd,jb),jdd=0,dd)
          END DO
          WRITE(3,*)

          WRITE(3,'(a)') '# (pcounts(jaa,jdd,jb),jdd=0,dd)'
          DO jaa = 0, aa
            WRITE(3,340) (pcounts(jaa,jdd,jb),jdd=0,dd)
          END DO
          WRITE(3,*)

          WRITE(3,120) &
          '# num_hbonds_22(jb), num_hbonds_aadd(jb) = ', &
             num_hbonds_22(jb), num_hbonds_aadd(jb)
          WRITE(3,*)
          WRITE(3,*)

! Additional checking: check for presence of over-coordinated molecules
          DO jaa = 0, aa
            DO jdd = 0, dd
              IF ( ( (jaa >=3) .OR. (jdd >= 3) ) .AND. &
                 ( ABS(counts(jaa,jdd,jb)) > eps ) ) THEN
                WRITE(10,*) 'name2 = ', name2
                WRITE(10,900) &
                'jaa, jdd, jb, z, counts, pcounts = ', &
                 jaa, jdd, jb, &
                 z(jb), counts(jaa,jdd,jb), pcounts(jaa,jdd,jb)
                CALL FLUSH(10)
              END IF
            END DO
          END DO

        END DO

        WRITE(3,220) &
          '# SUM(counts), SUM(pcounts), SUM(sum_counts) = ', &
             SUM(counts), SUM(pcounts), SUM(sum_counts)

        CALL FLUSH(3)
        CLOSE(3)

        m_check = SUM(counts)
        IF ( ABS(m_check - REAL(m,dp) ) > eps ) THEN
          WRITE(8,*) 'name2 = ', name2
          WRITE(8,*) 'm, m_check = ', m, m_check
          FLUSH(8)
        END IF

        jaa = 0
        jdd = 0
        my_file='hbond_' // TRIM(name2) //  '.00'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 1
        jdd = 0
        my_file='hbond_' // TRIM(name2) //  '.A'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 0
        jdd = 1
        my_file='hbond_' // TRIM(name2) //  '.D'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 1
        jdd = 1
        my_file='hbond_' // TRIM(name2) //  '.DA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 2
        jdd = 0
        my_file='hbond_' // TRIM(name2) //  '.AA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 0
        jdd = 2
        my_file='hbond_' // TRIM(name2) //  '.DD'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 2
        jdd = 1
        my_file='hbond_' // TRIM(name2) //  '.DAA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 1
        jdd = 2
        my_file='hbond_' // TRIM(name2) //  '.ADD'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        jaa = 2
        jdd = 2
        my_file='hbond_' // TRIM(name2) //  '.DDAA'
        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,800) z(jb), pcounts(jaa,jdd,jb)
        END DO
        CLOSE(3)

        WRITE(6,'(2a)') '# name2 = ', TRIM(name2)
        WRITE(6,210) '# rcut, qcut = ', rcut, qcut
        WRITE(6,300) '# f_avg, f, f0+1 = ', &
          f_avg, f, f0+1
        WRITE(6,100) '# m = ', m
        WRITE(6,220) &
          '# SUM(counts), SUM(pcounts) = ', &
             SUM(counts), SUM(pcounts)
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( sum_counts )
        DEALLOCATE( sum_pcounts )

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

100     FORMAT(1x,a,1(1x,i10))
120     FORMAT(1x,a,2(1x,f10.4))
200     FORMAT(1x,a,1(1x,i10),3(1x,f10.4))
210     FORMAT(1x,a,2(1x,f10.4))
220     FORMAT(1x,a,3(1x,f10.4))
230     FORMAT(1x,a,2(1x,i10))
300     FORMAT(1x,a,3(1x,i10))
330     FORMAT((1x,i10),2(1x,f10.4))
340     FORMAT(5(1x,f10.4))
700     FORMAT(1x,a,a23,a35)
710     FORMAT(3(1x,i5),10(1x,f10.4))
800     FORMAT(2(1x,f10.4))
900     FORMAT(1x,a,3(1x,i5),3(1x,f10.4))

        END SUBROUTINE get_hbond

        END MODULE get_hbond_mod
