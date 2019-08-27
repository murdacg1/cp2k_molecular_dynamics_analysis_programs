        PROGRAM combine_histograms

! Purpose: Line up / combine in overlapping histograms from n simulations and obtain pmf(z) - pmf(z_0) = -k*T*ln(p(z)/p(z_0)
! First read in n histograms that came from processing n CP2K z_s umbrella sampling simulation of water slab with a single solute molecule

! Program is run as:
!       ../combine_histograms.exe < combine_histograms.in >& combine_histograms.out &
!
! Program reads first the file combine_histograms.in:

! 0                           n1 (first dir)
! 21                          n2 (last dir)
! 100.0                       factor
! 0.1                         width (bin width for individual histograms) [Angstrom]
! slab_spce_spackman_NO2_z    dir1
! _300K_with_shift_umbrella   dir2
! 30                          dir3 (run time analyzed = equilibration + production, in ps)
! 300.0                       temperature [Kelvin]

        IMPLICIT NONE

        INTEGER, PARAMETER :: jmax=1000
        INTEGER, PARAMETER :: dp=8
        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: hart2kcal=627.51_dp, hart2kelvin=315773.3_dp, &
          kelvin2kcal=hart2kcal/hart2kelvin
        REAL(KIND=dp), PARAMETER :: eps=1.0e-12_dp

        INTEGER :: n1, n2, jn, z1, z2, zstride, jz
        REAL(KIND=dp) :: factor, width, temperature
        CHARACTER(len=200) dir1, dir2, dir3, my_file, string80, jn_char

        INTEGER :: jc_temp, counts_temp, jj
        REAL(KIND=dp) :: center_temp, ncounts_temp

        INTEGER, DIMENSION(:), ALLOCATABLE :: counts ! needed distribution
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: ncounts, pmf ! normalized distribution, pmf
        REAL(KIND=dp) :: area

        WRITE(6,'(a)') '# PROGRAM = combine_histograms'
        CALL FLUSH(6)

        READ(5,*) n1
        READ(5,*) n2
        READ(5,*) nrev
        ALLOCATE( vec_reverse(nr) )
        DO jrev = 1, nr
          READ(5,*) vecrev(jrev)
        END DO
        READ(5,*) factor
!       z1 = &
!         NINT(factor*( REAL(n1,dp)-0.5_dp*(REAL(n2,dp)-REAL(n1,dp)) ))
!       z2 = &
!         NINT(factor*( REAL(n2,dp)+0.5_dp*(REAL(n2,dp)-REAL(n1,dp)) ))
        z1 = NINT(-10.5_dp*factor)
        z2 = NINT(31.5_dp*factor)
        READ(5,*) width
        zstride = NINT(1.0_dp/width)
        READ(5,*) dir1
        READ(5,*) dir2
        READ(5,*) dir3
        dir1 = TRIM(dir1)
        dir2 = TRIM(dir2)
        dir3 = TRIM(dir3)
        READ(5,*) temperature
        temperature = kelvin2kcal*temperature

        WRITE(6,120) '# n1, n2 = ', n1, n2
        WRITE(6,210) '# factor = ', factor
        WRITE(6,130) '# z1, z2, zstride = ', z1, z2, zstride
        WRITE(6,'(a,1x,a,1x,a,1x,a)') &
          '# dir1, dir2, dir3 = ', dir1, dir2, dir3
        WRITE(6,210) '# width = ', width
        WRITE(6,210) '# temperature [kcal/mol] = ', temperature
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( counts(z1:z2) )
        ALLOCATE( ncounts(z1:z2) )
        ALLOCATE( pmf(z1:z2) )
        counts = 0
        ncounts = zero
        pmf = zero

        DO jn = n1, n2
          
          factor_rev = 1.0d0
          
!         jn_char = TRIM(CHAR(jn))
 
          SELECT CASE (jn)
          CASE (-6)
            jn_char = '-6'
          CASE (-5) 
            jn_char = '-5'
          CASE (-4) 
            jn_char = '-4'
          CASE (-3) 
            jn_char = '-3'
          CASE (-2) 
            jn_char = '-2'
          CASE (-1) 
            jn_char = '-1'
          CASE (0)
            jn_char = '0'
          CASE (1)
            jn_char = '1'
          CASE (2)
            jn_char = '2'
          CASE (3)
            jn_char = '3'
          CASE (4)
            jn_char = '4'
          CASE (5)
            jn_char = '5'
          CASE (6)
            jn_char = '6'
          CASE (7)
            jn_char = '7'
          CASE (8)
            jn_char = '8'
          CASE (9)
            jn_char = '9'
          CASE (10) 
            jn_char = '10'
          CASE (11) 
            jn_char = '11'
          CASE (12) 
            jn_char = '12'
          CASE (13) 
            jn_char = '13'
          CASE (14) 
            jn_char = '14'
          CASE (15) 
            jn_char = '15'
          CASE (16) 
            jn_char = '16'
          CASE (17) 
            jn_char = '17'
          CASE (18) 
            jn_char = '18'
          CASE (19) 
            jn_char = '19'
          CASE (20)
            jn_char = '20'
          CASE (21)
            jn_char = '21'
          CASE (22)
            jn_char = '22'
          CASE (23)
            jn_char = '23'
          CASE (24)
            jn_char = '24'
          CASE (25)
            jn_char = '25'
          CASE (26)
            jn_char = '26'
          END SELECT

          jn_char = TRIM(jn_char)

          WRITE(6,*) '# jn, jn_char = ', jn, jn_char
          CALL FLUSH(6)

          my_file = &
            '../' // TRIM(dir1) // TRIM(jn_char) // TRIM(dir2) // &
            '/BINZ_' // TRIM(dir3) // 'ps/z.bin'
          my_file = TRIM(my_file)
          WRITE(6,*) '# my_file = ', my_file
          CALL FLUSH(6)

          OPEN(3,file=TRIM(my_file),status='unknown')

          READ(3,*) string80

          DO jj = 1, jmax
            READ(3,*,END=999) &
              jc_temp, center_temp, counts_temp, ncounts_temp
            jz = NINT(factor*center_temp)
            IF ( ABS( (REAL(jz,dp)/factor) - center_temp ) > eps ) THEN
              WRITE(6,*) 'Problem: center_temp, REAL(jz,dp)/factor = ',&
                                   center_temp, REAL(jz,dp)/factor
              CALL FLUSH(6)
              STOP
            ELSE
              counts(jz) = counts(jz) + counts_temp
            END IF
          END DO

999       CLOSE(3)

        END DO

        area = REAL(sum(counts),dp)*width
        DO jz = z1, z2
          ncounts(jz) = REAL(counts(jz),dp)/area
          IF (counts(jz) > zero) &
            pmf(jz) = -temperature*LOG(ncounts(jz))
        END DO

        my_file = 'pmf.dat'
        OPEN(4,file=TRIM(my_file),status='unknown')
        WRITE(4,'(a)') '# z  counts  normalized_counts  pmf [kcal/mol]'
        CALL FLUSH(4)
        DO jz = z1, z2, zstride
          IF ( counts(jz) > 0 ) &
          WRITE(4,400) &
          REAL(jz,dp)/factor, counts(jz), ncounts(jz), pmf(jz)
          CALL FLUSH(4)
        END DO
        CLOSE(4)

        WRITE(6,210) '# bin width = ', width
        WRITE(6,210) '# area = ', area
        WRITE(6,210) '# area/width = ', area/width
        WRITE(6,110) '# sum(counts) = ', sum(counts)
        WRITE(6,210) '# sum(ncounts)*width = ', sum(ncounts)*width
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( counts )
        DEALLOCATE( ncounts )
        DEALLOCATE( pmf )

110     FORMAT(1x,a40,1(1x,i20))
120     FORMAT(1x,a40,2(1x,i20))
130     FORMAT(1x,a40,3(1x,i20))
210     FORMAT(1x,a40,1(1x,e20.10))
220     FORMAT(1x,a40,2(1x,e20.10))
330     FORMAT(3(1x,a40))
400     FORMAT(1x,e20.10,1x,i20,2(1x,e20.10))

        END PROGRAM combine_histograms
