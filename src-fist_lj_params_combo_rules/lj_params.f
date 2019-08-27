        PROGRAM lj_params

        USE common_data_mod

! Purpose: Read in Lennard-Jones parameters for a classical force field and print out the parameters
! (including those obtained with combining rules) in the form expected by the CP2K/Fist package.
! LJ parameters sigma and epsilon can be read in with various units. Then, the read in conversion factors
! allow the data to be converted to sigma = Angstroms and epsilon = kcal/mol.
!
! Program is run as:
!       lj_params.exe < lj_params.in >& lj_params.out &
!
! Note about Lennard-Jones combining rules (see: https://en.wikipedia.org/wiki/Combining_rules):
!   Lorentz-Berthelot combining rules (Jedlovszky group uses this):
!     sigma(i,j)   = arithmetic_avg (Lorentz   rule) =  0.5_dp*( sigma(i,i) + sigma(j,j) )      (sigma_rule   = 1)
!     epsilon(i,j) = geometric_avg  (Berthelot rule) =     sqrt( epsilon(i,i) * epsilon(j,j) )  (epsilon_rule = 2)
!   OPLS combining rules:
!     A(i,j) = geometric_avg = sqrt( A(i,i) * A(j,j) )  (A_rule = 2)
!     C(i,j) = geometric_avg = sqrt( C(i,i) * C(j,j) )  (C_rule = 2)
!
! Program reads first the file lj_params.in (Here the trans-formic acid  force field is from Jedlovszky and Turi, JPCA 103, 1999 and the water force field is SPC/E):
! 7                                                                    n               (number of kinds of sites, including zero values)
! 1                                                                    sigma_rule      (combining rule for sigma values)
! 2                                                                    epsilon_rule    (combining rule for epsilon values)
!  O  -0.8476      3.166     1.000     78.198   503.22271716452354     symbol   charge     sigma   sigma_conv     epsilon   epsilon_conv
!  H   0.4238      0.000     1.000      0.000   503.22271716452354
! C2   0.44469     3.727     1.000      0.376     4.184
! O2  -0.55296     3.180     1.000      0.392     4.184
! H2   0.43331     0.994     1.000      0.100     4.184
! O3  -0.43236     2.674     1.000      1.214     4.184
! H3   0.10732     0.800     1.000      0.020     4.184

        IMPLICIT NONE

        INTEGER :: i, j
        REAL(KIND=dp) :: sigma_temp, sigma_conv, epsilon_temp, epsilon_conv

        WRITE(6,'(a)') '# PROGRAM = lj_params'

        READ(5,*) n
        WRITE(6,100)  'n = ', n
        READ(5,*) sigma_rule
        WRITE(6,100)  'sigma_rule = ', sigma_rule
        READ(5,*) epsilon_rule
        WRITE(6,100)  'epsilon_rule = ', epsilon_rule
        ALLOCATE( symbol(n) )
        symbol = '  '
        ALLOCATE( charge(n) )
        charge = zero
        ALLOCATE( sigma(n,n) )
        sigma = zero
        ALLOCATE( epsilon(n,n) )
        epsilon = zero
        DO i = 1, n
          READ(5,*)  symbol(i), charge(i), sigma_temp, sigma_conv, epsilon_temp, epsilon_conv
          WRITE(6,200)  symbol(i), charge(i), sigma_temp, sigma_conv, epsilon_temp, epsilon_conv
          CALL FLUSH(6)
          sigma(i,i) = sigma_temp / sigma_conv
          epsilon(i,i) = epsilon_temp / epsilon_conv
          WRITE(6,300)  symbol(i), charge(i), sigma(i,i), epsilon(i,i)
        END DO

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DO i = 1, n-1
          DO j = i+1, n
            SELECT CASE (sigma_rule)
            CASE (1)
              sigma(i,j) = 0.5_dp*( sigma(i,i) + sigma(j,j) )
            CASE (2)
              sigma(i,j) = sqrt( sigma(i,i) * sigma(j,j) )
            END SELECT
            SELECT CASE (epsilon_rule)
            CASE (1)
              epsilon(i,j) = 0.5_dp*( epsilon(i,i) + epsilon(j,j) )
            CASE (2)
              epsilon(i,j) = sqrt( epsilon(i,i) * epsilon(j,j) )
            END SELECT
          END DO
        END DO

        DO i = 1, n
          WRITE(6,'(a)')           '      &CHARGE'
          WRITE(6,'(a,a2)')        '        ATOM ', symbol(i)
          WRITE(6,'(a,f10.5)')     '        CHARGE ', charge(i)
          WRITE(6,'(a)')           '      &END CHARGE'
        END DO

        WRITE(6,'(a)')             '      &NONBONDED'

        DO i = 1, n
          WRITE(6,'(a)')           '        &LENNARD-JONES'
          WRITE(6,'(a,a2,a,a2)')   '          ATOMS ', symbol(i), ' ', symbol(i)
          WRITE(6,'(a,f10.5)')     '          EPSILON [kcalmol] ', epsilon(i,i)
          WRITE(6,'(a,f10.5)')     '          SIGMA ', sigma(i,i)
          WRITE(6,'(a)')           '          RCUT 11.4'
          WRITE(6,'(a)')           '        &END LENNARD-JONES'
        END DO

        DO i = 1, n-1
          DO j = i+1, n
            WRITE(6,'(a)')         '        &LENNARD-JONES'
            WRITE(6,'(a,a2,a,a2)') '          ATOMS ', symbol(i), ' ', symbol(j)
            WRITE(6,'(a,f10.5)')   '          EPSILON [kcalmol] ', epsilon(i,j)
            WRITE(6,'(a,f10.5)')   '          SIGMA ', sigma(i,j)
            WRITE(6,'(a)')         '          RCUT 11.4'
            WRITE(6,'(a)')         '        &END LENNARD-JONES'
          END DO
        END DO

        WRITE(6,'(a)')             '      &END NONBONDED'

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( symbol )
        DEALLOCATE( charge )
        DEALLOCATE( sigma )
        DEALLOCATE( epsilon )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1(1x,a2),5(1x,f10.5))
300     FORMAT(1(1x,a2),3(1x,f10.5))

        END PROGRAM lj_params
