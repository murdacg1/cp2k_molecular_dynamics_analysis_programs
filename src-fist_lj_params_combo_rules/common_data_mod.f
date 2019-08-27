        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8
        INTEGER :: n, sigma_rule, epsilon_rule
        CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: symbol
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: charge
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: sigma, epsilon

        REAL(KIND=dp), PARAMETER :: zero=0.0_dp

        END MODULE common_data_mod
