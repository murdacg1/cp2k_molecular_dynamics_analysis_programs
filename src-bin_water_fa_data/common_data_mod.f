        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8
        INTEGER :: t, m, n, cc, c
        INTEGER, DIMENSION(:), ALLOCATABLE :: cc_numbers
        CHARACTER(len=80), DIMENSION(:), ALLOCATABLE :: cc_names
        CHARACTER(len=80), DIMENSION(:), ALLOCATABLE :: title_lines
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: full_array

        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: eps=1.0e-12_dp

        END MODULE common_data_mod


