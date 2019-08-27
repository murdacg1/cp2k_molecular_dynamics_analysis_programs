        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8

        INTEGER :: f_avg, f, m, a, c
        INTEGER, PARAMETER :: d=3
        REAL(KIND=dp) :: temperature
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: cell
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tim, pot
        CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: symbol
        REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
          coor, velo
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: com, velo_com

        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: fs2ps=1.0e-3_dp

        END MODULE common_data_mod
