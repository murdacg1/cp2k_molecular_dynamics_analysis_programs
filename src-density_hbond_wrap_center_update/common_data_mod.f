        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8
        REAL(KIND=dp), PARAMETER :: zero = 0.0_dp
        REAL(KIND=dp), PARAMETER :: bohr2ang=0.529177249_dp
        REAL(KIND=dp), PARAMETER :: eps = 1.0e-12_dp

        INTEGER :: f_avg, f, m, a, &
          b, &
          aa, dd, &
          bulk_versus_surface
        INTEGER, PARAMETER :: d=3
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: cell, z
        REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: coor
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: &
          pcounts_140, pcounts_150, pcounts_160, pcounts_wernet
        CHARACTER(len=1), DIMENSION(:,:), ALLOCATABLE :: symbol
        REAL(KIND=dp) :: dz, mass, z_gds, delta
        REAL(KIND=dp), DIMENSION(d) :: &
          center_old = zero, center = zero, center_update = zero
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          num_hbonds_22_140, num_hbonds_aadd_140, &
          num_hbonds_22_150, num_hbonds_aadd_150, &
          num_hbonds_22_160, num_hbonds_aadd_160, &
          num_hbonds_22_wernet, num_hbonds_aadd_wernet

        END MODULE common_data_mod
