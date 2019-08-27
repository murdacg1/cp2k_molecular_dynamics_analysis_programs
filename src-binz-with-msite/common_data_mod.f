        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8

        INTEGER :: f_avg, f, f0, &
          m_w, m_s, &
          a_w, a_s, &
          which
        REAL(KIND=dp) :: z_s
        REAL(KIND=dp) :: middle, width, separation
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass_w, mass_s
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: cell
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: z

        INTEGER, PARAMETER :: d=3
        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
 
        END MODULE common_data_mod
