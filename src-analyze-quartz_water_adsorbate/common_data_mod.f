        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8

        INTEGER :: f_avg, f, &
          n_si, &
          m, a, &
          m_surf, a_surf, &
          m_wat, a_wat, &
          m_sol, a_sol, &
          c, &
          n_dist
        INTEGER, PARAMETER :: d=3
        REAL(KIND=dp) :: dt, temperature
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: cell
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tim, pot
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass, mass_surf, mass_wat, mass_sol
        CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: symbol, symbol_surf, symbol_wat, symbol_sol
        REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: coor, coor_surf, coor_wat, coor_sol
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: &
          coor_wat_mol, &  ! 7: X, Y, Z, RR; angle_X, angle_Y, angle_Z
          coor_wat_mon     ! 4: r1, r2, r, theta
        INTEGER, DIMENSION(:), ALLOCATABLE :: atom1_array, atom2_array
        CHARACTER(len=80), DIMENSION(:), ALLOCATABLE :: filename12_array

        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: eps=1.0e-12_dp, eps2=1.0e-6_dp
        REAL(KIND=dp), PARAMETER :: bohr2ang=0.529177249_dp
        REAL(KIND=dp), PARAMETER :: au2debye=1.0_dp/0.39343029_dp
        REAL(KIND=dp), PARAMETER :: hart2kcal=627.51_dp, hart2cm=219474.63_dp, hart2kelvin=315773.3_dp
        REAL(KIND=dp), PARAMETER :: fs2ps=1.0e-3_dp, autime2fs=2.41888432650478e-2_dp

        END MODULE common_data_mod
