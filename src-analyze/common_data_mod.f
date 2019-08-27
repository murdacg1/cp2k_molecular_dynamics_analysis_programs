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
          coor, velo, forc, chrg, dipl, quad, nddo
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: com, velo_com
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: &
          coor_mol, &  ! 7: X, Y, Z, RR; angle_X, angle_Y, angle_Z
          coor_mon, &  ! 4: r1, r2, r, theta
          velo_mol, &  ! 7: vX, vY, vZ, vv; angle_vX, angle_vY, angle_vZ
          forc_mol, &  ! 7: fX, fY, fZ, ff; angle_fX, angle_fY, angle_fZ
          chrg_mol, &
          dipl_mol, &
          dipl_sys, &  ! total system dipole rel to oxygens
          dipl_cell, & ! total system dipole rel to cell
          quad_mol, &  ! quad_mol(0:f,m,7),
          nddo_o, &
          nddo_h       ! scp-nddo restraint: 1=ptot, 2=f
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          vacf_mol     ! vacf of molecule=oxygen (assumed to be same as center of mass)
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          msd_mol      ! msd of molecule=oxygen (assumed to be same as center of mass)

        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: eps=1.0e-12_dp
        REAL(KIND=dp), PARAMETER :: bohr2ang=0.529177249_dp
        REAL(KIND=dp), PARAMETER :: au2debye=1.0_dp/0.39343029_dp
        REAL(KIND=dp), PARAMETER :: hart2kcal=627.51_dp, hart2cm=219474.63_dp, hart2kelvin=315773.3_dp
        REAL(KIND=dp), PARAMETER :: fs2ps=1.0e-3_dp, autime2fs=2.41888432650478e-2_dp

        END MODULE common_data_mod
