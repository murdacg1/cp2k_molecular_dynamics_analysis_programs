        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8
        INTEGER, PARAMETER :: d=3

        CHARACTER(len=80) :: in_file, out_file, title1, title2
        INTEGER :: a, nx, ny, nz
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: origin
        INTEGER, DIMENSION(:), ALLOCATABLE :: jatomicnum
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: atomicnumeff
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: h
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: coor
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: cube
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: dat
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: z

        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: eps=1.0e-12_dp, eps2=1.0e-6_dp
        REAL(KIND=dp), PARAMETER :: bohr2ang=0.529177249_dp
        REAL(KIND=dp), PARAMETER :: au2debye=1.0_dp/0.39343029_dp
        REAL(KIND=dp), PARAMETER :: hart2kcal=627.51_dp, hart2cm=219474.63_dp, hart2kelvin=315773.3_dp
        REAL(KIND=dp), PARAMETER :: fs2ps=1.0e-3_dp, autime2fs=2.41888432650478e-2_dp

        END MODULE common_data_mod
