        MODULE common_data_mod

        INTEGER, PARAMETER :: dp=8

        INTEGER :: f_avg, f, f0, &
          m_wat, a_wat, &
          m_sol, a_sol, &
          a, &
          c
        INTEGER, PARAMETER :: d1=1, d3=3, d4=4
        REAL(KIND=dp) :: temperature
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tim
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass, mass_wat, mass_sol
        CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: symbol, symbol_wat, symbol_sol
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: velo
        REAL(KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: velo_wat, velo_sol
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: &
          ke, ke_slab, ke_wat, ke_sol, ke_com_sol, ke_vibrot_sol, u_com_sol
        REAL(KIND=dp), PARAMETER :: zero=0.0_dp
        REAL(KIND=dp), PARAMETER :: eps=1.0e-12_dp, eps2=1.0e-6_dp
        REAL(KIND=dp), PARAMETER :: bohr2ang=0.529177249_dp
        REAL(KIND=dp), PARAMETER :: au2debye=1.0_dp/0.39343029_dp
        REAL(KIND=dp), PARAMETER :: hart2kcal=627.51_dp, hart2cm=219474.63_dp
        REAL(KIND=dp), PARAMETER :: kcal2kJ=4.184_dp
!       REAL(KIND=dp), PARAMETER :: hart2kelvin=315773.3_dp
        REAL(KIND=dp), PARAMETER :: hart2kelvin=315774.647902944_dp ! cp2k internal values
        REAL(KIND=dp), PARAMETER :: fs2ps=1.0e-3_dp, autime2fs=2.41888432650478e-2_dp ! cp2k internal values
        REAL(KIND=dp), PARAMETER :: m_u2m_amu=1822.88848426455_dp ! cp2k internal values

        REAL(KIND=dp), PARAMETER :: au_velocity_to_ang_per_ps = bohr2ang/(autime2fs*fs2ps)

!       REAL(KIND=dp), PARAMETER :: mass_O=15.99491463_dp
!       REAL(KIND=dp), PARAMETER :: mass_H=1.0078250321_dp
!       REAL(KIND=dp), PARAMETER :: mass_C=12.0_dp
!       REAL(KIND=dp), PARAMETER :: mas_Cl=34.968852721_dp
!       REAL(KIND=dp), PARAMETER :: mass_Si=27.9769271_dp

! for best agreement with .ener file, use cp2k internal masses, see, eg, https://doxygen.cp2k.org/d5/de9/periodic__table_8f90_source.html

        REAL(KIND=dp), PARAMETER :: mass_O=15.9994_dp
        REAL(KIND=dp), PARAMETER :: mass_H=1.00794_dp
        REAL(KIND=dp), PARAMETER :: mass_C=12.0107_dp
        REAL(KIND=dp), PARAMETER :: mas_Cl=35.453_dp
        REAL(KIND=dp), PARAMETER :: mass_N=14.0067_dp
        REAL(KIND=dp), PARAMETER :: mass_M=zero
        REAL(KIND=dp), PARAMETER :: mass_S=32.065_dp
        REAL(KIND=dp), PARAMETER :: mass_Si=28.0855_dp

        END MODULE common_data_mod
