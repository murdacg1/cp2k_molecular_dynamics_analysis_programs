        PROGRAM analyze

        USE common_data_mod
        USE binning_energy_mod
        USE reading_mod
        USE writing_all_mod
        USE pos_mod
        USE binning_mod
        USE dist_mod

! Purpose: Analyze a CP2K simulation of [(H20)n) slab] + [single adsorbate]
! Program is run as:
!       analyze.exe < analyze.in >& analyze.out &
!
! It is assumed that the file slab-pos-1.xyz has the atoms organized as O,H,H (WCs); O,H,H (WCs); ...; Cl,H (WCs)
!
! Program reads first the file analyze.in:
! 25000         f_avg   final f frames are used for averaging (avg calculated from last 5ps of the run)
! 25000         f       total number of frames, not including the zeroth frame (total run=0.4fs*25k=10ps)
! 0.4           dt      time step [fs]
! 72            m_wat   number of water molecules
! 3             a_wat   number of atoms per water molecule (increase if have WCs)
! 1             m_sol   number of solute molecules
! 5             a_sol   number of atoms per solute molecule (increase if have WCs)
! 13.4724 15.5566 40.0              cell(d)   ABC values [Angstrom]
! 300.0                             temperature [Kelvin]
! -1237.4989259870 -38.7939203998   pot_a pot_b   pot_water_slab pot_adsorbate energies of the monomers [Hartree]
!                                                 (total energy of the identical peridic geo opt systems)
! 101           c   number of bins = c + 1
! 12                                n_dist   number of distances
! 120 121 Cl1_H2.dat                atom1_array(jn), atom2_array(jn), filename12_array(jn)   atom numbering based on vmd convention: 0, 1, ...
! 120 111 Cl1_O3.dat
! 121 111 H2_O3.dat
! 111 113 O3_H4.dat
! 111 99 O3_O5.dat
! 113 99 H4_O5.dat
! 111 112 O3_H6.dat
! 111 65 O3_O7.dat
! 99 100 O5_H8.dat
! 99 102 O5_O10.dat
! 120 103 Cl1_H11.dat
! 120 102 Cl1_O10.dat

        IMPLICIT NONE

        REAL(KIND=dp) :: pot_a, pot_b

        INTEGER :: &
          my_m, &
          my_a, &
          my_d, my_d1, my_d2, my_d3, &
          jn
        CHARACTER(len=80) my_file
        REAL(KIND=dp) :: &
          my_pot_a, my_pot_b, &
          my_factor

        WRITE(6,'(a)') '# PROGRAM = analyze'
        CALL FLUSH(6)

        READ(5,*) f_avg
        READ(5,*) f
        READ(5,*) dt
        READ(5,*) m_wat
        READ(5,*) a_wat
        READ(5,*) m_sol
        READ(5,*) a_sol
        m = 1 ! total number molecules
        a = m_wat*a_wat + m_sol*a_sol ! total number atoms
        ALLOCATE( cell(d) )
        cell = zero
        READ(5,*) cell(:)
        READ(5,*) temperature
        READ(5,*) pot_a, pot_b
        READ(5,*) c
        READ(5,*) n_dist
        ALLOCATE( atom1_array(n_dist), atom2_array(n_dist), filename12_array(n_dist) )
        DO jn = 1, n_dist
          READ(5,*) atom1_array(jn), atom2_array(jn), filename12_array(jn)
        END DO

        WRITE(6,100) '# f_avg = ', f_avg
        WRITE(6,100) '# f = ', f
        WRITE(6,200) '# dt = ', dt
        WRITE(6,100) '# m = ', m
        WRITE(6,100) '# a = ', a
        WRITE(6,100) '# m_wat = ', m_wat
        WRITE(6,100) '# a_wat = ', a_wat
        WRITE(6,100) '# m_sol = ', m_sol
        WRITE(6,100) '# a_sol = ', a_sol
        WRITE(6,100) '# d = ', d
        WRITE(6,300) '# cell = ', cell(:)
        WRITE(6,200) '# temperature [Kelvin] = ', temperature
        WRITE(6,210) '# pot_a, pot_b [Hartree] = ', pot_a, pot_b
        WRITE(6,100) '# number of bins c = ', c
        WRITE(6,100) '# number of distances n_dist = ', n_dist
        WRITE(6,'(a)') &
          '# atom1_array(jn), atom2_array(jn), filename12_array(jn):'
        DO jn = 1, n_dist
          WRITE(6,*) &
             atom1_array(jn), atom2_array(jn), filename12_array(jn)
        END DO
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( tim(0:f), pot(0:f) )
        ALLOCATE( symbol(a), symbol_wat(a), symbol_sol(a) )
        ALLOCATE( mass(a), mass_wat(a), mass_sol(a) )

!   Trajectory: read in each frame

        my_file='../slab-pos-1.xyz'
        my_d=d
        ALLOCATE( coor(0:f,m,a,my_d) )
        ALLOCATE( coor_wat(0:f,m_wat,a_wat,my_d) )
        ALLOCATE( coor_sol(0:f,m_sol,a_sol,my_d) )
        CALL reading(my_d, coor, coor_wat, coor_sol, &
          my_file)

! bin pot
        my_pot_a=zero
        my_pot_b=zero
        my_factor=1.0_dp
        my_file='pot.bin'
        CALL binning_energy(pot, &
          my_pot_a, my_pot_b, my_factor, &
          my_file)

! bin ecoh
        my_pot_a=pot_a
        my_pot_b=pot_b
        my_factor=hart2kcal
        my_file='ecoh.bin'
        CALL binning_energy(pot, &
          my_pot_a, my_pot_b, my_factor, &
          my_file)

! position magnitude and orientation of molecule=oxygen and
! monomer geo r1, r2, r=(r1+r2)/2, theta
        my_file='./coor_wat.dat'
        my_m=m_wat
        my_a=a_wat
        my_d1=d
        my_d2=7
        my_d3=4
        ALLOCATE( coor_wat_mol(0:f, my_m, my_d2) )
        ALLOCATE( coor_wat_mon(0:f, my_m, my_d3) )
        CALL pos(my_m, my_a, my_d1, my_d2, my_d3, coor_wat, &
          my_file)
        DEALLOCATE( coor_wat_mol )

! bin r_mon = r = (r1+r2)/2
        my_file='r_wat.bin'
        my_m=m_wat
        my_d=3
        CALL binning(my_m, my_d, coor_wat_mon, &
          my_file)

! bin theta_mon
        my_file='theta_wat.bin'
        my_m=m_wat
        my_d=4
        CALL binning(my_m, my_d, coor_wat_mon, &
          my_file)

        DEALLOCATE( coor_wat_mon )

! selected interatomic distances
        CALL dist

        DEALLOCATE( atom1_array, atom2_array, filename12_array )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( cell )
        DEALLOCATE( tim, pot )
        DEALLOCATE( symbol, symbol_wat, symbol_sol)
        DEALLOCATE( mass, mass_wat, mass_sol )
        DEALLOCATE( coor, coor_wat, coor_sol )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,1(1x,e20.10))
210     FORMAT(1x,a40,2(1x,e20.10))
300     FORMAT(1x,a40,3(1x,e20.10))
410     FORMAT(1x,a40,4(1x,e20.10))

        END PROGRAM analyze
