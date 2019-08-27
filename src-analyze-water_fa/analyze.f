        PROGRAM analyze

        USE common_data_mod
        USE binning_energy_mod
        USE reading_mod
        USE writing_all_mod
        USE pos_mod
        USE pos_sol_mod
        USE binning_mod
        USE dist_mod
        USE writing_dist_mod

! Purpose: Analyze a CP2K simulation of [(H20)n) slab] + [single adsorbate = sol = solute = fa = trans formic acid here]
! Program is run as:
!       analyze.exe < analyze.in >& analyze.out &
!
! It is assumed that the file slab-pos-1.xyz has the atoms organized as O,H,H (WCs); O,H,H (WCs); ...; Cl,H (WCs)
!
! Program reads first the file analyze.in:
! 10000         f_avg   final f frames are used for averaging (e.g., avg calculated from last 5ps of the run)
! 20000         f       total number of frames, not including the zeroth frame (total run=0.5fs*20k=10ps)
! 0.5           dt      time step [fs]
! 72            m_wat   number of water molecules
! 3             a_wat   number of atoms per water molecule (increase if have Wannier Centers or off-atomic sites in classical potential simulations)
! 1             m_sol   number of solute molecules
! 5             a_sol   number of atoms per solute molecule (this is for formic acid, change for differnt molecule--for, HCl will be 2; also increase if have Wannier Centers or off-atomic sites in classical potential simulations)
! 13.4724 15.5566 40.0              cell(d)   ABC values [Angstrom]
! 300.0                             temperature [Kelvin]
! -1237.4989259870 -38.7939203998   pot_a pot_b   pot_water_slab (equilibrated NVT sim of isolated water slab); pot_adsorbate energies of the adsorbate/solute (say from a geopt) [Hartree]
!                                                 (or: total energy of the identical periodic geo opt systems)
! 101           c   number of bins = c + 1
! 1                                 n_dist   number of distances
! 217 218 O217_H218_dist.dat        atom1_array(jn), atom2_array(jn), filename12_array(jn)   atom numbering based on vmd convention: 0, 1, ...
                                 !  some distances are already calculated here anyway by default, so this could be for additional distnces when Grotthus mechanism is active                                   
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

        REAL(KIND=dp) :: pi, rad2deg, deg2rad

        pi = ACOS(-1.0_dp)
        rad2deg = 180.0_dp/pi
        deg2rad = 1.0_dp/rad2deg

        WRITE(6,'(a)') '# PROGRAM = analyze'
        CALL FLUSH(6)

        READ(5,*) f_avg
        READ(5,*) f
        f0 = f - f_avg
        READ(5,*) dt
        READ(5,*) m_wat
        READ(5,*) a_wat
        READ(5,*) m_sol
        READ(5,*) a_sol
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
        WRITE(6,100) '# f0 = ', f0
        WRITE(6,200) '# dt = ', dt
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
        ALLOCATE( coor(0:f,a,my_d) )
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
        my_file='coor_wat.dat'
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

! coor_wat_mon(jf,jm,3) = r ! r = (r1+r2)/2
! coor_wat_mon(jf,jm,4) = theta

        my_file='r_wat.bin'
        my_m=m_wat
        my_d=3
        CALL binning(my_m, my_d, coor_wat_mon, &
          my_file)

        my_file='theta_wat.bin'
        my_m=m_wat
        my_d=4
        CALL binning(my_m, my_d, coor_wat_mon, &
          my_file)

        DEALLOCATE( coor_wat_mon )

! selected interatomic distances
        CALL dist

        DEALLOCATE( atom1_array, atom2_array, filename12_array )

! position magnitude and orientation of solute in the lab frame and
! solute monomer geometry
        my_file='coor_sol.dat'
        my_m=m_sol
        my_a=a_sol
        my_d1=d
        my_d2=10 ! XYZ and RR of C:four; and six angles
        ALLOCATE( coor_sol_mol(0:f, my_m, my_d2) )
! solute monomer geometry
! 2 3 O217_H218 O-H
! 1 4 C216_O219 C=O
! 1 5 C216_H220 C-H
! 1 2 C216_O217 C-O
! 3 2 1 H218_O217_C216 H-O-C
! 2 1 4 O217_C216_O219 O-C=O
! 5 1 4 H220_C216_O219 H-C=O
! dihedral angle
        my_d3=8 ! above four bonds and three angles and one dihedral angle
        ALLOCATE( coor_sol_mon(0:f, my_m, my_d3) )
        CALL pos_sol(my_m, my_a, my_d1, my_d2, my_d3, coor_sol, &
          my_file)

        my_m=m_sol

! coor_sol_mon(jf,jm,1) = roh
        my_d=1
        my_file='roh_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='O217_H218.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,2) = rc2o
        my_d=2
        my_file='rc2o_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='C216_O219.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,3) = rch
        my_d=3
        my_file='rch_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='C216_H220.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,4) = rco
        my_d=4
        my_file='rco_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='C216_O217.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,5) = ahoc
        my_d=5
        my_file='ahoc_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='H218_O217_C216.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,6) = aoc2o
        my_d=6
        my_file='aoc2o_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='O217_C216_O219.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,7) = ahc2o
        my_d=7
        my_file='ahc2o_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)
!       my_file='H220_C216_O219.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mon(jf,jm,8) = thcoh
        my_d=8
!       my_file='dihedral.dat'
!       CALL writing_dist(my_m, my_d, coor_sol_mon, my_file)
        coor_sol_mon(:,:,my_d) = ABS( coor_sol_mon(:,:,my_d) ) ! makes more sense to bin the absolute value of the dihedral angle
        my_file='thcoh_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mon, my_file)

! coor_sol_mol(jf,jm,3) = Z
        my_d=3
        my_file='Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

! coor_sol_mol(jf,jm,5) = aoh_Z
        my_d=5
        my_file='aoh_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)
!
        coor_sol_mol(:,:,my_d) = COS( coor_sol_mol(:,:,my_d)*deg2rad )
        my_file='cos_aoh_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

! coor_sol_mol(jf,jm,6) = ach_Z
        my_d=6
        my_file='ach_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)
!
        coor_sol_mol(:,:,my_d) = COS( coor_sol_mol(:,:,my_d)*deg2rad )
        my_file='cos_ach_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

! coor_sol_mol(jf,jm,7) = aco_Z
        my_d=7
        my_file='aco_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)
!
        coor_sol_mol(:,:,my_d) = COS( coor_sol_mol(:,:,my_d)*deg2rad )
        my_file='cos_aco_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

! coor_sol_mol(jf,jm,8) = an_Z
        my_d=8
        my_file='an_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)
!
        coor_sol_mol(:,:,my_d) = COS( coor_sol_mol(:,:,my_d)*deg2rad )
        my_file='cos_an_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

! coor_sol_mol(jf,jm,9) = ac2o_Z
        my_d=9
        my_file='ac2o_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)
!
        coor_sol_mol(:,:,my_d) = COS( coor_sol_mol(:,:,my_d)*deg2rad )
        my_file='cos_ac2o_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

!  coor_sol_mol(jf,jm,10) = abisector_Z
        my_d=10
        my_file='abisector_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)
!
        coor_sol_mol(:,:,my_d) = COS( coor_sol_mol(:,:,my_d)*deg2rad )
        my_file='cos_abisector_Z_sol.bin'
        CALL binning(my_m, my_d, coor_sol_mol, my_file)

        DEALLOCATE( coor_sol_mol )
        DEALLOCATE( coor_sol_mon )

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
