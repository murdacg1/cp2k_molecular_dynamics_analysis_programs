        PROGRAM analyze

        USE common_data_mod
        USE binning_energy_mod
        USE reading_mod
        USE writing_mod
        USE writing_all_mod
        USE pos_mod
        USE binning_mod
        USE info_mod

! Purpose: Analyze a CP2K simulation of (H20)n
! Program is run as:
!       analyze.exe < analyze.in >& analyze.out &
!
! It is assumed that the files run-01.xyz have the atoms organized as O,H,H; O,H,H; ...
!
! Program reads first the file analyze.in:
! 10000         f_avg (the final f frames are used for averaging) (avg calculated from last 5ps of the run)
! 20000         f (total number of frames, not including the zeroth frame) (total run=0.5fs*20k=10ps)
! 64            m (number of molecules per frame or cell)
! 3             a (number of atoms per molecule)
! 12.4138 12.4138 12.4138  cell(d) ABC values [Angstrom]
! 300.0		           temperature [Kelvin]
! -11.93952773617932       pot_mon, energy of the monomer [Hartree]
!                          (total energy of an isolated--cubic_box_length=30A--geometry-optimized
!                           water monomer at T=0K--geoopt calculation)
! 101           c (the number of bins = c + 1)
! 15.99491463 1.0078250321 1.0078250321  mass(a)

        IMPLICIT NONE

        REAL(KIND=dp) :: pot_mon

        INTEGER :: &
          my_m, my_m1, my_m2, &
          my_a, &
          my_d, my_d1, my_d2, my_d3, &
          my_d_chrg, my_d_dipl, my_m_dipl_cell
        CHARACTER(len=1) my_symbol
        CHARACTER(len=80) my_file
        REAL(KIND=dp) :: &
          my_pot_mon, &
          my_factor

        WRITE(6,'(a)') '# PROGRAM = analyze'
        CALL FLUSH(6)

        READ(5,*) f_avg
        READ(5,*) f
        f0 = f - f_avg
        READ(5,*) m
        READ(5,*) a
!       a_with_msite = a
        ALLOCATE( cell(d) )
        cell = zero
        READ(5,*) cell(:)
        READ(5,*) temperature
        READ(5,*) pot_mon
        READ(5,*) c
        ALLOCATE( mass(a) )
        mass = zero
        READ(5,*) mass(:)

        WRITE(6,100) '# f_avg = ', f_avg
        WRITE(6,100) '# f = ', f
        WRITE(6,100) '# f0+1 = ', f0+1
        WRITE(6,100) '# m = ', m
        WRITE(6,100) '# a = ', a
        WRITE(6,100) '# d = ', d
        WRITE(6,300) '# cell = ', cell(:)
        WRITE(6,200) '# temperature [Kelvin] = ', temperature
        WRITE(6,200) '# pot_mon [Hartree] = ', pot_mon
        WRITE(6,100) '# number of bins c = ', c
        WRITE(6,410) '# mass [amu] = ', mass(:), SUM(mass)
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( tim(0:f) )
        ALLOCATE( pot(0:f) )
        ALLOCATE( symbol(a) )

!   Trajectory: read in (no need to subtract COM position in each frame since calculating static quantity=NVT or re-thermalized=NVE )

        my_file='../run-01.xyz'
        my_d=d
        ALLOCATE( coor(0:f,m,a,my_d) )
        CALL reading(my_d, coor, &
          my_file)

! bin pot
        my_pot_mon=zero
        my_factor=1.0_dp
        my_file='pot.bin'
        my_m1=1
        my_m2=1
        CALL binning_energy(my_m1, my_m2, pot, &
          my_pot_mon, my_factor, &
          my_file)

! bin ecoh
        my_pot_mon=pot_mon
        my_factor=hart2kcal
        my_file='ecoh.bin'
        my_m1=1
        my_m2=m
        CALL binning_energy(my_m1, my_m2, pot, &
          my_pot_mon, my_factor, &
          my_file)

! position magnitude and orientation of molecule=oxygen and
! monomer geo r1, r2, r=(r1+r2)/2, theta
        my_file='coor.dat'
        my_d1=d
        my_d2=7
        my_d3=4
        ALLOCATE( coor_mol(0:f,m,my_d2) )
        ALLOCATE( coor_mon(0:f,m,my_d3) )
        CALL pos(my_d1, my_d2, my_d3, &
          my_file)

! msd_mol of molecule=oxygen
        my_file='msd_mol.dat'
        ALLOCATE ( msd_mol(0:f) )
        CALL msd(my_file)

        DEALLOCATE( coor_mol )

! d_msd_mol=self-diffusion coefficient from msd of molecule=oxygen
        CALL d_msd

        DEALLOCATE( msd_mol )

! bin r_mon = r = (r1+r2)/2
        my_file='r_mon.bin'
        my_d=3
        CALL binning(my_d, coor_mon, &
          my_file)

! bin theta_mon
        my_file='theta_mon.bin'
        my_d=4
        CALL binning(my_d, coor_mon, &
          my_file)

        DEALLOCATE( coor_mon )

!   Velocity (no need to subtract COM position in each frame since calculating static quantity=NVT or re-thermalized=NVE )

        my_file='../run-01_vel.xyz'
        my_d=d
        ALLOCATE( velo(0:f,m,a,my_d) )
        CALL reading(my_d, velo, &
          my_file)

! velocity magnitude and direction of molecule=oxygen
        my_file='velo_mol.dat'
        my_d1=3
        my_d2=7
        ALLOCATE( velo_mol(0:f,m,my_d2) )
        CALL info(my_d1, my_d2, velo, velo_mol, &
          my_file)

        DEALLOCATE( velo )

! bin velo_mol
        my_file='velo_mol.bin'
        my_d=4
        CALL binning(my_d, velo_mol, &
          my_file)

! vacf_mol=velocity autocorrelation function of molecule=oxygen
        my_file='vacf_mol.dat'
        ALLOCATE( vacf_mol(0:f) )
        CALL vacf(my_file)

        DEALLOCATE( velo_mol )

! d_vacf_mol=self-diffusion coefficient from velocity autocorrelation function of molecule=oxygen
        CALL d_vacf

        DEALLOCATE( vacf_mol )

!   Force
        my_file='../run-01_force.xyz'
        my_d=d
        ALLOCATE( forc(0:f,m,a,my_d) )
        CALL reading(my_d, forc, &
          my_file)

! force magnitude and direction of molecule=oxygen
        my_file='forc_mol.dat'
        my_d1=3
        my_d2=7
        ALLOCATE( forc_mol(0:f,m,my_d2) )
        CALL info(my_d1, my_d2, forc, forc_mol, &
          my_file)

        DEALLOCATE( forc )

! bin forc_mol
        my_file='forc_mol.bin'
        my_d=4
        CALL binning(my_d, forc_mol, &
          my_file)

        DEALLOCATE( forc_mol )

!   Charge
        my_file='../run-01_chrg.xyz'
        my_d_chrg=1
        ALLOCATE( chrg(0:f,m,a,my_d_chrg) )
        my_d=my_d_chrg
        CALL reading(my_d, chrg, &
          my_file)

! total charge of molecule=oxygen (should be CLOSE to 0.0au; charge transfer will modIFy this)
        my_file='chrg_mol.dat'
        ALLOCATE( chrg_mol(0:f,m,my_d_chrg) )
        CALL charge(my_d_chrg, &
          my_file)

! bin chrg_mol
        my_file='chrg_mol.bin'
        my_d=my_d_chrg
        CALL binning(my_d, chrg_mol, &
          my_file)

        DEALLOCATE( chrg_mol )

!   Dipole
        my_file='../run-01_dipl.xyz'
        my_d_dipl=d
        ALLOCATE( dipl(0:f,m,a,my_d_dipl) )
        my_d=my_d_dipl
        CALL reading(my_d, dipl, &
          my_file)

! total dipole of molecule=oxygen (should be CLOSE to the correct one for gas=1.86D or liquid=2.6-3.0D)
        my_file='dipl_mol.dat'
        my_d_dipl=4
        my_m_dipl_cell=1
        ALLOCATE( dipl_mol(0:f, m, my_d_dipl) )
        ALLOCATE( dipl_sys(0:f,my_m_dipl_cell,d) )
        ALLOCATE( dipl_cell(0:f,my_m_dipl_cell,d) )
        CALL dipole(my_d_chrg, my_d_dipl, my_m_dipl_cell, &
          my_file)

        DEALLOCATE( dipl )

        DEALLOCATE( cell )
        DEALLOCATE( mass )

        DEALLOCATE( coor )
        DEALLOCATE( chrg )

! DIPOLE_SYS
        my_file='DIPOLE.XYZ'
        my_m=my_m_dipl_cell
        my_d=d
        my_symbol='T' ! this is T=TOTAL system dipole
        CALL writing(my_m, my_d, dipl_sys, my_symbol, &
          my_file)

        DEALLOCATE( dipl_sys )

! DIPOLE_CELL
        my_file='DIPOLE_CELL.XYZ'
        my_m=my_m_dipl_cell
        my_d=d
        my_symbol='T' ! this is T=TOTAL system dipole
        CALL writing(my_m, my_d, dipl_cell, my_symbol, &
          my_file)

        DEALLOCATE( dipl_cell )

! bin dipl_mol
        my_file='dipl_mol.bin'
        my_d=4
        CALL binning(my_d, dipl_mol, &
          my_file)

        DEALLOCATE( dipl_mol )

!   Quadrupole
        my_file='../run-01_quad.xyz'
        my_d=6
        ALLOCATE( quad(0:f,m,a,my_d) )
        CALL reading(my_d, quad, &
          my_file)
        DEALLOCATE( quad )

!   SCP-NDDO restraint
        my_file='../run-01_nddo_restraint.xyz'
        my_d=2
        ALLOCATE( nddo(0:f,m,a,my_d) )
        CALL reading(my_d, nddo, &
          my_file)

! nddo array for oxygens
        my_file='nddo_o.dat'
        my_d=2
        my_factor=1.0_dp
        my_symbol='O'
        ALLOCATE( nddo_o(0:f,m,my_d) )
        CALL nddo_restraints(my_d, nddo, nddo_o, &
          my_factor, my_symbol, &
          my_file)

! bin nddo_ptot_o
        my_file='nddo_ptot_o.bin'
        my_d=1
        CALL binning(my_d, nddo_o, &
          my_file)

! bin nddo_f_o
        my_file='nddo_f_o.bin'
        my_d=2
        CALL binning(my_d, nddo_o, &
          my_file)

! nddo array for hydrogens
        my_file='nddo_h.dat'
        my_d=2
        my_factor=0.5_dp
        my_symbol='H'
        ALLOCATE( nddo_h(0:f,m,my_d) )
        CALL nddo_restraints(my_d, nddo, nddo_h, &
          my_factor, my_symbol, &
          my_file)

! bin nddo_ptot_h
        my_file='nddo_ptot_h.bin'
        my_d=1
        CALL binning(my_d, nddo_h, &
          my_file)

! bin nddo_f_h
        my_file='nddo_f_h.bin'
        my_d=2
        CALL binning(my_d, nddo_h, &
          my_file)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( pot )
        DEALLOCATE( tim )
        DEALLOCATE( symbol )

        DEALLOCATE( nddo )
        DEALLOCATE( nddo_o )
        DEALLOCATE( nddo_h )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,1(1x,e20.10))
300     FORMAT(1x,a40,3(1x,e20.10))
410     FORMAT(1x,a40,4(1x,e20.10))

        END PROGRAM analyze
