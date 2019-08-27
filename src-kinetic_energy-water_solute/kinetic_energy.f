        PROGRAM kinetic_energy

        USE common_data_mod
        USE reading_and_calc_ke_mod
        USE binning_mod

! Purvelocitye: Analyze a CP2K simulation of [(H20)n) slab] + [single adsorbate = sol = solute = fa = trans formic acid here]
! Program is run as:
!       kinetic_energy.exe < kinetic_energy.in >& kinetic_energy.out &
!
! It is assumed that the file slab-vel-1.xyz has the atoms organized as O,H,H (WCs); O,H,H (WCs); ...; Cl,H (WCs)
!
! Program reads first the file kinetic_energy.in:
! 10000         f_avg   final f frames are used for averaging (e.g., avg calculated from last 5ps of the run)
! 20000         f       total number of frames, not including the zeroth frame (total run=0.5fs*20k=10ps)
! 72            m_wat   number of water molecules
! 3             a_wat   number of atoms per water molecule (increase if have Wannier Centers or off-atomic sites in classical potential simulations)
! 1             m_sol   number of solute molecules
! 5             a_sol   number of atoms per solute molecule (this is for formic acid, change for differnt molecule--for, HCl will be 2; also increase if have Wannier Centers or off-atomic sites in classical potential simulations)
! 300.0         temperature   [Kelvin]
! 101           c   number of bins = c + 1

        IMPLICIT NONE

        INTEGER :: my_m, my_d
        CHARACTER(len=80) my_file

        REAL(KIND=dp) :: pi, rad2deg, deg2rad

        pi = ACOS(-1.0_dp)
        rad2deg = 180.0_dp/pi
        deg2rad = 1.0_dp/rad2deg

        WRITE(6,'(a)') '# PROGRAM = kinetic_energy'
        CALL FLUSH(6)

        READ(5,*) f_avg
        READ(5,*) f
        f0 = f - f_avg
        READ(5,*) m_wat
        READ(5,*) a_wat
        READ(5,*) m_sol
        READ(5,*) a_sol
        a = m_wat*a_wat + m_sol*a_sol ! total number atoms
        READ(5,*) temperature
        READ(5,*) c

        WRITE(6,100) '# f_avg = ', f_avg
        WRITE(6,100) '# f = ', f
        WRITE(6,100) '# f0 = ', f0
        WRITE(6,100) '# a = ', a
        WRITE(6,100) '# m_wat = ', m_wat
        WRITE(6,100) '# a_wat = ', a_wat
        WRITE(6,100) '# m_sol = ', m_sol
        WRITE(6,100) '# a_sol = ', a_sol
        WRITE(6,100) '# a = ', a
        WRITE(6,100) '# d3, d1 = ', d3, d1
        WRITE(6,200) '# temperature [Kelvin] = ', temperature
        WRITE(6,100) '# number of bins c = ', c
        WRITE(6,*) &
          '# bohr2ang, fs2ps, autime2fs, au_velocity_to_ang_per_ps = ', &
             bohr2ang, fs2ps, autime2fs, au_velocity_to_ang_per_ps
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( tim(0:f) )
        ALLOCATE( symbol(a), symbol_wat(a_wat), symbol_sol(a_sol) )
        ALLOCATE( mass(a), mass_wat(a_wat), mass_sol(a_sol) )

!   Trajectory: read in each frame and calculate ke

        my_file='../slab-vel-1.xyz'

        ALLOCATE( velo(0:f,a,d3) )
        ALLOCATE( velo_wat(0:f,m_wat,a_wat,d3) )
        ALLOCATE( velo_sol(0:f,m_sol,a_sol,d3) )

        ALLOCATE( ke(0:f,1,d1) )
        ALLOCATE( ke_slab(0:f,1,d1) )
        ALLOCATE( ke_wat(0:f,m_wat,d1) )
        ALLOCATE( ke_sol(0:f,m_sol,d1) )            ! ke_sol is total kinetic energy of the solute
        ALLOCATE( ke_com_sol(0:f,m_sol,d4) )        ! ke_com_sol is kinetic energy of the com of the solute, also get its components
        ALLOCATE( ke_vibrot_sol(0:f,m_sol,d1) )     ! ke_vibrot_sol is vib and rot kinetic energy of the solute
        ALLOCATE( u_com_sol(0:f,m_sol,d4) )         ! ke_com_sol is kinetic energy of the com of the solute, also get its components

        CALL reading_and_calc_ke(velo, velo_wat, velo_sol, &
          ke, ke_slab, ke_wat, ke_sol, ke_com_sol, ke_vibrot_sol, u_com_sol, &
          my_file)

        my_file='ke.bin'                              ! ke of whole simulation cell (can be compared to .ener file)
        my_m=1
        my_d=d1
        CALL binning(my_m, my_d, ke, &
          my_file)

        my_file='ke_slab.bin'                         ! water slab
        my_m=1
        my_d=d1
        CALL binning(my_m, my_d, ke_slab, &
          my_file)

        my_file='ke_wat.bin'                          ! all waters; should be equal to water slab
        my_m=m_wat
        my_d=d1
        CALL binning(my_m, my_d, ke_wat, &
          my_file)

        my_file='ke_sol.bin'                          ! total ke of solute
        my_m=m_sol
        my_d=d1
        CALL binning(my_m, my_d, ke_sol, &
          my_file)

        my_file='ke_com_sol.bin'                      ! ke: com trans of solute
        my_m=m_sol
        my_d=d4
        CALL binning(my_m, my_d, ke_com_sol, &
          my_file)

        my_file='ke_vibrot_sol.bin'                   ! vib and rot ke of solute
        my_m=m_sol
        my_d=d1
        CALL binning(my_m, my_d, ke_vibrot_sol, &
          my_file)

        my_file='ke_com_sol_x.bin'                    ! ke: com trans of solute, x-component 
        my_m=m_sol
        my_d=1
        CALL binning(my_m, my_d, ke_com_sol, &
          my_file)

        my_file='ke_com_sol_y.bin'                    ! ke: com trans of solute, y-component
        my_m=m_sol
        my_d=2
        CALL binning(my_m, my_d, ke_com_sol, &
          my_file)

        my_file='ke_com_sol_z.bin'                    ! ke: com trans of solute, z-component
        my_m=m_sol
        my_d=3
        CALL binning(my_m, my_d, ke_com_sol, &
          my_file)

        my_file='u_com_sol.bin'                      ! speed u: com trans of solute
        my_m=m_sol
        my_d=d4
        CALL binning(my_m, my_d, u_com_sol, &
          my_file)

        my_file='u_com_sol_x.bin'                    ! speed u: com trans of solute, x-component 
        my_m=m_sol
        my_d=1
        CALL binning(my_m, my_d, u_com_sol, &
          my_file)

        my_file='u_com_sol_y.bin'                    ! speed u: com trans of solute, y-component
        my_m=m_sol
        my_d=2
        CALL binning(my_m, my_d, u_com_sol, &
          my_file)

        my_file='u_com_sol_z.bin'                    ! speed u: com trans of solute, z-component
        my_m=m_sol
        my_d=3
        CALL binning(my_m, my_d, u_com_sol, &
          my_file)

        DEALLOCATE( ke )
        DEALLOCATE( ke_slab )
        DEALLOCATE( ke_wat )
        DEALLOCATE( ke_sol )
        DEALLOCATE( ke_com_sol )
        DEALLOCATE( ke_vibrot_sol )
        DEALLOCATE( u_com_sol )

        DEALLOCATE( tim )
        DEALLOCATE( symbol, symbol_wat, symbol_sol)
        DEALLOCATE( mass, mass_wat, mass_sol )
        DEALLOCATE( velo, velo_wat, velo_sol )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,1(1x,e20.10))

        END PROGRAM kinetic_energy
