        PROGRAM analyze

        USE common_data_mod
        USE reading_mod
        USE minus_com_mod
        USE writing_all_mod

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
          my_m, &
          my_a, &
          my_d
        CHARACTER(len=1) my_symbol
        CHARACTER(len=80) my_file
        REAL(KIND=dp) :: &
          my_pot_mon, &
          my_factor

        WRITE(6,'(a)') '# PROGRAM = analyze'
        CALL FLUSH(6)

        READ(5,*) f_avg
        READ(5,*) f
        READ(5,*) m
        READ(5,*) a
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

!       GOTO 10

!   Trajectory: read in and subtract off COM position in each frame and write out COM-less data

        my_file='../run-01.xyz'
        my_d=d
        ALLOCATE( coor(0:f,m,a,my_d) )
        CALL reading(my_d, coor, &
          my_file)

        ALLOCATE( com(0:f,d+1) ) ! position of com
        my_file='./run-01_com.xyz'
        CALL minus_com(coor, com, &
          my_file)
        DEALLOCATE( com )

        my_file='./run-01.xyz'
        my_m=m
        my_a=a
        my_d=d
        CALL writing_all(my_m, my_a, my_d, coor, &
     &    my_file)

        DEALLOCATE( coor )

!10      CONTINUE

!   Velocity: read in and subtract off COM velocity in each frame and write out COM-less data

        my_file='../run-01_vel.xyz'
        my_d=d
        ALLOCATE( velo(0:f,m,a,my_d) )
        CALL reading(my_d, velo, &
          my_file)

        ALLOCATE( velo_com(0:f,d+1) ) ! position of com
        my_file='./run-01_vel_com.xyz'
        CALL minus_com(velo, velo_com, &
          my_file)
        DEALLOCATE( velo_com )

        my_file='./run-01_vel.xyz'
        my_m=m
        my_a=a
        my_d=d
        CALL writing_all(my_m, my_a, my_d, velo, &
     &    my_file)

        DEALLOCATE( velo )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( cell )
        DEALLOCATE( mass )
        DEALLOCATE( tim )
        DEALLOCATE( pot )
        DEALLOCATE( symbol )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,1(1x,e20.10))
300     FORMAT(1x,a40,3(1x,e20.10))
410     FORMAT(1x,a40,4(1x,e20.10))

        END PROGRAM analyze
