        PROGRAM density

        USE common_data_mod
        USE reading_mod
        USE get_density_mod

! Purpose: density for CP2K simulation of (H20)n
! Program is run as:
!       ../density.exe < density.in >& density.out &
!
! It is assumed that the files run-01*.xyz have the atoms organized as O,H,H; O,H,H; ...
!
! Program reads first the file analyze_run.in:
! 20000         f_avg (the final f frames are used for averaging) (avg calculated from last 5ps of the run)
! 29531         f (total number of frames, not including the zeroth frame) (total run=0.5fs*20k=10ps)
! 216           m (number of molecules per frame or cell)
! 3             a (number of atoms per molecule)
! 15.0 15.0 71.44  cell ABC values [Angstrom]
! 0.10             dz [Angstrom]
! 15.99491463 1.0078250321 1.0078250321  mass1 mass2 mass3

        IMPLICIT NONE

        REAL(KIND=dp) :: mass1, mass2, mass3
        CHARACTER(len=1) name1
        CHARACTER(len=80) my_file

        WRITE(6,'(a)') '# PROGRAM = density'
        CALL FLUSH(6)

        ALLOCATE( cell(d) )

        READ(5,*) f_avg
        READ(5,*) f
        READ(5,*) m
        READ(5,*) a
        READ(5,*) cell(:)
        READ(5,*) dz
        READ(5,*) mass1, mass2, mass3

        b = int( cell(d) / dz )
        mass = mass1 + mass2 + mass3

        WRITE(6,200) '# f_avg, f = ', f_avg, f
        WRITE(6,100) '# m = ', m
        WRITE(6,100) '# a = ', a
        WRITE(6,310) '# cell = ', cell(:)
        WRITE(6,310) '# dz = ', dz
        WRITE(6,100) '# b = ', b
        WRITE(6,410) '# mass1, mass2, mass3, mass [amu] = ', &
          mass1, mass2, mass3, mass

        ALLOCATE( coor(0:f,m,a,d) )
        ALLOCATE( symbol(m,a) )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        my_file='../run-01.xyz'
        my_file=TRIM(my_file)
        CALL reading(my_file)

        name1='O'
        my_file='density.dat'
        CALL get_density(name1, my_file)

        DEALLOCATE( cell )

        DEALLOCATE( coor )
        DEALLOCATE( symbol )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,2(1x,i10))
110     FORMAT(1x,a40,1(1x,e20.10))
310     FORMAT(1x,a40,3(1x,e20.10))
410     FORMAT(1x,a40,4(1x,e20.10))

        END PROGRAM density
