        PROGRAM msite

        USE common_data_mod
        USE reading_mod
        USE get_msite_mod
        USE writing_mod

! Purpose: Determine and write Dang-Change potential M site for CP2K simulation of (H20)n

! Program is run as:
!       ../msite.exe < msite.in >& msite.out &
!
! It is assumed that the in file water_no_msite.xyz  will contain some frames organized as O,H,H; O,H,H; ...
!          and that the out file water_yes_msite.xyz will contain some frames organized as O,H,H,M; O,H,H,M; ...

! Program reads first the file msite.in:

! 0         f_avg (the final f frames are used for averaging) (avg calculated from last 5ps of the run)
! 0         f (total number of frames, not including the zeroth frame) (total run=0.5fs*20k=10ps)
! 216           m (number of molecules per frame or cell)
! 3             a (number of atoms per molecule)
! 15.0 15.0 71.44  cell ABC values [Angstrom]
! 0.215            dist of M site from O site along the bisector [Angstrom]
! 15.99491463 1.0078250321 1.0078250321 0.0000000000 mass1 mass2 mass3 mass4

        IMPLICIT NONE

        REAL(KIND=dp) :: mass1, mass2, mass3, mass4

        INTEGER :: my_d
        CHARACTER(len=80) my_file

        WRITE(6,'(a)') '# PROGRAM = msite'
        CALL FLUSH(6)

        ALLOCATE( cell(d) )

        READ(5,*) f_avg
        READ(5,*) f
        READ(5,*) m
        READ(5,*) a
        READ(5,*) cell(:)
        READ(5,*) dist
        READ(5,*) mass1, mass2, mass3, mass4

        mass = mass1 + mass2 + mass3 + mass4

        WRITE(6,200) '# f_avg, f = ', f_avg, f
        WRITE(6,100) '# m = ', m
        WRITE(6,100) '# a = ', a
        WRITE(6,310) '# cell = ', cell(:)
        WRITE(6,110) '# dist = ', dist
        WRITE(6,510) '# mass1, mass2, mass3, mass4, mass [amu] = ', &
          mass1, mass2, mass3, mass4, mass

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        my_d = d

        ALLOCATE( symbol(m,a+1) )
        ALLOCATE( molecule(m,a+1) )
        ALLOCATE( coor(0:f,m,a+1,my_d) )

        my_file = './water_no_msite.xyz'
        my_file = TRIM(my_file)
        CALL reading(my_d, coor, my_file)

        CALL get_msite

        my_file = './water_yes_msite.xyz'
        my_file = TRIM(my_file)
        CALL writing(my_d, coor, my_file)

        DEALLOCATE( cell )
        DEALLOCATE( coor )
        DEALLOCATE( symbol )
        DEALLOCATE( molecule )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,2(1x,i10))
110     FORMAT(1x,a40,1(1x,e20.10))
310     FORMAT(1x,a40,3(1x,e20.10))
510     FORMAT(1x,a40,5(1x,e20.10))

        END PROGRAM msite
