        PROGRAM hbond

        USE common_data_mod
        USE reading_mod
        USE get_hbond_mod

! Purpose: hbond for CP2K simulation of (H20)n
! Program is run as:
!       ../hbond.exe < hbond.in >& hbond.out &
!
! It is assumed that the files run-01*.xyz have the atoms organized as O,H,H; O,H,H; ...
!
! Program reads first the file hbond.in:
! 10000         f_avg (the final f frames are used for averaging) (avg calculated from last 5ps of the run)
! 20000         f (total number of frames, not including the zeroth frame) (total run=0.5fs*20k=10ps)
! 64            m (number of molecules per frame or cell)
! 3             a (number of atoms per molecule)
! 12.4138 12.4138 12.4138  cell ABC values [Angstrom]
! 2 2           aa dd (max number acceptors and donors)   

        IMPLICIT NONE

        CHARACTER(len=10) name1
        CHARACTER(len=80) my_file

        WRITE(6,'(a)') '# PROGRAM = hbond'
        CALL FLUSH(6)

        ALLOCATE( cell(d) )

        READ(5,*) f_avg
        READ(5,*) f
        READ(5,*) m
        READ(5,*) a
        READ(5,*) cell(:)
        READ(5,*) aa, dd

        WRITE(6,200) '# f_avg, f = ', f_avg, f
        WRITE(6,200) '# m, a = ', m, a
        WRITE(6,310) '# cell = ', cell(:)
        WRITE(6,200) '# aa, dd = ', aa, dd

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( coor(0:f,m,a,d) )
        ALLOCATE( symbol(m,a) )

        my_file = '../run-01.xyz'
        my_file = TRIM(my_file)
        WRITE(6,'(2a)') 'my_file = ', my_file
        CALL FLUSH(6)

        CALL reading(my_file)

        name1='160'
        my_file='hbond_' // TRIM(name1) // '.dat'
        CALL get_hbond(name1, my_file)

        name1='150'
        my_file='hbond_' // TRIM(name1) // '.dat'
        CALL get_hbond(name1, my_file)

        name1='140'
        my_file='hbond_' // TRIM(name1) // '.dat'
        CALL get_hbond(name1, my_file)

        name1='wernet'
        my_file='hbond_' // TRIM(name1) // '.dat'
        CALL get_hbond(name1, my_file)

        DEALLOCATE( cell )

        DEALLOCATE( coor )
        DEALLOCATE( symbol )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,2(1x,i10))
110     FORMAT(1x,a40,1(1x,e20.10))
310     FORMAT(1x,a40,3(1x,e20.10))

        END PROGRAM hbond
