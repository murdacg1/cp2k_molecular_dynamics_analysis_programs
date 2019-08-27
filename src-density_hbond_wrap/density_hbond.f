        PROGRAM density_hbond

        USE common_data_mod
        USE reading_mod
        USE get_density_mod
        USE get_hbond_mod

! Purpose: density for CP2K simulation of (H20)n

! Program is run as:
!       ../density_hbond.exe < density_hbond.in >& density_hbond.out &
!
! It is assumed that the file run-01.xyz have the atoms organized as O,H,H; O,H,H; ...

! This slab program will work for three cases depending on the input:

! (1a) First run: unknown z_gds and delta (set to large numbers ~ 999.0): dz other than cell(d): b /= 0, b /= 1 
! (1b) First run: unknown z_gds and delta (set to large numbers ~ 999.0): dz = cell(d): single z bin: b = 0
!      (for testing or for runnning on the bulk)
! (2)  Second run: known z_gds and delta (from fit to first run): now segregate into bulk--jb = 0 and surface--jb = 1: b = 1

! Program reads first the file density_hbond.in:

! 20000         f_avg (the final f frames are used for averaging) (avg calculated from last 5ps of the run)
! 29531         f (total number of frames, not including the zeroth frame) (total run=0.5fs*20k=10ps)
! 216           m (number of molecules per frame or cell)
! 3             a (number of atoms per molecule)
! 15.0 15.0 71.44  cell ABC values [Angstrom]
! 0.10 OR ... OR 71.44    dz [Angstrom] (cell(d)=71.44 means sum everything in one z bin)
! 15.99491463 1.0078250321 1.0078250321  mass1 mass2 mass3
! 4 4           aa dd (max number acceptors and donors)
! 999.0 OR 14.61         z_gds (z of gibbs dividing surface, Ang) (from tanh fit in gnuplot)
! 999.0 OR 1.007         delta (width of the interfacial region, Ang) (from tanh fit in gnuplot)

        IMPLICIT NONE

        REAL(KIND=dp) :: mass1, mass2, mass3
        INTEGER :: jaa, jdd
        INTEGER :: jb
        CHARACTER(len=1) name1
        CHARACTER(len=10) name2
        CHARACTER(len=80) my_file

        WRITE(6,'(a)') '# PROGRAM = density_hbond'
        CALL FLUSH(6)

        ALLOCATE( cell(d) )

        READ(5,*) f_avg
        READ(5,*) f
        READ(5,*) m
        READ(5,*) a
        READ(5,*) cell(:)
        READ(5,*) dz
        READ(5,*) mass1, mass2, mass3
        READ(5,*) aa, dd
        READ(5,*) z_gds
        READ(5,*) delta

        mass = mass1 + mass2 + mass3

        IF ( (z_gds >= 900.0_dp) .OR. (delta >= 900.0_dp) ) THEN
          bulk_versus_surface = 0
          b = INT( (cell(d)/2.0_dp) / dz )
        ELSE
          bulk_versus_surface = 1
          b = 1
          dz = cell(d)/2.0_dp
        END IF

        ALLOCATE( z(0:b) )
        z = zero
        DO jb = 0, b
          z(jb) = (REAL(jb,dp)+0.5_dp)*dz
        END DO

        WRITE(6,200) '# f_avg, f = ', f_avg, f
        WRITE(6,100) '# m = ', m
        WRITE(6,100) '# a = ', a
        WRITE(6,310) '# cell = ', cell(:)
        WRITE(6,110) '# dz = ', dz
        WRITE(6,100) '# b = ', b
        WRITE(6,410) '# mass1, mass2, mass3, mass [amu] = ', &
          mass1, mass2, mass3, mass
        WRITE(6,200) '# aa, dd = ', aa, dd
        WRITE(6,110) 'z_gds = ', z_gds
        WRITE(6,110) 'delta = ', delta

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( coor(0:f,m,a,d) )
        ALLOCATE( symbol(m,a) )

        my_file = '../ANALYZE/run-01.xyz'
        my_file = TRIM(my_file)
        WRITE(6,'(2a)') 'my_file = ', my_file
        CALL FLUSH(6)

        CALL reading(my_file)

        name1='O'
        my_file='density.dat'
        CALL get_density(name1, my_file)

        name2='140'
        my_file='hbond_' // TRIM(name2) // '.dat'
        ALLOCATE( pcounts_140(0:aa,0:dd,0:b) )
        ALLOCATE( num_hbonds_22_140(0:b) )
        ALLOCATE( num_hbonds_aadd_140(0:b) )
        CALL get_hbond(name2, pcounts_140, &
          num_hbonds_22_140, num_hbonds_aadd_140, my_file)

!        GOTO 999

        name2='150'
        my_file='hbond_' // TRIM(name2) // '.dat'
        ALLOCATE( pcounts_150(0:aa,0:dd,0:b) )
        ALLOCATE( num_hbonds_22_150(0:b) )
        ALLOCATE( num_hbonds_aadd_150(0:b) )
         CALL get_hbond(name2, pcounts_150, &
          num_hbonds_22_150, num_hbonds_aadd_150, my_file)

        name2='160'
        my_file='hbond_' // TRIM(name2) // '.dat'
        ALLOCATE( pcounts_160(0:aa,0:dd,0:b) )
        ALLOCATE( num_hbonds_22_160(0:b) )
        ALLOCATE( num_hbonds_aadd_160(0:b) )
         CALL get_hbond(name2, pcounts_160, &
          num_hbonds_22_160, num_hbonds_aadd_160, my_file)

        name2='wernet'
        my_file='hbond_' // TRIM(name2) // '.dat'
        ALLOCATE( pcounts_wernet(0:aa,0:dd,0:b) )
        ALLOCATE( num_hbonds_22_wernet(0:b) )
        ALLOCATE( num_hbonds_aadd_wernet(0:b) )
         CALL get_hbond(name2, pcounts_wernet, &
          num_hbonds_22_wernet, num_hbonds_aadd_wernet, my_file)

        IF ( (bulk_versus_surface == 1) .OR. (b == 0) ) THEN ! b = 0 means bulk calculation

          WRITE(6,'(a)') &
          '         Hbond cutoffs:     160  150  140  Wernet; etc.'
          WRITE(6,*)
          WRITE(6,*)
 
          DO jb = 0, b

            WRITE(6,510) 'jb, z(jb) = ', jb, z(jb)
            SELECT CASE (jb)
            CASE (0)
              WRITE(6,'(a)') '         Bulk region'
            CASE (1)
              WRITE(6,'(a)') '         Surface region'
            END SELECT
            WRITE(6,*)

            WRITE(6,'(a55,a25)') &
            '     jdd =                   0                        1', &
            '                        2'
            WRITE(6,*)
            DO jaa = 0, 2
              WRITE(6,530) 'jaa = ', jaa, &
                (pcounts_160(jaa,jdd,jb), &
                 pcounts_150(jaa,jdd,jb), &
                 pcounts_140(jaa,jdd,jb), &
                 pcounts_wernet(jaa,jdd,jb),jdd=0,2)
              WRITE(6,*)
            END DO
            WRITE(6,*)
            WRITE(6,*)
            FLUSH(6)

          END DO

          WRITE(6,'(a)') &
          'num_hbonds_22:    160  150  140  Wernet'
          WRITE(6,'(a)') &
          'num_hbonds_aadd:  160  150  140  Wernet'
          WRITE(6,*)
          WRITE(6,*)
 
          DO jb = 0, b

            SELECT CASE (jb)
            CASE (0)
              WRITE(6,'(a)') '         Bulk region'
            CASE (1)
              WRITE(6,'(a)') '         Surface region'
            END SELECT
            WRITE(6,*)

            WRITE(6,600) &
              num_hbonds_22_160(jb), num_hbonds_22_150(jb), &
              num_hbonds_22_140(jb), num_hbonds_22_wernet(jb)
            WRITE(6,600) &
              num_hbonds_aadd_160(jb), num_hbonds_aadd_150(jb), &
              num_hbonds_aadd_140(jb), num_hbonds_aadd_wernet(jb)
            WRITE(6,*)
            WRITE(6,*)
            FLUSH(6)

          END DO
 
        END IF

!999     CONTINUE

        DEALLOCATE( z )

        DEALLOCATE( cell )

        DEALLOCATE( coor )
        DEALLOCATE( symbol )

        DEALLOCATE( pcounts_140 )
        DEALLOCATE( pcounts_150 )
        DEALLOCATE( pcounts_160 )
        DEALLOCATE( pcounts_wernet )

        DEALLOCATE( num_hbonds_22_140 )
        DEALLOCATE( num_hbonds_aadd_140 )
        DEALLOCATE( num_hbonds_22_150 )
        DEALLOCATE( num_hbonds_aadd_150 )
        DEALLOCATE( num_hbonds_22_160 )
        DEALLOCATE( num_hbonds_aadd_160 )
        DEALLOCATE( num_hbonds_22_wernet )
        DEALLOCATE( num_hbonds_aadd_wernet )

100     FORMAT(1x,a40,1(1x,i10))
200     FORMAT(1x,a40,2(1x,i10))
110     FORMAT(1x,a40,1(1x,e20.10))
310     FORMAT(1x,a40,3(1x,e20.10))
410     FORMAT(1x,a40,4(1x,e20.10))
510     FORMAT(1x,a20,1(1x,i10),3(1x,f10.4))
520     FORMAT(1x,a10,3(17x,i1))
530     FORMAT(1x,a10,1(1x,i10),3(5x,4(1x,f6.1)))
600     FORMAT(4(1x,f10.1))

        END PROGRAM density_hbond
