        PROGRAM avg_cubes

        USE common_data_mod
        USE reading_writing_mod

! Purpose: Read in cube file cube(x,y,z) (cube = v_hartree, ELECTRON_DENSITY, or TOTAL_DENSITY). Average over (x,y). Write out avg(z).
! Program is run as:
!       avg_cubes.exe < avg_cubes.in >& avg_cubes.out &
!
! Program reads first the file avg_cubes.in:
! v_hartree.cube   in_file
! v_hartree_z.dat  out_file

        IMPLICIT NONE

        INTEGER :: iflag
        CHARACTER(len=80) in_pos_file, out_pos_file, in_vel_file, out_vel_file

        WRITE(6,'(a)') '# PROGRAM = avg_cubes'
        CALL FLUSH(6)

        READ(5,*) in_file
        READ(5,*) out_file

        WRITE(6,'(a)') '# in_file = ', in_file
        WRITE(6,'(a)') '# out_file = ', out_file

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        CALL reading_writing

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        END PROGRAM avg_cubes
