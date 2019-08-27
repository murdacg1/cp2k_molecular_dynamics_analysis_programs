        PROGRAM com_drift

        USE common_data_mod
        USE reading_writing_mod

! Purpose: Correct the com drift problem. For each frame:
! (1) re-center the pos file (center = center of Si atoms in 0th frame)
! (2) remove com velocity from the vel file
! Program is run as:
!       com_drift.exe < com_drift.in >& com_drift.out &
!
! Program reads first the file com_drift.in:
! 6000                       f               total number of frames, not including the zeroth frame (total run=0.4fs*25k=10ps)
! 127                        a               number of atoms
! quartz-pos-1.xyz           in_pos_file
! quartz-pos-1_shifted.xyz   out_pos_file
! quartz-vel-1.xyz           in_vel_file
! quartz-vel-1_shifted.xyz   out_vel_file

        IMPLICIT NONE

        INTEGER :: iflag
        CHARACTER(len=80) in_pos_file, out_pos_file, in_vel_file, out_vel_file

        WRITE(6,'(a)') '# PROGRAM = com_drift'
        CALL FLUSH(6)

        READ(5,*) f
        READ(5,*) a
        READ(5,*) in_pos_file
        READ(5,*) out_pos_file
        READ(5,*) in_vel_file
        READ(5,*) out_vel_file

        WRITE(6,100) '# f = ', f
        WRITE(6,100) '# a = ', a
        WRITE(6,'(a)') '# in_pos_file = ', in_pos_file
        WRITE(6,'(a)') '# out_pos_file = ', out_pos_file
        WRITE(6,'(a)') '# in_vel_file = ', in_vel_file
        WRITE(6,'(a)') '# out_vel_file = ', out_vel_file
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        iflag = 1
        CALL reading_writing(iflag, in_pos_file, out_pos_file)

        iflag = 2
        CALL reading_writing(iflag, in_vel_file, out_vel_file)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1x,a40,1(1x,i10))

        END PROGRAM com_drift
