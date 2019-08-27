        PROGRAM bin_times

        USE common_data_mod
        USE binning_mod

! Purpose: Bin ionization times
! Program is run as:
!       bin_times.exe < bin_times.in >& bin_times.out &
!
! Program reads first the file bin_times.in:
! times.dat     data_file (second column is the data)
! 11            n        (number of entries)
! 101           c        (number of bins = c + 1)

        IMPLICIT NONE

        INTEGER :: n, jn, i
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: times
        CHARACTER(len=80) data_file, my_file

        WRITE(6,'(a)') '# PROGRAM = bin_times'

        READ(5,*) data_file
        READ(5,*) n
        READ(5,*) c

        WRITE(6,'(a)') '# data_file = ', TRIM(data_file)
        WRITE(6,100) '# n = ', n
        WRITE(6,100) '# number of bins c = ', c
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        ALLOCATE( times(n) )

        OPEN(3,file=TRIM(data_file),status='unknown')
        WRITE(6,'(a)') 'i, jn, times(jn):'
        DO jn = 1, n
          READ(3,*) i, times(jn)
          WRITE(6,*) i, jn, times(jn)
        END DO
        CLOSE(3)
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

! bin times
        my_file='times.bin'
        CALL binning(n, times, my_file)

        DEALLOCATE( times)

100     FORMAT(1x,a40,1(1x,i10))

        END PROGRAM bin_times
