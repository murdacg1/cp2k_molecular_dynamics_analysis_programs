        PROGRAM binning_program

        USE common_data_mod
        USE readit_mod
        USE binit_mod

! Purpose: This is a general program to read in a file of vectors in columns, focus on some columns: get avg, sd, and histogram
!
! Program is run as:
!       binning_program.exe < binning_program.in >& binning_program.out &
!
! Program reads first the file binning_program.in:
! coor_sol.dat  datafile
! 1             t          number of title lines
! 721211        m          total number of lines
! 21            n          total number of columns
! 6             cc         number of columns to bin; followed by jcc=1..cc:
! 6  Z_sol      cc_numbers cc_names(jcc)
! 14 roh_sol
! 11 an_Z_sol
! 10 aco_Z_sol
! 9  ach_Z_sol
! 12 ac2o_Z_sol
! 101           c          number of bins

        IMPLICIT NONE

        INTEGER :: my_t, my_m, my_n, my_d
        CHARACTER(len=80) datafile, my_file

        INTEGER :: jcc

        WRITE(6,'(a)') '# PROGRAM = binning_program'
        CALL FLUSH(6)

        READ(5,*) datafile
        READ(5,*) t
        READ(5,*) m
        READ(5,*) n
        READ(5,*) cc
        ALLOCATE( cc_numbers(cc) )
        ALLOCATE( cc_names(cc) )
        DO jcc = 1, cc
          READ(5,*) cc_numbers(jcc), cc_names(jcc)
        END DO
        READ(5,*) c

        WRITE(6,*) 'datafile = ', TRIM(datafile)
        WRITE(6,*) 't = ', t
        WRITE(6,*) 'm = ', m
        WRITE(6,*) 'n = ', n
        WRITE(6,*) 'cc = ', cc
        DO jcc = 1, cc
          WRITE(6,*) 'jcc, cc_numbers(jcc), cc_names(jcc) = ', jcc, cc_numbers(jcc), TRIM(cc_names(jcc))
        END DO
        WRITE(6,*) 'c = ', c

        ALLOCATE( title_lines(t) )
        ALLOCATE( full_array(m, n) )

        my_t = t
        my_m = m
        my_n = n
        my_file = TRIM(datafile)
        CALL readit(my_t, title_lines, my_m, my_n, full_array, datafile)
        DEALLOCATE( title_lines )

        DO jcc = 1, cc
          my_d = cc_numbers(jcc)
          my_file = TRIM(cc_names(jcc)) // '.bin'
          CALL binit(my_d, my_m, my_n, full_array, my_file)
        END DO

        DEALLOCATE( cc_numbers )
        DEALLOCATE( cc_names )
        DEALLOCATE( full_array )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        END PROGRAM binning_program
