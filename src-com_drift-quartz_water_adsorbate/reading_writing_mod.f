        MODULE reading_writing_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading_writing(iflag, in_file, out_file)

        INTEGER :: iflag
        CHARACTER(len=80) :: in_file, out_file

        INTEGER :: jf, ja, jd
        CHARACTER(len=80) :: title

        CHARACTER(len=2) :: symbol_temp
        CHARACTER(len=3), DIMENSION(:), ALLOCATABLE :: symbol
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: atom_array
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass

        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          middle_orig, middle, shift, middle_check

        INTEGER :: n_si
        INTEGER :: jf_temp
        INTEGER :: a_temp
        REAL(KIND=dp) :: mass_tot
        REAL(KIND=dp) :: shift_mag
        REAL(KIND=dp) :: tim, pot
        CHARACTER(len=3) char3_temp
        CHARACTER(len=8) char8_temp
        CHARACTER(len=5) char5_temp

        WRITE(6,*) &
         '# SUBROUTINE = reading_writing; ', &
         'd, f, a, iflag, in_file, out_file = ', &
          d, f, a, iflag, TRIM(in_file), '   ', TRIM(out_file)
        CALL FLUSH(6)

        ALLOCATE( symbol(a), atom_array(a,d), mass(a) )
        symbol = '  '
        atom_array = zero
        mass = zero

        ALLOCATE( middle_orig(d), middle(d), shift(d), middle_check(d) )

        SELECT CASE (iflag)

        CASE (1)

          OPEN(3,file=TRIM(in_file),status='unknown')
          DO jf = 0, 0
            READ(3,*) a_temp
            IF ( ABS(REAL(a-a_temp,dp)) > eps ) THEN
              WRITE(6,*) 'STOPPING: a, a_temp = ', a, a_temp
              STOP
            END IF
            READ(3,200) &
              char3_temp, jf_temp, &
              char8_temp, tim, &
              char5_temp, pot
!           WRITE(6,100) a
!           CALL FLUSH(6)
!           WRITE(6,200) &
!             char3_temp, jf_temp, &
!             char8_temp, tim, &
!             char5_temp, pot
!           CALL FLUSH(6)
            n_si = 0
            middle_orig = zero
            DO ja = 1, a
              READ(3,*) symbol_temp, atom_array(ja,:) ! , (atom_array(ja,jd),jd=1,d)
              symbol(ja) = ' ' // symbol_temp
!             WRITE(6,300) symbol(ja), atom_array(ja,:) ! , (atom_array(ja,jd),jd=1,d)
!             CALL FLUSH(6)
              IF ( symbol(ja) == ' Si') THEN
                n_si = n_si + 1
                middle_orig(:) = middle_orig(:) + atom_array(ja,:)
              END IF
            END DO
            middle_orig = middle_orig / REAL(n_si,dp) ! original center of Si atoms
            WRITE(6,*) 'a, jf, n_si, middle_orig = ', &
                        a, jf, n_si, middle_orig
            CALL FLUSH(6)
          END DO
          CLOSE(3)

          OPEN(3,file=TRIM(in_file),status='unknown')
          OPEN(4,file=TRIM(out_file),status='unknown')
          OPEN(7,file='shift_pos.xyz',status='unknown')
          DO jf = 0, f
            READ(3,*) a_temp
            READ(3,200) &
              char3_temp, jf_temp, &
              char8_temp, tim, &
              char5_temp, pot
            WRITE(4,100) a
            WRITE(4,200) &
              char3_temp, jf_temp, &
              char8_temp, tim, &
              char5_temp, pot
            middle = zero
            DO ja = 1, a
              READ(3,*) symbol_temp, atom_array(ja,:)
              symbol(ja) = ' ' // symbol_temp
              IF ( symbol(ja) == ' Si') THEN
                middle(:) = middle(:) + atom_array(ja,:)
              END IF
            END DO
            middle = middle / REAL(n_si,dp) ! current center of Si atoms
!           WRITE(6,*) 'a, jf, n_si, middle = ', &
!                       a, jf, n_si, middle
!           CALL FLUSH(6)
            shift = middle_orig - middle
            DO ja = 1, a
              WRITE(4,300) symbol(ja), atom_array(ja,:) + shift(:)
            END DO
            shift_mag = zero
            DO jd = 1, d
              shift_mag = shift_mag + shift(jd)**2
            END DO
            shift_mag = SQRT(shift_mag)
            WRITE(7,500) jf, shift(:), shift_mag
          END DO
          CLOSE(3)
          CLOSE(4)
          CLOSE(7)

        CASE (2)

          OPEN(3,file=TRIM(in_file),status='unknown')
          DO jf = 0, 0
            READ(3,*) a_temp
            IF ( ABS(REAL(a-a_temp,dp)) > eps ) THEN
              WRITE(6,*) 'STOPPING: a, a_temp = ', a, a_temp
              STOP
            END IF 
            READ(3,200) &
              char3_temp, jf_temp, &
              char8_temp, tim, &
              char5_temp, pot
            mass_tot = zero
            middle = zero
            DO ja = 1, a
              READ(3,*) symbol_temp, atom_array(ja,:)
              symbol(ja) = ' ' // symbol_temp
              SELECT CASE (symbol(ja))
              CASE (' Si')
                mass(ja) = 27.9769271_dp
              CASE ('  O')
                mass(ja) = 15.99491463_dp
              CASE ('  H')
                mass(ja) = 1.0078250321_dp
              CASE (' Cl')
                mass(ja) = 34.968852721_dp
              CASE ('  S')
                mass(ja) = 31.97207070_dp
              CASE ('  X')
                mass(ja) = 0.0000000000_dp
              END SELECT
              mass_tot = mass_tot + mass(ja)
              middle(:) = middle(:) +  mass(ja) * atom_array(ja,:)
            END DO
            middle = middle / mass_tot ! original velocity of COM
            WRITE(6,*) 'a, jf, mass_tot, middle = ', &
                        a, jf, mass_tot, middle
            CALL FLUSH(6)
          END DO
          CLOSE(3)

          OPEN(3,file=TRIM(in_file),status='unknown')
          OPEN(4,file=TRIM(out_file),status='unknown')
          OPEN(7,file='shift_vel.xyz',status='unknown')
          DO jf = 0, f
            READ(3,*) a_temp
            READ(3,200) &
              char3_temp, jf_temp, &
              char8_temp, tim, &
              char5_temp, pot
            WRITE(4,100) a
            WRITE(4,200) &
              char3_temp, jf_temp, &
              char8_temp, tim, &
              char5_temp, pot
            middle = zero
            DO ja = 1, a
              READ(3,*) symbol_temp, atom_array(ja,:)
              symbol(ja) = ' ' // symbol_temp
              middle(:) = middle(:) + mass(ja) * atom_array(ja,:)
            END DO
            middle = middle / mass_tot ! current velocity of COM
            DO ja = 1, a
              WRITE(4,300) symbol(ja), atom_array(ja,:) - middle(:)
            END DO
            shift_mag = zero
            DO jd = 1, d
              shift_mag = shift_mag + middle(jd)**2
            END DO
            shift_mag = SQRT(shift_mag)
            WRITE(7,500) jf, middle(:), shift_mag
            middle_check = zero
            DO ja = 1, a
              middle_check(:) = middle_check(:) + mass(ja) * ( atom_array(ja,:) - middle(:) )
            END DO
            middle_check = middle_check / mass_tot ! current velocity of COM
            WRITE(8,*) jf, middle_check(:)
          END DO
          CLOSE(3)
          CLOSE(4)
          CLOSE(7)

        END SELECT

        DEALLOCATE( symbol, atom_array, mass ) 
        DEALLOCATE( middle_orig, middle, shift, middle_check )

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(i8)
200     FORMAT(1x,a3,1x,i8,a8,1x,f12.3,a5,1x,f20.10)
300     FORMAT(a,1x,3(f20.10))
500     FORMAT(1(1x,i10),4(1x,f20.10))

        END SUBROUTINE reading_writing

        END MODULE reading_writing_mod
