        MODULE binning_z_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE binning_z(my_array, &
          my_file)

        REAL(KIND=dp), DIMENSION(0:f) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: v ! length of vector to bin
        INTEGER :: jf, jc, jv
        INTEGER :: c ! number of bins
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: vector
        INTEGER, DIMENSION(:), ALLOCATABLE :: counts ! needed distribution
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: ncounts ! normalized distribution
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: center ! the center of the distribution
        REAL(KIND=dp) :: low_orig, low ! low of bins
        REAL(KIND=dp) :: high_orig, high ! high of bins
        REAL(KIND=dp) :: mean ! mean of bins (abs)
        REAL(KIND=dp) :: sd ! sd of mean of bins (abs)
        REAL(KIND=dp) :: area

        WRITE(6,'(2a)') '# SUBROUTINE = binning_z; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        v = f_avg ! length of vector to bin
        ALLOCATE( vector(v) )
        vector = zero
        jv = 0
        DO jf = f0+1, f
          jv = jv + 1
          vector(jv) = my_array(jf)
        END DO

        low_orig = minval(vector)
        high_orig = maxval(vector)
        mean = sum(vector)/REAL(v,dp)

        sd = zero
        DO jv = 1, v
          sd = sd + (vector(jv)-mean)**2
        END DO
        sd = SQRT(sd/REAL(v,dp))

        low = REAL(ANINT(low_orig/width),dp)*width
        high = REAL(ANINT(high_orig/width),dp)*width

        c = ANINT((high-low)/width)
        ALLOCATE( counts(-1:c+1) )
        ALLOCATE( ncounts(-1:c+1) )
        ALLOCATE( center(-1:c+1) )
        counts = 0
        ncounts = zero
        center = zero

        DO jc = -1, c+1
          center(jc) = low + (REAL(jc,dp)+0.5_dp)*width
        END DO
 
        DO jv = 1, v
          jc = NINT( (vector(jv)-low) / width )
!         WRITE(6,400) 'jv, vector(jv)-low = ', jv, vector(jv)-low
!         WRITE(6,400) 'jc, center(jc)-low = ', jc, center(jc)-low
!         WRITE(6,*)
!         CALL FLUSH(6)
          IF ( (jc >= -1) .AND. (jc <= c+1) ) &
            counts(jc) = counts(jc) + 1
        END DO

        area = REAL(sum(counts),dp)*width
        DO jc = -1, c+1
          ncounts(jc) = REAL(counts(jc),dp)/area
        END DO

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(a)') '# jc  center  counts  normalized_counts'

        DO jc = -1, c+1
          IF ( counts(jc) > 0 ) &
            WRITE(3,300) jc, center(jc), counts(jc), ncounts(jc)
        END DO

        CLOSE(3)

        WRITE(6,110) '# samples = ', v
        WRITE(6,110) '# bins: c = ', c
        WRITE(6,200) '# bin width = ', width
        WRITE(6,200) '# z_s = ', z_s
        WRITE(6,220) '# low_orig, low = ', low_orig, low
        WRITE(6,220) '# high_orig, high = ',high_orig, high
        WRITE(6,200) '# mean = ', mean
        WRITE(6,200) '# st. dev. = ', sd
        WRITE(6,200) '# area = ', area
        WRITE(6,200) '# area/width = ', area/width
        WRITE(6,110) '# sum(counts) = ', sum(counts)
        WRITE(6,200) '# sum(ncounts)*width = ', sum(ncounts)*width

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        DEALLOCATE( vector )
        DEALLOCATE( counts )
        DEALLOCATE( ncounts )
        DEALLOCATE( center )

110     FORMAT(1x,a40,1(1x,i10))
120     FORMAT(1x,a40,2(1x,i10))
200     FORMAT(1x,a40,1(1x,e20.10))
220     FORMAT(1x,a40,2(1x,e20.10))
300     FORMAT(1x,i6,1x,e20.10,1x,i20,1x,e20.10)
400     FORMAT(1x,a40,1(1x,i10),1(1x,f10.3))

        END SUBROUTINE binning_z

        END MODULE binning_z_mod
