        MODULE info_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE info(my_d1, my_d2, my_array1, my_array2, &
          my_file)

        INTEGER :: my_d1, my_d2
        REAL(KIND=dp), DIMENSION(0:f,m,a,my_d1) :: my_array1
        REAL(KIND=dp), DIMENSION(0:f,m,my_d2) :: my_array2
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, jd1, jd2
        REAL(KIND=dp) &
          vv, &
          X, Y, Z, &
          angle_X, angle_Y, angle_Z

        REAL(KIND=dp) :: pi, rad2deg

        WRITE(6,'(2a)') '# SUBROUTINE = info; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        pi = ACOS(-1.0_dp)
        rad2deg = 180.0_dp/pi

        my_array2 = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(2a)') &
          '#  1f  2tim  3m  ', &
          '4X  5Y  6Z  7RR  8angle_X  9angle_Y  10angle_Z'

        DO jf = 0, f
          DO jm = 1, m

            vv = zero
            DO jd1 = 1, my_d1
              my_array2(jf,jm,jd1) = my_array1(jf,jm,1,jd1)
              vv = vv + my_array2(jf,jm,jd1)**2
            END DO
            vv = SQRT(vv)
            my_array2(jf,jm,4) = vv
            X = my_array2(jf,jm,1)
            Y = my_array2(jf,jm,2)
            Z = my_array2(jf,jm,3)

            angle_X = rad2deg*ACOS(X/vv)
            angle_Y = rad2deg*ACOS(Y/vv)
            angle_Z = rad2deg*ACOS(Z/vv)

            my_array2(jf,jm,5) = angle_X
            my_array2(jf,jm,6) = angle_Y
            my_array2(jf,jm,7) = angle_Z

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (my_array2(jf,jm,jd2),jd2=1,my_d2)

          END DO
        END DO

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,e20.10),1(1x,i4),7(1x,e20.10))

        END SUBROUTINE info

        END MODULE info_mod
