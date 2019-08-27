        MODULE pos_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE pos(my_m, my_a, my_d1, my_d2, my_d3, my_array, &
          my_file)

        INTEGER :: my_m, my_a, my_d1, my_d2, my_d3
        REAL(KIND=dp), DIMENSION(0:f, my_m, my_a, my_d1) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, jd1, jd2, jd3, f_counter
        REAL(KIND=dp) :: &
          x1, y1, z1, &
          x2, y2, z2, &
          r1, r2, r, theta, &
          r_avg, rsq_avg, r_rmsd, &
          theta_avg, thetasq_avg, theta_rmsd, &
          vv, &
          angle_X, angle_Y, angle_Z, &
          bisector_X, bisector_Y, bisector_Z, bisector

        REAL(KIND=dp) :: pi, rad2deg

        WRITE(6,'(2a)') '# SUBROUTINE = pos; filename = ', &
          TRIM(my_file)
        WRITE(6,*) 'my_m, my_a, my_d1, my_d2, my_d3 = ', &
          my_m, my_a, my_d1, my_d2, my_d3
        CALL FLUSH(6)

        pi = ACOS(-1.0_dp)
        rad2deg = 180.0_dp/pi

        coor_wat_mol = zero
        coor_wat_mon = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(3a)') &
          '#  1f  2tim  3m  ', &
          '4X  5Y  6Z  7RR  8angle_X  9angle_Y  10angle_Z  ', &
          '11r1  12r2  13r  14theta'
        CALL FLUSH(3) 

        f_counter = 0

        r_avg = zero
        rsq_avg = zero
        r_rmsd = zero
        theta_avg = zero
        thetasq_avg = zero
        theta_rmsd = zero

        DO jf = 0, f

          f_counter = f_counter + 1

          DO jm = 1, my_m

            vv = zero
            DO jd1 = 1, my_d1
              coor_wat_mol(jf,jm,jd1) = my_array(jf,jm,1,jd1)
              vv = vv + coor_wat_mol(jf,jm,jd1)**2
            END DO
            vv = SQRT(vv)
            coor_wat_mol(jf,jm,4) = vv

            x1 = my_array(jf,jm,2,1) - my_array(jf,jm,1,1)
            y1 = my_array(jf,jm,2,2) - my_array(jf,jm,1,2)
            z1 = my_array(jf,jm,2,3) - my_array(jf,jm,1,3)
            x1 = x1 - cell(1)*ANINT(x1/cell(1))
            y1 = y1 - cell(2)*ANINT(y1/cell(2))
            z1 = z1 - cell(3)*ANINT(z1/cell(3))

            x2 = my_array(jf,jm,3,1) - my_array(jf,jm,1,1)
            y2 = my_array(jf,jm,3,2) - my_array(jf,jm,1,2)
            z2 = my_array(jf,jm,3,3) - my_array(jf,jm,1,3)
            x2 = x2 - cell(1)*ANINT(x2/cell(1))
            y2 = y2 - cell(2)*ANINT(y2/cell(2))
            z2 = z2 - cell(3)*ANINT(z2/cell(3))

            r1 = SQRT(x1**2 + y1**2 + z1**2)
            r2 = SQRT(x2**2 + y2**2 + z2**2)
            r = (r1+r2)/2.0_dp
            theta = rad2deg * ACOS( (x1*x2+y1*y2+z1*z2) / &
              (r1*r2) )

            coor_wat_mon(jf,jm,1) = r1
            coor_wat_mon(jf,jm,2) = r2
            coor_wat_mon(jf,jm,3) = r
            coor_wat_mon(jf,jm,4) = theta

            bisector_X = (x1+x2)/2.0_dp
            bisector_Y = (y1+y2)/2.0_dp
            bisector_Z = (z1+z2)/2.0_dp
            bisector = &
              SQRT(bisector_X**2 + bisector_Y**2 + bisector_Z**2)

            angle_X = rad2deg*ACOS(bisector_X/bisector)
            angle_Y = rad2deg*ACOS(bisector_Y/bisector)
            angle_Z = rad2deg*ACOS(bisector_Z/bisector)

            coor_wat_mol(jf,jm,5) = angle_X
            coor_wat_mol(jf,jm,6) = angle_Y
            coor_wat_mol(jf,jm,7) = angle_Z

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (coor_wat_mol(jf,jm,jd2),jd2=1,my_d2), &
              (coor_wat_mon(jf,jm,jd3),jd3=1,my_d3)
            CALL FLUSH(3) 

            r_avg = r_avg + r
            rsq_avg = rsq_avg + r**2
            theta_avg = theta_avg + theta
            thetasq_avg = thetasq_avg + theta**2

          END DO
        END DO

        CLOSE(3)

        r_avg = r_avg/REAL(f_counter*my_m,dp)
        rsq_avg = rsq_avg/REAL(f_counter*my_m,dp)
        r_rmsd = SQRT(rsq_avg - r_avg**2)
        theta_avg = theta_avg/REAL(f_counter*my_m,dp)
        thetasq_avg = thetasq_avg/REAL(f_counter*my_m,dp)
        theta_rmsd = SQRT(thetasq_avg - theta_avg**2)

        WRITE(6,200) '# r_avg  r_rmsd [Ang] = ', &
          r_avg, r_rmsd
        WRITE(6,200) '# theta_avg  theta_rmsd [deg] = ', &
          theta_avg, theta_rmsd

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i10),1(1x,e20.10),1(1x,i4),11(1x,e20.10))
200     FORMAT(1x,a40,2(1x,e20.10))

        END SUBROUTINE pos

        END MODULE pos_mod
