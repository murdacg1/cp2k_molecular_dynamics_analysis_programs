        MODULE pos_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE pos(my_d1, my_d2, my_d3, &
          my_file)

        INTEGER :: my_d1, my_d2, my_d3
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, jd2, jd3, f_counter
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
        CALL FLUSH(6)

        pi = ACOS(-1.0_dp)
        rad2deg = 180.0_dp/pi

        coor_mol = zero
        coor_mon = zero

        OPEN(3,file=TRIM(my_file),status='unknown')

        WRITE(3,'(3a)') &
          '#  1f  2tim  3m  ', &
          '4X  5Y  6Z  7RR  8angle_X  9angle_Y  10angle_Z  ', &
          '11r1  12r2  13r  14theta'

        f_counter = 0

        r_avg = zero
        rsq_avg = zero
        r_rmsd = zero
        theta_avg = zero
        thetasq_avg = zero
        theta_rmsd = zero

        DO jf = f0+1, f

          f_counter = f_counter + 1

          DO jm = 1, m

            vv = zero
            DO jd2 = 1, 3
              coor_mol(jf,jm,jd2) = coor(jf,jm,1,jd2)
              vv = vv + coor_mol(jf,jm,jd2)**2
            END DO
            vv = SQRT(vv)
            coor_mol(jf,jm,4) = vv

            x1 = coor(jf,jm,2,1) - coor(jf,jm,1,1)
            y1 = coor(jf,jm,2,2) - coor(jf,jm,1,2)
            z1 = coor(jf,jm,2,3) - coor(jf,jm,1,3)
            x1 = x1 - cell(1)*ANINT(x1/cell(1))
            y1 = y1 - cell(2)*ANINT(y1/cell(2))
            z1 = z1 - cell(3)*ANINT(z1/cell(3))

            x2 = coor(jf,jm,3,1) - coor(jf,jm,1,1)
            y2 = coor(jf,jm,3,2) - coor(jf,jm,1,2)
            z2 = coor(jf,jm,3,3) - coor(jf,jm,1,3)
            x2 = x2 - cell(1)*ANINT(x2/cell(1))
            y2 = y2 - cell(2)*ANINT(y2/cell(2))
            z2 = z2 - cell(3)*ANINT(z2/cell(3))

            r1 = SQRT(x1**2 + y1**2 + z1**2)
            r2 = SQRT(x2**2 + y2**2 + z2**2)
            r = (r1+r2)/2.0d0
            theta = rad2deg * ACOS( (x1*x2+y1*y2+z1*z2) / &
              (r1*r2) )

            coor_mon(jf,jm,1) = r1
            coor_mon(jf,jm,2) = r2
            coor_mon(jf,jm,3) = r
            coor_mon(jf,jm,4) = theta

            bisector_X = (x1+x2)/2.0d0
            bisector_Y = (y1+y2)/2.0d0
            bisector_Z = (z1+z2)/2.0d0
            bisector = &
              SQRT(bisector_X**2 + bisector_Y**2 + bisector_Z**2)

            angle_X = rad2deg*ACOS(bisector_X/bisector)
            angle_Y = rad2deg*ACOS(bisector_Y/bisector)
            angle_Z = rad2deg*ACOS(bisector_Z/bisector)

            coor_mol(jf,jm,5) = angle_X
            coor_mol(jf,jm,6) = angle_Y
            coor_mol(jf,jm,7) = angle_Z

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (coor_mol(jf,jm,jd2),jd2=1,my_d2), &
              (coor_mon(jf,jm,jd3),jd3=1,my_d3)

            r_avg = r_avg + r
            rsq_avg = rsq_avg + r**2
            theta_avg = theta_avg + theta
            thetasq_avg = thetasq_avg + theta**2

          END DO
        END DO

        CLOSE(3)

        r_avg = r_avg/REAL(f_counter*m,dp)
        rsq_avg = rsq_avg/REAL(f_counter*m,dp)
        r_rmsd = SQRT(rsq_avg - r_avg**2)
        theta_avg = theta_avg/REAL(f_counter*m,dp)
        thetasq_avg = thetasq_avg/REAL(f_counter*m,dp)
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
