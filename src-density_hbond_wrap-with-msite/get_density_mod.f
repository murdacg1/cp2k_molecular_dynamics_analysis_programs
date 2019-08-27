        MODULE get_density_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE get_density(name1, my_file)

        CHARACTER(len=1) name1
        CHARACTER(len=80) my_file

        INTEGER :: f0, jf, jm, ja, jd, jb
        REAL(KIND=dp) :: z_temp, volume, dvolume, contrib, &
          factor, m_check

        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: &
          den, sd
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: &
          den_frame

        REAL(KIND=dp), PARAMETER :: Na = 6.02214179e23_dp, ang2cm = 1.0e-8_dp

        WRITE(6,'(2a)') '# SUBROUTINE = get_density; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        contrib = 1.0_dp
 
        ALLOCATE( den(0:b) )
        ALLOCATE( sd(0:b) )
        ALLOCATE( den_frame(0:f,0:b) )

        den = zero
        sd = zero
        den_frame = zero

        volume = 1.0_dp
        DO jd = 1, d
          volume = volume*cell(jd)*ang2cm
        END DO

        dvolume = volume * dz/cell(d)
        IF ( b /= 0 ) dvolume = 2.0_dp * dvolume ! because of the wrapping below
        factor = ( SUM(mass) ) / ( dvolume*Na )

        f0 = f - f_avg
        DO jf = f0+1, f

          DO jm = 1, m

            z_temp = coor(jf,jm,1,d)

            IF ( z_temp <= cell(d)/2.0_dp ) THEN ! use both surfaces
              z_temp = -z_temp + cell(d)/2.0_dp
            ELSE
              z_temp =  z_temp - cell(d)/2.0_dp
            END IF

            IF ( bulk_versus_surface == 0 ) THEN
              jb = INT(z_temp/dz) ! NINT
              IF ( ABS(dz - cell(d) ) < eps ) jb = 0 ! bulk only calc
            ELSE
              IF ( z_temp < (z_gds-2.0_dp*delta) ) THEN
                jb = 0 ! bulk region
              ELSE
                jb = 1 ! surface region
              END IF
            END IF

            IF  ( (jb >= 0) .AND. (jb <= b) ) THEN

              den_frame(jf,jb) = den_frame(jf,jb) + contrib

            ELSE

              WRITE(7,*) 'jf, jm, jb, b = ', jf, jm, jb, b
              WRITE(7,*) 'coor(jf,jm,1,d), z_temp = ', &
                coor(jf,jm,1,d), z_temp
              WRITE(7,*)
              CALL FLUSH(7)

            END IF

          END DO

        END DO

        m_check = zero
        DO jb = 0, b
          DO jf = f0+1, f
            den(jb) = den(jb) + den_frame(jf,jb)
          END DO
          m_check = m_check + den(jb)
          den(jb) = den(jb) * factor / REAL(f_avg,dp)
        END DO
        m_check = m_check / REAL(f_avg,dp)
        IF ( ABS(m_check - REAL(m,dp) ) > eps ) THEN
          WRITE(7,*) 'm, m_check = ', m, m_check
          FLUSH(7)
        END IF

        DO jb = 0, b
          DO jf = f0+1, f
            den_frame(jf,jb) = den_frame(jf,jb) * factor
            sd(jb) = sd(jb) + ( den(jb) - den_frame(jf,jb) )**2
          END DO
          sd(jb) = dsqrt( sd(jb) / REAL(f_avg,dp) )
        END DO

        WRITE(6,'(2a)') '# name1 = ', name1
        WRITE(6,100) '# m = ', m
        WRITE(6,110) '# m_check = ', m_check
        WRITE(6,100) '# a = ', a
        WRITE(6,110) '# volume [A^3] = ', volume / ang2cm**3
        WRITE(6,110) '# dvolume [A^3] = ', dvolume / ang2cm**3
        WRITE(6,110) '# contrib = ', contrib
        WRITE(6,110) '# factor = ', factor
        WRITE(6,300) '# f_avg, f, f0+1 = ', &
          f_avg, f, f0+1
        WRITE(6,100) '# b = ', b
        WRITE(6,110) '# dz = ', dz

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,320) z(jb), den(jb), sd(jb)
          CALL FLUSH(3)
        END DO
        CLOSE(3)

        DEALLOCATE( den )
        DEALLOCATE( sd )
        DEALLOCATE( den_frame )

100     FORMAT(1x,a40,1(1x,i10))
110     FORMAT(1x,a40,1(1x,f20.10))
200     FORMAT(1x,a40,2(1x,i10))
210     FORMAT(1x,a40,2(1x,f20.10))
300     FORMAT(1x,a40,3(1x,i10))
320     FORMAT(3(1x,e20.10))

        END SUBROUTINE get_density

        END MODULE get_density_mod
