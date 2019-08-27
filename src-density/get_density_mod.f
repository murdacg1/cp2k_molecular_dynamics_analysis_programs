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
          z, den, sd
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: &
          den_frame
 
        REAL(KIND=dp), PARAMETER :: Na=6.02214179e23_dp, ang2cm=1.0e-8_dp

        WRITE(6,'(2a)') '# SUBROUTINE = get_density; filename = ', &
          TRIM(my_file)
        CALL FLUSH(6)

        ALLOCATE( z(0:b) )
        ALLOCATE( den(0:b) )
        ALLOCATE( sd(0:b) )
        ALLOCATE( den_frame(0:f,0:b) )
 
        z = zero
        den = zero
        sd = zero
        den_frame = zero

        DO jb = 0, b
          z(jb) = (REAL(jb,dp)+0.5_dp)*dz
        END DO

        volume = 1.0_dp
        DO jd = 1, d
          volume = volume*cell(jd)*ang2cm
        END DO

        dvolume = volume * dz/cell(d)
        factor = (mass) / (dvolume*Na)
        contrib = 1.0_dp

        f0 = f - f_avg
        DO jf = f0+1, f
          DO jm = 1, m
            DO ja = 1, a
              IF ( symbol(jm,ja) /= name1 ) CYCLE
              z_temp = coor(jf,jm,ja,d)
              jb = NINT(z_temp/dz)
              IF (jb <= b) den_frame(jf,jb) = den_frame(jf,jb) + contrib
            END DO
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

        DO jb = 0, b
          DO jf = f0+1, f
            den_frame(jf,jb) = den_frame(jf,jb) * factor
            sd(jb) = sd(jb) + ( den(jb) - den_frame(jf,jb) )**2
          END DO
          sd(jb) = dsqrt( sd(jb) / REAL(f_avg,dp) )
        END DO

        OPEN(3,file=TRIM(my_file),status='unknown')
        DO jb = 0, b
          WRITE(3,320) z(jb), den(jb), sd(jb)
        END DO
        CLOSE(3)

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

        DEALLOCATE( z )
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
