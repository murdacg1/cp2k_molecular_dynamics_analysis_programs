        MODULE reading_writing_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE reading_writing

        INTEGER :: ja, jx, jy, jz

        WRITE(6,*) &
         '# SUBROUTINE = reading_writing; ', &
         'in_file, out_file = ', &
          TRIM(in_file), '   ', TRIM(out_file)
        CALL FLUSH(6)

        OPEN(1,file=TRIM(in_file),status='unknown')
        READ(1,*) title1
        READ(1,*) title2
        ALLOCATE( origin(d), h(d,d) )
        READ(1,*) a, origin(:)
        READ(1,*) nx, h(1,:)
        READ(1,*) ny, h(2,:)
        READ(1,*) nz, h(3,:)
        ALLOCATE( jatomicnum(a), atomicnumeff(a), coor(a,d) )
        DO ja = 1, a
          READ(1,*) &
            jatomicnum(ja), atomicnumeff(ja), coor(ja,:)
        END DO
        ALLOCATE( cube(nz,ny,nx) )
        DO jx = 1, nx
          DO jy = 1, ny
            READ(1,*) &
              (cube(jz,jy,jx), jz=1,nz)
          END DO
        END DO
        CLOSE(1)

        OPEN(2,file=TRIM(out_file),status='unknown')
        ALLOCATE( dat(nz) )
        ALLOCATE( z(nz) )
        dat = zero
        DO jz = 1, nz
          DO jy = 1, ny
            DO jx = 1, nx
              dat(jz) = dat(jz) + cube(jz,jy,jx)
            END DO
          END DO
          dat(jz) = dat(jz)/REAL(nx*ny,dp)
!         z(jz) = &
!           origin(3) + &
!           REAL(jx-1,dp)*h(1,3) + &
!           REAL(jy-1,dp)*h(2,3) + &
!           REAL(jz-1,dp)*h(3,3)
          z(jz) = origin(3) + REAL(jz-1,dp)*h(3,3)
          z(jz) = bohr2ang*z(jz)
          WRITE(2,200) z(jz), dat(jz)
        END DO
        CLOSE(2)

        DEALLOCATE( origin, h )
        DEALLOCATE( jatomicnum, atomicnumeff, coor )
        DEALLOCATE( cube )
        DEALLOCATE( dat )
        DEALLOCATE( z ) 

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

200     FORMAT(2(1x,e13.5))

        END SUBROUTINE reading_writing

        END MODULE reading_writing_mod
