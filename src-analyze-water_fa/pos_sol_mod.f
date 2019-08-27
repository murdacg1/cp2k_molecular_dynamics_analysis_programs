        MODULE pos_sol_mod

        USE common_data_mod

        IMPLICIT NONE

        CONTAINS

        FUNCTION cross(a, b)

        REAL(KIND=dp), DIMENSION(3) :: cross
        REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

        END FUNCTION cross

        SUBROUTINE pos_sol(my_m, my_a, my_d1, my_d2, my_d3, my_array, &
          my_file)

        INTEGER :: my_m, my_a, my_d1, my_d2, my_d3
        REAL(KIND=dp), DIMENSION(0:f, my_m, my_a, my_d1) :: my_array
        CHARACTER(len=80) my_file

        INTEGER :: jf, jm, jd1, jd2, jd3
        REAL(KIND=dp), DIMENSION(my_d1) :: &
          roh_vec, roh_unit, &
          rc2o_vec, rc2o_unit, &
          rch_vec, rch_unit, &
          rco_vec, rco_unit, &
          rbisector_vec, rbisector_unit, &
          n1_vec, n1_unit, &
          n2_vec, n2_unit, &
          n_vec, n_unit
        REAL(KIND=dp) :: &
          roh, &
          rc2o, &
          rch, &
          rco, &
          ahoc, &
          aoc2o, &
          ahc2o, &
          thcoh, &
          Z, &
          aoh_Z, &
          ach_Z, &
          aco_Z, &
          an_Z, &
          rbisector, &
          ac2o_Z, &
          abisector_Z, &
          RR, &
          n1, n2, n
        REAL(KIND=dp) :: pi, rad2deg

        WRITE(6,'(2a)') '# SUBROUTINE = pos_sol; filename = ', &
          TRIM(my_file)
        WRITE(6,*) 'my_m, my_a, my_d1, my_d2, my_d3 = ', &
          my_m, my_a, my_d1, my_d2, my_d3

        pi = ACOS(-1.0_dp)
        rad2deg = 180.0_dp/pi

        coor_sol_mol = zero
        coor_sol_mon = zero

       OPEN(3,file=TRIM(my_file),status='unknown')

       WRITE(3,'(6a)') &
          '#          1f         2tim   3m', &
          '           4X           5Y           6Z          7RR', &
          '        8aoh_Z      9ach_Z      10aco_Z', &
          '       11an_Z      12ac2o_Z   13abisector_Z', &
          '    14roh       15rc2o      16rch        17rco', &
          '       18ahoc      19aoc2o     20ahc2o      21thcoh'

        DO jf = f0, f

          DO jm = 1, my_m

! 2 3 O217_H218 O-H
            roh_vec = my_array(jf,jm,3,:) - my_array(jf,jm,2,:)
            roh_vec(:) = roh_vec(:) - cell(:)*ANINT(roh_vec(:)/cell(:)) ! do not forget to appy PBCs
            roh = SQRT(DOT_PRODUCT(roh_vec, roh_vec))
            roh_unit = roh_vec/roh
! 1 4 C216_O219 C=O
            rc2o_vec = my_array(jf,jm,4,:) - my_array(jf,jm,1,:)
            rc2o_vec(:) = rc2o_vec(:) - cell(:)*ANINT(rc2o_vec(:)/cell(:))
            rc2o = SQRT(DOT_PRODUCT(rc2o_vec, rc2o_vec))
            rc2o_unit = rc2o_vec/rc2o
! 1 5 C216_H220 C-H
            rch_vec = my_array(jf,jm,5,:) - my_array(jf,jm,1,:)
            rch_vec(:) = rch_vec(:) - cell(:)*ANINT(rch_vec(:)/cell(:))
            rch = SQRT(DOT_PRODUCT(rch_vec, rch_vec))
            rch_unit = rch_vec/rch
! 1 2 C216_O217 C-O
            rco_vec = my_array(jf,jm,2,:) - my_array(jf,jm,1,:)
            rco_vec(:) = rco_vec(:) - cell(:)*ANINT(rco_vec(:)/cell(:))
            rco = SQRT(DOT_PRODUCT(rco_vec, rco_vec))
            rco_unit = rco_vec/rco
! 3 2 1 H218_O217_C216 H-O-C
            ahoc = ACOS(DOT_PRODUCT(roh_unit, -rco_unit))*rad2deg
! 2 1 4 O217_C216_O219 O-C=O
            aoc2o = ACOS(DOT_PRODUCT(rco_unit, rc2o_unit))*rad2deg
! 5 1 4 H220_C216_O219 H-C=O
            ahc2o = ACOS(DOT_PRODUCT(rch_unit, rc2o_unit))*rad2deg
! 5 1 2 3 H220_C216_O217_H218 H-C-O-H dihedral angle: dihedral angle is angles between two planes
            n1_vec = cross(rco_unit, rch_unit) ! normal for plane 1
            n1 = SQRT(DOT_PRODUCT(n1_vec, n1_vec))
            n1_unit = n1_vec/n1
!           WRITE(6,*) 'n1_vec = ', n1_vec
!           WRITE(6,*) 'n1 = ', n1
!           WRITE(6,*) 'n1_unit = ', n1_unit
            n2_vec = cross(roh_unit, -rco_unit) ! normal for plane 2
            n2 = SQRT(DOT_PRODUCT(n2_vec, n2_vec))
            n2_unit = n2_vec/n2
            thcoh = ACOS(DOT_PRODUCT(n1_unit, n2_unit))*rad2deg 

            coor_sol_mon(jf,jm,1) = roh
            coor_sol_mon(jf,jm,2) = rc2o
            coor_sol_mon(jf,jm,3) = rch
            coor_sol_mon(jf,jm,4) = rco
            coor_sol_mon(jf,jm,5) = ahoc
            coor_sol_mon(jf,jm,6) = aoc2o
            coor_sol_mon(jf,jm,7) = ahc2o
            coor_sol_mon(jf,jm,8) = thcoh

            Z = my_array(jf,jm,1,3)
            aoh_Z = ACOS(DOT_PRODUCT(roh_unit, Z_unit))*rad2deg
            ach_Z = ACOS(DOT_PRODUCT(rch_unit, Z_unit))*rad2deg
            aco_Z = ACOS(DOT_PRODUCT(rco_unit, Z_unit))*rad2deg
            n_vec = cross(rc2o_unit, rco_unit) ! normal for plane of molecule
            n = SQRT(DOT_PRODUCT(n_vec, n_vec))
            n_unit = n_vec/n
            an_Z = ACOS(DOT_PRODUCT(n_unit, Z_unit))*rad2deg
            ac2o_Z = ACOS(DOT_PRODUCT(rc2o_unit, Z_unit))*rad2deg
            rbisector_vec =  rco_vec + rc2o_vec
            rbisector = SQRT(DOT_PRODUCT(rbisector_vec, rbisector_vec))
            rbisector_unit = rbisector_vec/rbisector
            abisector_Z = ACOS(DOT_PRODUCT(rbisector_unit,Z_unit))*rad2deg

            coor_sol_mol(jf,jm,:) = my_array(jf,jm,1,:)
            RR = SQRT( SUM( coor_sol_mol(jf,jm,:)**2 ) )
            coor_sol_mol(jf,jm,4) = RR
            coor_sol_mol(jf,jm,5) = aoh_Z
            coor_sol_mol(jf,jm,6) = ach_Z
            coor_sol_mol(jf,jm,7) = aco_Z
            coor_sol_mol(jf,jm,8) = an_Z
            coor_sol_mol(jf,jm,9) = ac2o_Z
            coor_sol_mol(jf,jm,10) = abisector_Z

            WRITE(3,100) &
              jf, tim(jf), jm, &
              (coor_sol_mol(jf,jm,jd2),jd2=1,my_d2), &
              (coor_sol_mon(jf,jm,jd3),jd3=1,my_d3)

          END DO

        END DO

        CALL FLUSH(3)

        CLOSE(3)

        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

100     FORMAT(1(1x,i12),1(1x,f12.6),1(1x,i4),18(1x,f12.6))
200     FORMAT(1x,a40,2(1x,f12.6))

        END SUBROUTINE pos_sol

        END MODULE pos_sol_mod
