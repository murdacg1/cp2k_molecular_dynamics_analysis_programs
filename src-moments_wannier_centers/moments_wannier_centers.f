        program moments_wannier_centers

! Purpose: Calculate tesseral multipole moments as defined in
! AJ Stone's gdma program.  Input is a file containing positions of atoms
! (O, H charges = +6, +1) and positions of the the four MLWFCs per molecule
! on midbonds and lone pairs (charge = -2).  Thus any change in dipole and
! higher moments is due to shifts in the MLWFCs due to changes in electronic density
! (and a smaller changes due to small changes in postions of the nuclei).
! For now, position multipoles at O atom only.
! Thus, everything is as in CJM & Will Kuo Science paper
! on the BLYP interface, only differences:
! (a) now analyzing a BLYP-D simulation
! (b) examining also the quadrupole and higher multipoles
!
! Moment components Qlm with up to
! l=0,1,2,3,4 (charge,dipole,quadrupole,octupole,hexadecapole)
! are obtained.  The multipole moments Qlm are defined as intergrals
! over the regular solid harmonics Rlm:
!     Q_{lm} = \int d^3 {\bf r'} \rho({\bf r'}) R_{lm}(r-r')
! Ql, the magnitudes of the moments, are independent of the axis
! system, and thus also rotationally invariant:
!     |Q_l| = \sqrt ( \sum_m |Q_{lm}|^2 )
! Some moments in cartesian form are also printed,
! as well as principal values of the quadrupole obtained
! through a diagonalization
!
! Program is run as:
!       ./moments_wannier_centers.x < quartz-1.wfc >& quartz-wfc.mom &
!
! Program reads first the file moments.in:
! 1             nf (number of frames, avg and sd will be taken over these frames)
! 8             nm (number of molecules per frame)
! 3             na (number atoms per molecule)
! 4             nc (number Wannier centers per molecule)

        implicit none

        integer nf_max, nm_max, na_max, nc_max
        parameter (nf_max=74, nm_max=216,
     &             na_max=3, nc_max=4)
        integer nf, nm, na, nc
        character*80 title1,title2
        character*1 symbol
        real*8 coords(nf_max,nm_max,na_max+2*nc_max,3),
     &         charge(nf_max,nm_max,na_max+2*nc_max)
        real*8 temp(3), charge_temp
        real*8 atoms(nm_max,na_max,3), centers(nm_max*nc_max,3)
        integer f, m, a, c, i, j
        integer mc, m_closest
        integer centers_counter(nm_max)
        real*8 dist(nm_max)
        integer ccm
        real*8 dist_array(nm_max,2*nc_max)

        real*8 x,y,z

        real*8 pi
        real*8 bohr2ang
        parameter (bohr2ang=0.529177249d0)
        real*8 au2debye
        parameter (au2debye=1.0d0/0.39343029d0)
        real*8 const, rad2deg

        real*8 Q0(nf_max,nm_max),
     &         Q1(nf_max,nm_max),
     &         Q2(nf_max,nm_max),
     &         Q3(nf_max,nm_max),
     &         Q4(nf_max,nm_max)

        real*8 Q00(nf_max,nm_max),
     &         Q10(nf_max,nm_max),
     &         Q11c(nf_max,nm_max),Q11s(nf_max,nm_max),
     &         Q20(nf_max,nm_max),
     &         Q21c(nf_max,nm_max),Q21s(nf_max,nm_max),
     &         Q22c(nf_max,nm_max),Q22s(nf_max,nm_max),
     &         Q30(nf_max,nm_max),
     &         Q31c(nf_max,nm_max),Q31s(nf_max,nm_max),
     &         Q32c(nf_max,nm_max),Q32s(nf_max,nm_max),
     &         Q33c(nf_max,nm_max),Q33s(nf_max,nm_max),
     &         Q40(nf_max,nm_max),
     &         Q41c(nf_max,nm_max),Q41s(nf_max,nm_max),
     &         Q42c(nf_max,nm_max),Q42s(nf_max,nm_max),
     &         Q43c(nf_max,nm_max),Q43s(nf_max,nm_max),
     &         Q44c(nf_max,nm_max),Q44s(nf_max,nm_max)

        real*8 R10,R11c,R11s,
     &         R20,R21c,R21s,R22c,R22s,
     &         R30,R31c,R31s,R32c,R32s,R33c,R33s,
     &         R40,R41c,R41s,R42c,R42s,R43c,R43s,R44c,R44s

        real*8
     & theta(nf_max,nm_max,3,3),
     & theta_trace(nf_max,nm_max),
     & theta_det(nf_max,nm_max),
     & theta_standard(nf_max,nm_max,3)

        real*8 x1,y1,z1,x2,y2,z2
        real*8 r_OH1(nf_max,nm_max),
     &         r_OH2(nf_max,nm_max),
     &         theta_HOH(nf_max,nm_max)

        real*8 bisector_XX,bisector_YY,bisector_ZZ,bisector
        real*8 angle_XX(nf_max,nm_max),
     &         angle_YY(nf_max,nm_max),
     &         angle_ZZ(nf_max,nm_max)
        real*8 XX(nf_max,nm_max),
     &         YY(nf_max,nm_max),
     &         ZZ(nf_max,nm_max)
        real*8 ZZ_temp1, ZZ_temp2

        real *8 c3, c2, c1, c0, aa, bb, cc, CCC, BBB, AAA
        real*8 trace_old, det_old
        real*8 trace_new, det_new

        real*8 aL(3)
        data aL/15.00d0, 15.00d0, 71.44d0/
        real*8 dr(3)

        do i = 1, 3
          aL(i) = aL(i)/bohr2ang
        end do

        const=au2debye*bohr2ang
        pi=dacos(-1.0d0)
        rad2deg=180.0d0/pi

        open(1,file='moments.in',status='unknown')
        read(1,*) nf
        read(1,*) nm
        read(1,*) na
        read(1,*) nc
        close(1)

        write(6,*) 'nf = ', nf
        write(6,*) 'nm = ', nm
        write(6,*) 'na = ', na
        write(6,*) 'nc = ', nc
        call flush(6)

! read
! note: all coords in Angstroms, thus convert to bohr to get au for all quantities
! later convert to Angstrom, Debye, Debye*Angstrom, etc.

        do f = 1, nf ! frames
          do m = 1, nm ! molecules
            do a = 1, na+2*nc ! atoms+2*centers
               charge(f,m,a)=0.0d0
               do i = 1, 3
               coords(f,m,a,i)=0.0d0
              end do
            end do
          end do
        end do

        do f = 1, nf ! frames
          read(5,*) title1
          read(5,*) title2
          write(6,*) title1
          write(6,*) title2
          call flush(6)

          do m = 1, nm ! molecules
            do a = 1, na ! atoms
              read(5,*) symbol, (temp(i),i=1,3)
              if (symbol .eq. 'O') charge(f,m,a)=6.0d0
              if (symbol .eq. 'H') charge(f,m,a)=1.0d0
              do i = 1, 3
                atoms(m,a,i) = temp(i)/bohr2ang
!               atoms(m,a,i) =  atoms(m,a,i)-
!    &            aL(i)*anint(atoms(m,a,i)/aL(i))
                coords(f,m,a,i) = atoms(m,a,i)
              end do
            end do
          end do

          mc = 0
          do m = 1, nm ! molecules
            do c = 1, nc ! Wannier centers: note these are scrambled
              a = na + c
              mc = mc + 1
              read(5,*) symbol, (temp(i),i=1,3)
              if (symbol .eq. 'X') charge(f,m,a)=-2.0d0
              do i = 1, 3
                centers(mc,i) = temp(i)/bohr2ang
!               centers(mc,i) =  centers(mc,i) -
!    &            aL(i)*anint(centers(mc,i)/aL(i))
              end do
            end do
          end do

! for each Wannier center WC, assign it to the molecule that has the smallest Oxygen-WC distance

          do m = 1, nm ! keep track of where each WC has gone using a counter
            centers_counter(m) = 0
          end do

          do mc = 1, nm*nc ! loop over WCs

            do m = 1, nm ! loop over Oxygens
              dist(m)=0.0d0
              do i = 1, 3
                dr(i) = atoms(m,1,i) - centers(mc,i)
                dr(i) = dr(i) - aL(i)*anint(dr(i)/aL(i))
                dist(m) = dist(m) + dr(i)**2
              end do
              dist(m) = dsqrt(dist(m)) ! vector of distances to each Oxygen
            end do

            m_closest=1
            do m = 2, nm
              if (dist(m) .lt. dist(m_closest)) m_closest = m
            end do

            centers_counter(m_closest) = centers_counter(m_closest) + 1
            dist_array(m_closest,centers_counter(m_closest)) =
     &        dist(m_closest)
            a = na + centers_counter(m_closest) ! counter used to assign the WC to the proper atomic position in the molecule
            do i =1, 3
              coords(f,m_closest,a,i) = centers(mc,i)
            end do

          end do

! write out the molecules, check the centers_counter(m) to be sure it always equals nc (otherwise something got clobbered)
          do m = 1, nm ! molecules
            do a = 1, na ! atoms
              write(6,4001)
     &          f, m, a,
     &          charge(f,m,a), (coords(f,m,a,i),i=1,3),
     &          centers_counter(m)
            end do
            do ccm = 1, centers_counter(m) ! centers
              a = na + ccm
              write(6,4002)
     &          f, m, a,
     &          charge(f,m,a), (coords(f,m,a,i),i=1,3),
     &          ccm,
     &          dist_array(m,ccm)
              if (centers_counter(m) .ne. 4)
     &          write(20,*)
     &            'f, m, centers_counter(m) = ',
     &             f, m, centers_counter(m)
            end do
            write(6,*)
            call flush(6)
          end do

        end do

! calculate the multipole moments around the OXY atom

        do f = 1, nf

          do m = 1, nm

            Q00(f,m)  = 0.0d0

            Q10(f,m)  = 0.0d0
            Q11c(f,m) = 0.0d0
            Q11s(f,m) = 0.0d0

            Q20(f,m)  = 0.0d0
            Q21c(f,m) = 0.0d0
            Q21s(f,m) = 0.0d0
            Q22c(f,m) = 0.0d0
            Q22s(f,m) = 0.0d0

            Q30(f,m)  = 0.0d0
            Q31c(f,m) = 0.0d0
            Q31s(f,m) = 0.0d0
            Q32c(f,m) = 0.0d0
            Q32s(f,m) = 0.0d0
            Q33c(f,m) = 0.0d0
            Q33s(f,m) = 0.0d0

            Q40(f,m)  = 0.0d0
            Q41c(f,m) = 0.0d0
            Q41s(f,m) = 0.0d0
            Q42c(f,m) = 0.0d0
            Q42s(f,m) = 0.0d0
            Q43c(f,m) = 0.0d0
            Q43s(f,m) = 0.0d0
            Q44c(f,m) = 0.0d0
            Q44s(f,m) = 0.0d0

            do a = 1, na+nc

              charge_temp = charge(f,m,a)
              x = coords(f,m,a,1) - coords(f,m,1,1)
              y = coords(f,m,a,2) - coords(f,m,1,2)
              z = coords(f,m,a,3) - coords(f,m,1,3)
              x = x - aL(1)*anint(x/aL(1))
              y = y - aL(2)*anint(y/aL(2))
              z = z - aL(3)*anint(z/aL(3))

              Q00(f,m) =Q00(f,m) +charge_temp

              Q10(f,m) =Q10(f,m) +charge_temp*R10 (z,y,x)
              Q11c(f,m)=Q11c(f,m)+charge_temp*R11c(z,y,x)
              Q11s(f,m)=Q11s(f,m)+charge_temp*R11s(z,y,x)

              Q20(f,m) =Q20(f,m) +charge_temp*R20 (z,y,x)
              Q21c(f,m)=Q21c(f,m)+charge_temp*R21c(z,y,x)
              Q21s(f,m)=Q21s(f,m)+charge_temp*R21s(z,y,x)
              Q22c(f,m)=Q22c(f,m)+charge_temp*R22c(z,y,x)
              Q22s(f,m)=Q22s(f,m)+charge_temp*R22s(z,y,x)

              Q30(f,m) =Q30(f,m) +charge_temp*R30 (z,y,x)
              Q31c(f,m)=Q31c(f,m)+charge_temp*R31c(z,y,x)
              Q31s(f,m)=Q31s(f,m)+charge_temp*R31s(z,y,x)
              Q32c(f,m)=Q32c(f,m)+charge_temp*R32c(z,y,x)
              Q32s(f,m)=Q32s(f,m)+charge_temp*R32s(z,y,x)
              Q33c(f,m)=Q33c(f,m)+charge_temp*R33c(z,y,x)
              Q33s(f,m)=Q33s(f,m)+charge_temp*R33s(z,y,x)

              Q40(f,m) =Q40(f,m) +charge_temp*R40 (z,y,x)
              Q41c(f,m)=Q41c(f,m)+charge_temp*R41c(z,y,x)
              Q41s(f,m)=Q41s(f,m)+charge_temp*R41s(z,y,x)
              Q42c(f,m)=Q42c(f,m)+charge_temp*R42c(z,y,x)
              Q42s(f,m)=Q42s(f,m)+charge_temp*R42s(z,y,x)
              Q43c(f,m)=Q43c(f,m)+charge_temp*R43c(z,y,x)
              Q43s(f,m)=Q43s(f,m)+charge_temp*R43s(z,y,x)
              Q44c(f,m)=Q44c(f,m)+charge_temp*R44c(z,y,x)
              Q44s(f,m)=Q44s(f,m)+charge_temp*R44s(z,y,x)

            end do

            Q0(f,m) = dsqrt(Q00(f,m)**2)
            Q1(f,m) = dsqrt(Q10(f,m)**2 +
     &                      Q11c(f,m)**2+Q11s(f,m)**2)
            Q2(f,m) = dsqrt(Q20(f,m)**2 +
     &                      Q21c(f,m)**2+Q21s(f,m)**2 +
     &                      Q22c(f,m)**2+Q22s(f,m)**2)
            Q3(f,m) = dsqrt(Q30(f,m)**2 +
     &                      Q31c(f,m)**2+Q31s(f,m)**2 +
     &                      Q32c(f,m)**2+Q32s(f,m)**2 +
     &                      Q33c(f,m)**2+Q33s(f,m)**2)
            Q4(f,m) = dsqrt(Q40(f,m)**2 +
     &                      Q41c(f,m)**2+Q41s(f,m)**2 +
     &                      Q42c(f,m)**2+Q42s(f,m)**2 +
     &                      Q43c(f,m)**2+Q43s(f,m)**2 +
     &                      Q44c(f,m)**2+Q44s(f,m)**2)

            theta(f,m,1,1) = -(1.0d0/2.0d0)*Q20(f,m) +
     &                        (1.0d0/2.0d0)*dsqrt(3.0d0)*Q22c(f,m)
            theta(f,m,1,2) =  (1.0d0/2.0d0)*dsqrt(3.0d0)*Q22s(f,m)
            theta(f,m,1,3) =  (1.0d0/2.0d0)*dsqrt(3.0d0)*Q21c(f,m)
            theta(f,m,2,2) = -(1.0d0/2.0d0)*Q20(f,m) -
     &                        (1.0d0/2.0d0)*dsqrt(3.0d0)*Q22c(f,m)
            theta(f,m,2,3) =  (1.0d0/2.0d0)*dsqrt(3.0d0)*Q21s(f,m)
            theta(f,m,3,3) =  Q20(f,m)
            theta(f,m,2,1) = theta(f,m,1,2)
            theta(f,m,3,1) = theta(f,m,1,3)
            theta(f,m,3,2) = theta(f,m,2,3)

            theta_trace(f,m) =
     & theta(f,m,1,1) +
     & theta(f,m,2,2) +
     & theta(f,m,3,3)

            theta_det(f,m) =
     & theta(f,m,1,1)*theta(f,m,2,2)*theta(f,m,3,3) +
     & 2.0d0*theta(f,m,1,2)*theta(f,m,2,3)*theta(f,m,1,3) -
     & theta(f,m,1,1)*theta(f,m,2,3)**2 -
     & theta(f,m,2,2)*theta(f,m,1,3)**2 -
     & theta(f,m,3,3)*theta(f,m,1,2)**2

! geometry of monomer

            x1 = coords(f,m,2,1) - coords(f,m,1,1)
            y1 = coords(f,m,2,2) - coords(f,m,1,2)
            z1 = coords(f,m,2,3) - coords(f,m,1,3)
            x1 = x1 - aL(1)*anint(x1/aL(1))
            y1 = y1 - aL(2)*anint(y1/aL(2))
            z1 = z1 - aL(3)*anint(z1/aL(3))

            x2 = coords(f,m,3,1) - coords(f,m,1,1)
            y2 = coords(f,m,3,2) - coords(f,m,1,2)
            z2 = coords(f,m,3,3) - coords(f,m,1,3)
            x2 = x2 - aL(1)*anint(x2/aL(1))
            y2 = y2 - aL(2)*anint(y2/aL(2))
            z2 = z2 - aL(3)*anint(z2/aL(3))

            r_OH1(f,m) = dsqrt(x1**2 + y1**2 + z1**2)
            r_OH2(f,m) = dsqrt(x2**2 + y2**2 + z2**2)
            theta_HOH(f,m) = dacos( (x1*x2+y1*y2+z1*z2) /
     & (r_OH1(f,m)*r_OH2(f,m)) )
            r_OH1(f,m) = r_OH1(f,m)*bohr2ang
            r_OH2(f,m) = r_OH2(f,m)*bohr2ang
            theta_HOH(f,m) = theta_HOH(f,m)*rad2deg

! orientation of monomer

            bisector_XX = (x1+x2)/2.0d0
            bisector_YY = (y1+y2)/2.0d0
            bisector_ZZ = (z1+z2)/2.0d0
            bisector =
     & dsqrt(bisector_XX**2 + bisector_YY**2 + bisector_ZZ**2)
            angle_XX(f,m) = dacos(bisector_XX/bisector)
            angle_YY(f,m) = dacos(bisector_YY/bisector)
            angle_ZZ(f,m) = dacos(bisector_ZZ/bisector)
            XX(f,m) = coords(f,m,1,1)*bohr2ang
            YY(f,m) = coords(f,m,1,2)*bohr2ang
            ZZ_temp1 = coords(f,m,1,3) ! ZZ coordinate seems scrambled in file; thus fix it
            if (ZZ_temp1 .lt. 0.0d0) ZZ_temp2 = ZZ_temp1 + aL(3)/2.0d0
            if (ZZ_temp1 .gt. 0.0d0) ZZ_temp2 = ZZ_temp1 - aL(3)/2.0d0
            ZZ(f,m) = ZZ_temp2*bohr2ang
            angle_XX(f,m) = angle_XX(f,m)*rad2deg
            angle_YY(f,m) = angle_YY(f,m)*rad2deg
            angle_ZZ(f,m) = angle_ZZ(f,m)*rad2deg

            write(6,*)
            write(6,*)
            write(6,*) 'frame = ', f, 'molecule = ', m
            write(6,*) 'moments calculated from wannier centers file'
            write(6,*) '(dipole also in Debye)'
            write(6,*) '(1 bohr = 0.529177249 Angstrom)'

            write(6,*) 'Ang, deg:'
            write(6,350) 'r_OH1, r_OH2, theta_HOH = ',
     &                    r_OH1(f,m), r_OH2(f,m), theta_HOH(f,m)
            write(6,350) 'XX, YY, ZZ = ',
     &                    XX(f,m), YY(f,m), ZZ(f,m)
            write(6,350) 'angle_XX, angle_YY, angle_ZZ = ',
     &                    angle_XX(f,m), angle_YY(f,m), angle_ZZ(f,m)

            write(6,*)'In au:'
            write(6,300)
     & ' |Q0|  = ', Q0(f,m),
     & '  Q00  = ', Q00(f,m)

            write(6,310)
     & ' |Q1|  = ',        Q1(f,m),
     & '  Q11c = mu_x = ', Q11c(f,m),
     & '  Q11s = mu_y = ', Q11s(f,m),
     & '  Q10  = mu_z = ', Q10(f,m)

            write(6,320)
     & ' |Q2|  = ', Q2(f,m),
     & '  Q20  = ', Q20(f,m),
     & '  Q21c = ', Q21c(f,m),
     & '  Q21s = ', Q21s(f,m),
     & '  Q22c = ', Q22c(f,m),
     & '  Q22s = ', Q22s(f,m)

            write(6,330)
     & ' |Q3|  = ', Q3(f,m),
     & '  Q30  = ', Q30(f,m),
     & '  Q31c = ', Q31c(f,m),
     & '  Q31s = ', Q31s(f,m),
     & '  Q32c = ', Q32c(f,m),
     & '  Q32s = ', Q32s(f,m),
     & '  Q33c = ', Q33c(f,m),
     & '  Q33s = ', Q33s(f,m)

            write(6,340)
     & ' |Q4|  = ', Q4(f,m),
     & '  Q40  = ', Q40(f,m),
     & '  Q41c = ', Q41c(f,m),
     & '  Q41s = ', Q41s(f,m),
     & '  Q42c = ', Q42c(f,m),
     & '  Q42s = ', Q42s(f,m),
     & '  Q43c = ', Q43c(f,m),
     & '  Q43s = ', Q43s(f,m),
     & '  Q44c = ', Q44c(f,m),
     & '  Q44s = ', Q44s(f,m)

            write(6,*)'In Debye:'
            write(6,310)
     & ' |Q1|  = ',        Q1(f,m)*au2debye,
     & '  Q11c = mu_x = ', Q11c(f,m)*au2debye,
     & '  Q11s = mu_y = ', Q11s(f,m)*au2debye,
     & '  Q10  = mu_z = ', Q10(f,m)*au2debye

            write(6,*)'In au, quadrupole theta(i,j) in cartesian form:'
            write(6,*)   '    x         y         z'
            write(6,350) 'x', (theta(f,m,1,j), j=1,3)
            write(6,350) 'y', (theta(f,m,2,j), j=1,3)
            write(6,350) 'z', (theta(f,m,3,j), j=1,3)

            write(6,*)'In D*A, quadrupole theta(i,j) in cartesian form:'
            write(6,*)   '    x         y         z'
            write(6,350) 'x', (theta(f,m,1,j)*const, j=1,3)
            write(6,350) 'y', (theta(f,m,2,j)*const, j=1,3)
            write(6,350) 'z', (theta(f,m,3,j)*const, j=1,3)


! calculate principal values of the cartesian traceless quadrupole matrix theta(i,j)
! using the Viete method/subroutine (from my MS thesis work)
! c3*lmabda**3 + c2*lmabda**2 + c1*lmabda + c0 = 0
! or
! lmabda**3 + aa*lmabda**2 + bb*lmabda + cc = 0

            c3 = -1.0d0
            c2 = theta_trace(f,m)
            c1 =
     & - theta(f,m,1,1)*theta(f,m,2,2)
     & - theta(f,m,2,2)*theta(f,m,3,3)
     & - theta(f,m,1,1)*theta(f,m,3,3)
     & + theta(f,m,1,2)**2
     & + theta(f,m,2,3)**2
     & + theta(f,m,1,3)**2
            c0 = theta_det(f,m)
            aa = c2/c3
            bb = c1/c3
            cc = c0/c3
            trace_old = aa*c3
            det_old = c0*c3
            write(6,3002) 'DA, (DA)^3: trace_old, det_old = ',
     &                     trace_old*const, det_old*const**3
            call cubicn(aa,bb,cc,CCC,BBB,AAA)
            if (CCC .eq. 999.0d0) write(21,*) 'Roots complex!'
            write(6,*) 'AAA, CCC, BBB = ', AAA, CCC, BBB
            trace_new = AAA+CCC+BBB
            det_new = AAA*CCC*BBB
            write(6,3002) 'DA, (DA)^3: trace_new, det_new = ',
     &                     trace_new*const, det_new*const**3
            theta_standard(f,m,1)=AAA ! this ordering agrees with Laaksonen standard orientation moments: Qxx, Qyy, Qzz
            theta_standard(f,m,2)=CCC
            theta_standard(f,m,3)=BBB

! write it out in the most convenient units

            write(9,900) f, m,
     & XX(f,m), YY(f,m), ZZ(f,m), ! A=Angstroms
     & angle_XX(f,m), angle_YY(f,m), angle_ZZ(f,m), ! deg=degrees
     & r_OH1(f,m), r_OH2(f,m), theta_HOH(f,m), ! A, deg
     & Q00(f,m), ! au
     & Q1(f,m)*au2debye, ! D=Debye
     & Q11c(f,m)*au2debye, Q11s(f,m)*au2debye, Q10(f,m)*au2debye,
     & Q2(f,m)*const, theta_trace(f,m)*const, !DA
     & theta(f,m,1,1)*const, theta(f,m,1,2)*const, theta(f,m,1,3)*const,
     & theta(f,m,2,2)*const, theta(f,m,2,3)*const, theta(f,m,3,3)*const,
     & (theta_standard(f,m,i)*const,i=1,3)

          end do

        end do

! write it out in the most convenient units

         open(2,file='ref4-wc.dip_quad',status='unknown')

        do f = 1, nf ! frames
          write(2,*) 'frame = ', f, 'molecules = ', nm

          write(2,*) 'Dipole moments [D]:',
     &' f, m,',
     &' |p|, px, py, pz'
          do m = 1, nm ! molecules
            write(2,5001)
     & f, m,
     & Q1(f,m)*au2debye,
     & Q11c(f,m)*au2debye, Q11s(f,m)*au2debye, Q10(f,m)*au2debye
          end do

          write(2,*) 'Quadrupole moments [DA]:',
     &' f, m,',
     &' Qxx, Qxy, Qxz, Qyy, Qyz, Qzz;',
     &' standard orientation: Qxx, Qyy, Qzz'
          do m = 1, nm ! molecules
            write(2,5002)
     & f, m,
     & theta(f,m,1,1)*const, theta(f,m,1,2)*const, theta(f,m,1,3)*const,
     & theta(f,m,2,2)*const, theta(f,m,2,3)*const, theta(f,m,3,3)*const,
     & (theta_standard(f,m,i)*const,i=1,3)

          end do

        end do

100     format(1(1x,i4),3(1x,f11.6))
200     format(1(1x,i4),4(1x,f11.6))
250     format(6e13.5)
300     format(2(2x,a,e12.6))
310     format(4(2x,a,e12.6))
320     format(6(2x,a,e12.6))
330     format(8(2x,a,e12.6))
340     format(10(2x,a,e12.6))
350     format(2x,a,3f10.4)
900     format(2(1x,i4),25(1x,f10.4))

2001    format(1(4x,a,1x,f10.4))
2002    format(2(4x,a,1x,f10.4))
2003    format(3(4x,a,1x,f10.4))
2004    format(4(4x,a,1x,f10.4))

3001    format(1(1x,i4),2(1x,f10.4))
3002    format(1(1x,a),2(1x,f10.4))

4001    format(3(1x,i4),1(1x,f4.1),3(1x,f10.4),1(4x,i4))
4002    format(3(1x,i4),1(1x,f4.1),3(1x,f10.4),1(1x,i4),1(1x,f10.4))

5001    format(2(1x,i4),4(1x,f12.6))
5002    format(2(1x,i4),9(1x,f12.6))

        end


! Stone Rl0, Rlmc, Rlms operators


        real*8 function R10(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R10  = z

        end


        real*8 function R11c(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R11c = x

        end


        real*8 function R11s(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R11s = y

        end


        real*8 function R20(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z

        r2=x2+y2+z2

        R20  = (1.0d0/2.0d0) * (3.0d0*z2 - r2)

        end


        real*8 function R21c(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R21c = dsqrt(3.0d0) * x*z

        end


        real*8 function R21s(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R21s = dsqrt(3.0d0) * y*z

        end


        real*8 function R22c(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z

        x2=x*x

        y2=y*y

        R22c = dsqrt(3.0d0/4.0d0) * (x2 - y2)

        end


        real*8 function R22s(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R22s = dsqrt(3.0d0) * x*y

        end


        real*8 function R30(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,z3,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z
        z3=z2*z

        r2=x2+y2+z2

        R30  = (1.0d0/2.0d0) * (5.0d0*z3 - 3.0d0*z*r2)

        end


        real*8 function R31c(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z

        r2=x2+y2+z2

        R31c = dsqrt(3.0d0/8.0d0) * x * (5.0d0*z2 - r2)

        end


        real*8 function R31s(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z

        r2=x2+y2+z2

        R31s = dsqrt(3.0d0/8.0d0) * y * (5.0d0*z2 - r2)

        end


        real*8 function R32c(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z

        x2=x*x

        y2=y*y

        R32c = dsqrt(15.0d0/4.0d0) * z * (x2 - y2)

        end


        real*8 function R32s(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        R32s = dsqrt(15.0d0) * x*y*z

        end


        real*8 function R33c(z,y,x)
        implicit none
        real*8 x,x2,x3,
     &         y,y2,
     &         z


        x2=x*x
        x3=x2*x

        y2=y*y

        R33c = dsqrt(5.0d0/8.0d0) * (x3 - 3.0d0*x*y2)

        end


        real*8 function R33s(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,y3,
     &         z

        x2=x*x

        y2=y*y
        y3=y2*y

        R33s = dsqrt(5.0d0/8.0d0) * (3.0d0*x2*y-y3)

        end


        real*8 function R40(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,z3,z4,
     &         r2,r4

        x2=x*x

        y2=y*y

        z2=z*z
        z3=z2*z
        z4=z3*z

        r2=x2+y2+z2
        r4=r2*r2

        R40
     &          = (1.0d0/8.0d0) * (35.0d0*z4 - 30.0d0*z2*r2 + 3.0d0*r4)

        end


        real*8 function R41c(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,z3,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z
        z3=z2*z

        r2=x2+y2+z2

        R41c = (5.0d0/8.0d0) * (7.0d0*x*z3 - 3.0d0*x*z*r2)

        end


        real*8 function R41s(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,z3,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z
        z3=z2*z

        r2=x2+y2+z2

        R41s = (5.0d0/8.0d0) * (7.0d0*y*z3 - 3.0d0*y*z*r2)

        end


        real*8 function R42c(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z

        r2=x2+y2+z2

        R42c = (5.0d0/16.0d0) * (x2 - y2) * (7.0d0*z2 - r2)

        end


        real*8 function R42s(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,
     &         z,z2,
     &         r2

        x2=x*x

        y2=y*y

        z2=z*z

        r2=x2+y2+z2

        R42s = (5.0d0/4.0d0) * x*y * (7.0d0*z2 - r2)

        end


        real*8 function R43c(z,y,x)
        implicit none
        real*8 x,x2,x3,
     &         y,y2,
     &         z

        x2=x*x
        x3=x2*x

        y2=y*y

        R43c = (35.0d0/8.0d0) * z * (x3 - 3.0d0*x*y2)

        end


        real*8 function R43s(z,y,x)
        implicit none
        real*8 x,x2,
     &         y,y2,y3,
     &         z

        x2=x*x

        y2=y*y
        y3=y2*y

        R43s = (35.0d0/8.0d0) * z * (3.0d0*x2*y - y3)

        end


        real*8 function R44c(z,y,x)
        implicit none
        real*8 x,x2,x3,x4,
     &         y,y2,y3,y4,
     &         z

        x2=x*x
        x3=x2*x
        x4=x3*x

        y2=y*y
        y3=y2*y
        y4=y3*y

        R44c = (35.0d0/64.0d0) * (x4 - 6.0d0*x2*y2 + y4)

        end


        real*8 function R44s(z,y,x)
        implicit none
        real*8 x,x2,x3,
     &         y,y2,y3,
     &         z

        x2=x*x
        x3=x2*x

        y2=y*y
        y3=y2*y

        R44s = (35.0d0/4.0d0) * (x3*y - x*y3)

        end


! second moments to compare with cp2k printouts


        real*8 function Rxx(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxx  = x*x

        end


        real*8 function Rxy(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxy  = x*y

        end


        real*8 function Rxz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxz  = x*z

        end


        real*8 function Ryy(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Ryy  = y*y

        end


        real*8 function Ryz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Ryz  = y*z

        end


        real*8 function Rzz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rzz  = z*z

        end


! third moments to compare with cp2k printouts


        real*8 function Rxxx(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxxx  = x*x*x

        end


        real*8 function Rxxy(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxxy  = x*x*y

        end


        real*8 function Rxxz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxxz  = x*x*z

        end


        real*8 function Rxyy(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxyy  = x*y*y

        end


        real*8 function Rxyz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxyz  = x*y*z

        end


        real*8 function Rxzz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rxzz  = x*z*z

        end


        real*8 function Ryyy(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Ryyy  = y*y*y

        end


        real*8 function Ryyz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Ryyz  = y*y*z

        end


        real*8 function Ryzz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Ryzz  = y*z*z

        end


        real*8 function Rzzz(z,y,x)
        implicit none
        real*8 x,
     &         y,
     &         z

        Rzzz  = z*z*z

        end


! Francois Viete's method used in my ms thesis work

        subroutine cubicn(a,b,c,CCC,BBB,AAA)
        implicit none
        integer nn,k
        double precision pi,epsilon
        double precision a,b,c
        double precision Q,Q3,R,R2,theta,yr_1,yr_2,yr_3,arr(3),
     &                   CCC,BBB,AAA
        double precision AA,AA2,BB,sgnR,absR,absAA
        double complex i,yc_1,yc_2,yc_3
        data i / (0.0d0,1.0d0) /

        pi=dacos(-1.0d0)
        nn=3
        Q=(a*a-3.0d0*b)/9.0d0
        R=(2.0d0*a*a*a-9.0d0*a*b+27.0d0*c)/54.0d0
        R2=R*R
        Q3=Q*Q*Q

        if (R2 .lt. Q3) then

          print *, 'All three roots real.'
          theta=dacos(R/dsqrt(Q3))
          yr_1=-2.0d0*(dsqrt(Q))*dcos(theta/3.0d0)-a/3.0d0
          yr_2=-2.0d0*(dsqrt(Q))*dcos((theta+2.0d0*pi)/3.0d0)-a/3.0d0
          yr_3=-2.0d0*(dsqrt(Q))*dcos((theta-2.0d0*pi)/3.0d0)-a/3.0d0

          arr(1)=yr_1
          arr(2)=yr_2
          arr(3)=yr_3
!         print *, 'cubicn: beta2 = ',beta2
!         print *, 'unsorted arr = ', (arr(k), k=1,nn)
          call piksrt(nn,arr)
!         print *, 'sorted arr = ', (arr(k), k=1,nn)
          CCC=arr(1)
          BBB=arr(2)
          AAA=arr(3)
!         print *, 'CCC, BBB, AAA = ', CCC, BBB, AAA

        else

          print *, 'Not all three roots real.'
          sgnR=1.0d0
          absR=dsqrt(R2)
          if (R .lt. 0.0d0) sgnR=-1.0d0
          AA=-sgnR*( absR+dsqrt(R2-Q3) )**(1.0d0/3.0d0)
          AA2=AA*AA
          B=0.0d0
          absAA=dsqrt(AA2)
          if (absAA .gt. 0.d0) BB=Q/AA
          yc_1=(AA+BB)-a/3.0d0
          yc_2=-0.50d0*(AA+BB)-a/3.0d0+i*dsqrt(3.0d0)/2.0d0*(AA-BB)
          yc_3=-0.50d0*(AA+BB)-a/3.0d0-i*dsqrt(3.0d0)/2.0d0*(AA-BB)
!         print *, 'yc_1                        yc_2                    yc_3'
!         print *, yc_1, yc_2, yc_3
          CCC=999.0d0
          BBB=999.0d0
          AAA=999.0d0

        end if

        end


        subroutine piksrt(nn,arr)
        implicit none
        integer nn
        double precision arr(nn)
        integer i,j
        double precision a

        do j=2,nn
          a=arr(j)
          do i=j-1,1,-1
            if (arr(i) .le. a) goto 10
            arr(i+1)=arr(i)
          end do
          i=0
10        arr(i+1)=a
        end do

        end
