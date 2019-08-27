        PROGRAM binz

        USE common_data_mod
        USE reading_mod
        USE binning_z_mod

! Purpose: Process a CP2K z_s umbrella sampling simulation of (H20)n slab with a single solute molecule (last molecule in the xyz file)
! to obtain the histogram of the solute z (COM = M3) and later line up the histograms to obtain the PMF

! Program is run as:
!       ../binz.exe < binz.in >& binz.out &
!
! It is assumed that the file run-01.xyz have the atoms organized as O,H,H(,M) ; O,H,H(,M); ...; N,O,O(,M),M3

! Program reads first the file binz.in:

! 450           f_avg (the final f frames are used for averaging) (avg calculated from last 20ps of the run)
! 900           f (total number of frames, not including the zeroth frame) (total run=1.0fs*30k=30ps)
! 215           m_w (number of water molecules per frame or cell)
! 1             m_s (number of solute molecules per frame or cell)
! 4             a_w (number of atoms per water molecule, last one is the M site if more than 3)
! 4             a_s (number of atoms per solute molecule, last one is the M site if more than 3)
! 1             which (which atom coordinate to do the binning on)
! 15.0 15.0 71.44  cell ABC values [Angstrom]
! 0.0           z_s (z position of the umbrella external potential) [ANGSTROM]
! 15.99491463     1.0078250321  1.0078250321  0.0  mass_w(:)
! 14.0030740052  15.99491463   15.99491463    0.0  mass_s(:)
! 35.72         middle (middle of the slab) [Angstrom]
! 0.1           width (bin width for the histogram) [Angstrom]
! 1.0           separation (separation between neighboring runs / umbrella potentials / histograms) [Angstrom]

        IMPLICIT NONE

        INTEGER :: my_d
        CHARACTER(len=80) my_file

        WRITE(6,'(a)') '# PROGRAM = binz'
        CALL FLUSH(6)

        READ(5,*) f_avg
        READ(5,*) f
        f0 = f - f_avg
        READ(5,*) m_w
        READ(5,*) m_s
        READ(5,*) a_w
        READ(5,*) a_s
        READ(5,*) which
        ALLOCATE( cell(d) )
        cell = zero
        READ(5,*) cell(:)
        READ(5,*) z_s
        ALLOCATE( mass_w(a_w) )
        mass_w = zero
        READ(5,*) mass_w(:)
        ALLOCATE( mass_s(a_s) )
        mass_s = zero
        READ(5,*) mass_s(:)
        READ(5,*) middle
        READ(5,*) width
        READ(5,*) separation

        WRITE(6,200) '# f_avg, f, f0+1 = ', f_avg, f, f0+1
        WRITE(6,100) '# m_w, m_s = ', m_w, m_s
        WRITE(6,100) '# a_w, a_s = ', a_w, a_s
        WRITE(6,105) '# which = ', which
        WRITE(6,310) '# cell = ', cell(:)
        WRITE(6,110) '# z_s = ', z_s
        WRITE(6,410) '# mass_w, SUM(mass_w) = ', &
                        mass_w, SUM(mass_w)
        WRITE(6,410) '# mass_s, SUM(mass_s) = ', &
                        mass_s, SUM(mass_s)
        WRITE(6,110) '# middle [Angstrom] = ', middle
        WRITE(6,110) '# width [Angstrom] = ', width
        WRITE(6,110) '# separation [Angstrom] = ', separation
        WRITE(6,*)
        WRITE(6,*)
        CALL FLUSH(6)

        my_file = '../run-01.xyz'
        my_file = TRIM(my_file)
        my_d = d
        ALLOCATE( z(0:f) )
        CALL reading(my_d, z, &
          my_file) ! reading and getting com of NO2 for each frame

        my_file='z.bin'
        my_d=1
        CALL binning_z(z, &
          my_file)

        DEALLOCATE( cell )
        DEALLOCATE( mass_w )
        DEALLOCATE( mass_s )
        DEALLOCATE( z )

105     FORMAT(1x,a40,1(1x,i10))
100     FORMAT(1x,a40,2(1x,i10))
200     FORMAT(1x,a40,3(1x,i10))
110     FORMAT(1x,a40,1(1x,e20.10))
310     FORMAT(1x,a40,3(1x,e20.10))
410     FORMAT(1x,a40,5(1x,e20.10))
510     FORMAT(1x,a20,1(1x,i10),3(1x,f10.4))
520     FORMAT(1x,a10,3(17x,i1))
530     FORMAT(1x,a10,1(1x,i10),3(5x,4(1x,f6.1)))
540     FORMAT(1x,a40,5(1x,e20.10))
600     FORMAT(4(1x,f10.1))

        END PROGRAM binz
