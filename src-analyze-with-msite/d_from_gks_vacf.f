        program d_from_gks_vacf

        implicit none
 
        integer f_avg_dummy, f_dummy, m_dummy, c_dummy
        real*8 pot_mon_dummy
        real*8 cell_dummy(3)

        integer a, j
        real*8 temperature
        real*8 mass(3)
        real*8 mass_total
        integer nt, idummy, i, nc
        real*8 dt, cx_velo_mol_0, cx_velo_mol
        real*8 d_vacf_mol, d_vacf_mol_1, d_vacf_mol_2

        real*8 zero
        parameter (zero=0.0d0)
        real*8 fs2ps
        parameter (fs2ps=1.0d-3)

        real*8 kb, amu2kg
        real*8 s2ps, m2ang
        parameter (kb=1.3806503d-23, amu2kg=1.66053886d-27)
        parameter (s2ps=1.0d12, m2ang=1.0d10)

        write(6,'(a)') '# program = d_from_gks_vacf'

        open(10,file='analyze.in',status='unknown')
        read(10,*) f_avg_dummy
        read(10,*) f_dummy
        read(10,*) m_dummy
        read(10,*) a
        read(10,*) (cell_dummy(j),j=1,3)
        read(10,*) temperature
        read(10,*) pot_mon_dummy
        read(10,*) c_dummy
        read(10,*) (mass(j),j=1,a)
        close(10)

        mass_total = zero
        do j = 1, a
          mass_total = mass_total + mass(j)
        end do

        open(10,file='cx_vel.param',status='unknown')
        read(10,*) nt,dt,idummy ! dt in ps
        close(10)

        d_vacf_mol = zero

!       do i=1,aint(nt/2.0d0)
        do i=0,300 ! nt
          read(5,*,end=100) cx_velo_mol
!!!       cx_velo_mol = cx_velo_mol/(fs2ps**2) ! convert from A**2/fs**2 to A**2/ps**2 !!! cx_velo_mol already in A**2/ps**2
          if (i .eq. 0) cx_velo_mol_0 = cx_velo_mol
          cx_velo_mol = cx_velo_mol/cx_velo_mol_0 ! normalize
          nc=i ! count
          d_vacf_mol = d_vacf_mol + cx_velo_mol ! integrate
        end do

100     continue

        d_vacf_mol_1 = cx_velo_mol_0*d_vacf_mol*dt/3.0d0
        mass_total = mass_total*amu2kg
        d_vacf_mol_2=((kb*(m2ang/s2ps)**2)*temperature/mass_total)*d_vacf_mol*dt

        write(6,*) 'mass_total, ((kb*(m2ang/s2ps)**2)*temperature/mass_total), dt = '
        write(6,*) mass_total, ((kb*(m2ang/s2ps)**2)*temperature/mass_total), dt
        write(6,900) '# nc  d_vacf_mol = ', nc, d_vacf_mol
        write(6,900) '# nc  d_vacf_mol_1 = ', nc, d_vacf_mol_1
        write(6,900) '# nc  d_vacf_mol_2 = ', nc, d_vacf_mol_2

        write(6,*)
        write(6,*)

900     format(1(1x,a),1(1x,i12),1(1x,e20.10))

        end
