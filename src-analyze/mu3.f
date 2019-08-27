      integer NP
      parameter (NP=61000)
      real*8 r(3,NP)
      character*2 sy(NP)
      character*2 syt
      real*8 p(3,NP)
      real*8 q(NP)
      integer i,j,k,n
      real*8 charge, r_coc(3), d(3), mu_q(3), mu_p(3), mu(3) ! center of charge, charge, dipole: of the cell

      real*8 bohr2ang
      parameter (bohr2ang=0.529177249d0)
      real*8 au2debye
      parameter (au2debye=1.0d0/0.39343029d0)

      integer f_avg,f,m,a,c
      real*8 cell(3)
      real*8 pot_mon

        read(5,*) f_avg
        read(5,*) f
        read(5,*) m
        read(5,*) a
        read(5,*) (cell(k),k=1,3)
        read(5,*) temperature
        read(5,*) pot_mon
        read(5,*) c

        write(6,*) '# f_avg = ', f_avg
        write(6,*) '# f = ', f
        write(6,*) '# m = ', m
        write(6,*) '# a = ', a
        write(6,*) '# cell = ', (cell(k),k=1,3)
        do k=1,3
          cell(k)=cell(k)/bohr2ang
        end do
        write(6,*) '# cell = ', (cell(k),k=1,3)
        write(6,*) '# temperature [Kelvin] = ', temperature
        write(6,*) '# pot_mon [Hartree] = ', pot_mon
        write(6,*) '# number of bins c = ', c
        write(6,*)
        write(6,*)
        call flush(6)

      open(10,FILE='./run-01.xyz')
      open(11,FILE='../run-01_chrg.xyz')
      open(12,FILE='../run-01_dipl.xyz')

      open(16,FILE='CHARGE3_CELL.XYZ_temp')
      open(13,FILE='DIPOLE3_CELL_Q.XYZ_temp')
      open(14,FILE='DIPOLE3_CELL_P.XYZ_temp')
      open(15,FILE='DIPOLE3_CELL.XYZ_temp')

! for each frame
100   read(10,*,END=200) n
      read(10,*)
      do i=1,n
        read(10,*) sy(i),(r(k,i),k=1,3)
        do k=1,3
          r(k,i) = r(k,i)/bohr2ang
!         r(k,i) = r(k,i) - cell(k)*anint(r(k,i)/cell(k))
        end do
      enddo
      read(11,*)
      read(11,*)
      do i=1,n
        read(11,*) syt,q(i)
      enddo
      read(12,*)
      read(12,*)
      do i=1,n
        read(12,*) syt,(p(k,i),k=1,3)
      enddo

! total charge of cell
      charge=0.0d0
      do i=1,n
        charge=charge+q(i)
      enddo

! center of charge of cell
      do k=1,3
        r_coc(k)=0.0d0
        do i=1,n
          r_coc(k)=r_coc(k)+q(i)*r(k,i)
        enddo
        r_coc(k) = r_coc(k)/charge
!       r_coc(k) = r_coc(k) - cell(k)*anint(r_coc(k)/cell(k))
      enddo

      do k=1,3
        mu_q(k)=0.0d0
        mu_p(k)=0.0d0
        mu(k)=0.0d0
        do i=1,n
          d(k) = r(k,i) - r_coc(k)
!         d(k) = d(k) - cell(k)*anint(d(k)/cell(k))
          mu_q(k)=mu_q(k)+q(i)*d(k) ! relative to center of charge
          mu_p(k)=mu_p(k)+p(k,i)
          mu(k)=mu_q(k)+mu_p(k)
        enddo
      enddo

      write(16,'(1e20.10)') charge
      write(13,'(3e20.10)') (mu_q(k)*au2debye,k=1,3)
      write(14,'(3e20.10)') (mu_p(k)*au2debye,k=1,3)
      write(15,'(3e20.10)') (mu(k)*au2debye,k=1,3)

      goto 100
! for each frame

200   continue

      close(10)
      close(11)
      close(12)

      close(16)
      close(13)
      close(14)
      close(15)

      end
