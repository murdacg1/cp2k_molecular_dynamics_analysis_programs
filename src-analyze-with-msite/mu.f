      integer NP
      parameter (NP=61000)
      real*8 r(3,NP)
      character*2 sy(NP)
      character*2 syt
      real*8 p(3,NP)
      real*8 q(NP)
      integer i,j,k,n
      real*8 mu_q(3), mu_p(3), mu(3) ! dipole of the cell

      real*8 bohr2ang
      parameter (bohr2ang=0.529177249d0)
      real*8 au2debye
      parameter (au2debye=1.0d0/0.39343029d0)

      open(10,FILE='./run-01.xyz')
      open(11,FILE='../run-01_chrg.xyz')
      open(12,FILE='../run-01_dipl.xyz')

      open(13,FILE='DIPOLE_CELL_Q_gks_check.XYZ')
      open(14,FILE='DIPOLE_CELL_P_gks_check.XYZ')
      open(15,FILE='DIPOLE_CELL_gks_check.XYZ')

100   read(10,*,END=200) n
      read(10,*)
      do i=1,n
        read(10,*) sy(i),(r(k,i),k=1,3)
        do k=1,3
          r(k,i)=r(k,i)/bohr2ang
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

      do k=1,3
        mu_q(k)=0.0d0
        mu_p(k)=0.0d0
        mu(k)=0.0d0
        do i=1,n
          mu_q(k)=mu_q(k)+q(i)*r(k,i)
          mu_p(k)=mu_p(k)+p(k,i)
          mu(k)=mu_q(k)+mu_p(k)
        enddo
      enddo

      write(13,'(3e20.10)') (mu_q(k)*au2debye,k=1,3)
      write(14,'(3e20.10)') (mu_p(k)*au2debye,k=1,3)
      write(15,'(3e20.10)') (mu(k)*au2debye,k=1,3)

      goto 100

200   continue

      close(10)
      close(11)
      close(12)

      close(13)
      close(14)
      close(15)

      end
