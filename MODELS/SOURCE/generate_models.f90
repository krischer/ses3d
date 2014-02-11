!=============================================================================
! Model generation for ses3d.
! Last modified: 5 February 2014 by Andreas Fichtner
!=============================================================================

program generate_models
use variables_gm
implicit none
include 'mpif.h'

	!======================================================================
	! initialise mpi and open logfile
	!======================================================================

	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_COMM_WORLD,p,ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)

	call int2str(my_rank,mrc)
	open(unit=99,file='logfile'//junk(1:len_trim(junk)),action='write')

	!======================================================================
	! Make sure all required directories are in existance.
	! Let's hope this works for all compilers...
	!======================================================================
    if (my_rank==0) then
        !- Create the necessary directories.
        call system('mkdir -p ../../DATA/COORDINATES')
        call system('mkdir -p ../MODELS')
        ! Block other processes until the directory has been created.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    else
        ! Block other processes until the directory has been created.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endif

	!======================================================================
	! process 0 reads input and broadcasts it to the other processes
	!======================================================================

	if (my_rank==0) then

		write(*,*) 'process 0 reading setupt file'

		write(99,*) '----------------------------------'
		write(99,*) '- setup file ---------------------'
		write(99,*) '----------------------------------'

		open (UNIT=15,FILE='../../INPUT/setup',STATUS='OLD',ACTION='READ')

		!- model geometry, anisotropy and dissipation ----------------

		read(15,*)junk
		read(15,*)x_min
		write(99,*) 'global minimum theta-extension: xmin=',x_min
		read(15,*)x_max
		write(99,*) 'global maximum theta-extension: xmax=',x_max
		read(15,*)y_min
		write(99,*) 'global minimum phi-extension: ymin=',y_min
		read(15,*)y_max
		write(99,*) 'global maximum phi-extension: ymax=',y_max
		read(15,*)z_min
		write(99,*) 'global minimum z-extension: zmin=',z_min
		read(15,*)z_max
		write(99,*) 'global maximum z-extension: zmax=',z_max
		read(15,*)is_diss
		write(99,*) 'dissipation on: ', is_diss
		read(15,*)model_type
		write(99,*) 'model type: ', model_type

		x_min=x_min*pi/180
    	x_max=x_max*pi/180

    	y_min=y_min*pi/180
    	y_max=y_max*pi/180

		!- computational setup, parallelisation ----------------------

		read(15,*)junk
		read(15,*) nx
		write(99,*) 'total number of elements in x-direction:', nx
		read(15,*) ny
		write(99,*) 'total number of elements in y-direction:', ny
		read(15,*) nz
		write(99,*) 'total number of elements in z-direction:', nz
		read(15,*) lpd
		write(99,*) 'LAGRANGE polynomial degree:', lpd
		read(15,*) px
		write(99,*) 'number of processors in x direction:', px
		read(15,*) py
		write(99,*) 'number of processors in y direction:', py
		read(15,*) pz
		write(99,*) 'number of processors in z direction:', pz

		write(99,*) '----------------------------------'

		close(15)

	endif

	call mpi_bcast(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(lpd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(x_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(x_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(y_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(y_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(z_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(z_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(model_type,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(is_diss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(px,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(py,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(pz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	!======================================================================
	! determine collocation points (knots)
	!======================================================================

	if (lpd==2) then

		knots(0)=-1.0
		knots(1)=0.0
		knots(2)=1.0

	elseif (lpd==3) then

		knots(0)=-1.0
		knots(1)=-0.4472135954999579
		knots(2)=0.4472135954999579
		knots(3)=1.0

	elseif (lpd==4) then

		knots(0)=-1.0
		knots(1)=-0.6546536707079772
		knots(2)=0.0
		knots(3)=0.6546536707079772
		knots(4)=1.0

	elseif (lpd==5) then

		knots(0)=-1.0
		knots(1)=-0.7650553239294647
		knots(2)=-0.2852315164806451
		knots(3)=0.2852315164806451
		knots(4)=0.7650553239294647
		knots(5)=1.0

	elseif (lpd==6) then

		knots(0)=-1.0
		knots(1)=-0.8302238962785670
		knots(2)=-0.4688487934707142
		knots(3)=0.0
		knots(4)=0.4688487934707142
		knots(5)=0.8302238962785670
		knots(6)=1.0

	elseif (lpd==7) then

		knots(0)=-1.0
		knots(1)=-0.8717401485096066
		knots(2)=-0.5917001814331423
		knots(3)=-0.2092992179024789
		knots(4)=0.2092992179024789
		knots(5)=0.5917001814331423
		knots(6)=0.8717401485096066
		knots(7)=1.0

	endif

	!======================================================================
	! make boxfile containing the parameters for parallel computing
	!======================================================================

	if (my_rank==0) then

		write(*,*) 'writing boxfile'

	endif

	allocate(min_index(1:px*py*pz,1:3),STAT=status)		! minimum element indeces
	allocate(max_index(1:px*py*pz,1:3),STAT=status)		! maximum element indeces

	dz=(z_max-z_min)/(nz+pz)					! width of the elements in z direction
	dy=(y_max-y_min)/(ny+py)					! width of the elements in y direction
	dx=(x_max-x_min)/(nx+px)					! width of the elements in x direction

	bnx=floor(real((nx+1)/px))				! number of elements per processor in x direction
	bny=floor(real((ny+1)/py))				! number of elements per processor in y direction
	bnz=floor(real((nz+1)/pz))				! number of elements per processor in z direction

	if (my_rank==0) then

		open(unit=11,file='../MODELS/boxfile',action='write')

		write(11,*) "- format description -"
		write(11,*) "total number of processes"
		write(11,*) "number of processes in x direction"
		write(11,*) "number of processes in y direction"
		write(11,*) "number of processes in z direction"
		write(11,*) "single index (npz-1)*py*px+(npy-1)*px+npx"
		write(11,*) "multi index (npx, npy, npz)"
		write(11,*) "index boundaries x (x_min, x_max)"
		write(11,*) "index boundaries y (y_min, y_max)"
		write(11,*) "index boundaries z (z_min, z_max)"
		write(11,*) "physical boundaries x"
		write(11,*) "physical boundaries y"
		write(11,*) "phyiscal boundaries z"

		write(11,*) '-------------------'
		write(11,*) px*py*pz
		write(11,*) px
		write(11,*) py
		write(11,*) pz
		write(11,*) '-------------------'

	endif

	do k=1,pz,1
		do j=1,py,1
			do i=1,px,1

				! indices (multiple and single) of the processor boxes (mi=my_rank+1)
				mi=(k-1)*py*px+(j-1)*px+i

				if (my_rank==0) then

					write(11,*) mi
					write(11,*) i,j,k

				endif

				! index boundaries of the processor boxes
				min_index(mi,1)=(i-1)*bnx
				min_index(mi,2)=(j-1)*bny
				min_index(mi,3)=(k-1)*bnz

				if (i/=px) then
					max_index(mi,1)=i*bnx
				else
					max_index(mi,1)=nx
				endif

				if (j/=py) then
					max_index(mi,2)=j*bny
				else
					max_index(mi,2)=ny
				endif

				if (k/=pz) then
					max_index(mi,3)=k*bnz
				else
					max_index(mi,3)=nz
				endif

				if (my_rank==0) then

					write(11,*) min_index(mi,1), max_index(mi,1)
					write(11,*) min_index(mi,2), max_index(mi,2)
					write(11,*) min_index(mi,3), max_index(mi,3)

					! physical boundaries of the processor boxes
					write(11,*) x_min+(min_index(mi,1)+i-1)*dx, x_min+(max_index(mi,1)+i)*dx
					write(11,*) y_min+(min_index(mi,2)+j-1)*dy, y_min+(max_index(mi,2)+j)*dy
					write(11,*) z_min+(min_index(mi,3)+k-1)*dz, z_min+(max_index(mi,3)+k)*dz
					write(11,*) '-------------------'

				endif

				if (my_rank==(mi-1)) then

					nx_loc=max_index(my_rank+1,1)-min_index(my_rank+1,1)
					ny_loc=max_index(my_rank+1,2)-min_index(my_rank+1,2)
					nz_loc=max_index(my_rank+1,3)-min_index(my_rank+1,3)

					x_min_loc=x_min+(min_index(my_rank+1,1)+i-1)*dx
					x_max_loc=x_min+(max_index(my_rank+1,1)+i)*dx

					y_min_loc=y_min+(min_index(my_rank+1,2)+j-1)*dy
					y_max_loc=y_min+(max_index(my_rank+1,2)+j)*dy

					z_min_loc=z_min+(min_index(my_rank+1,3)+k-1)*dz
					z_max_loc=z_min+(max_index(my_rank+1,3)+k)*dz

				endif

			enddo
		enddo
	enddo

	if (my_rank==0) then

		close(unit=11)

	endif

	!======================================================================
	! make coordinates and write coordinate files
	!======================================================================

	write(99,*) '- indices ----------------------'
	write(99,*) nx_loc, ny_loc, nz_loc
	write(99,*) '- phys. model dimensions -------'
	write(99,*) x_min_loc*180/pi, x_max_loc*180/pi
	write(99,*) y_min_loc*180/pi, y_max_loc*180/pi
	write(99,*) z_min_loc, z_max_loc

	allocate(x(0:nx_loc,0:lpd))
	allocate(y(0:ny_loc,0:lpd))
	allocate(z(0:nz_loc,0:lpd))

	allocate(XX(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(YY(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(ZZ(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))

	open(unit=101,file='../../DATA/COORDINATES/xco_'//mrc(1:len_trim(mrc)),status='unknown')
	open(unit=102,file='../../DATA/COORDINATES/yco_'//mrc(1:len_trim(mrc)),status='unknown')
	open(unit=103,file='../../DATA/COORDINATES/zco_'//mrc(1:len_trim(mrc)),status='unknown')

	do i=0,nx_loc
		do n=0,lpd
			x(i,n)=x_min_loc+i*dx+0.5*(1+knots(n))*dx
			XX(i,:,:,n,:,:)=x_min_loc+i*dx+0.5*(1+knots(n))*dx
			write(101,*) x(i,n)
		enddo
	enddo

	do j=0,ny_loc
		do n=0,lpd
			y(j,n)=y_min_loc+j*dy+0.5*(1+knots(n))*dy
			YY(:,j,:,:,n,:)=y_min_loc+j*dy+0.5*(1+knots(n))*dy
			write(102,*) y(j,n)
		enddo
	enddo

	do k=0,nz_loc
		do n=0,lpd
			z(k,n)=z_max_loc-k*dz-0.5*(1+knots(n))*dz
			ZZ(:,:,k,:,:,n)=z_max_loc-k*dz-0.5*(1+knots(n))*dz
			write(103,*) z(k,n)
		enddo
	enddo

	close(unit=101)
	close(unit=102)
	close(unit=103)

	XX=XX*180/pi
	YY=YY*180/pi

	!======================================================================
	! allocate memory and generate models
	!======================================================================

	allocate(rhoinv(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(mu(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(lambda(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(A(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(B(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(C(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(Q(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))

    allocate(rho(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(cs(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))
	allocate(cp(0:nx_loc,0:ny_loc,0:nz_loc,0:lpd,0:lpd,0:lpd))

	!=============================================================
	! implement various 1D Earth models
	!=============================================================

	!- homogeneous model -----------------------------------------

	if (model_type==1) then
		
		call homogeneous
		
	!- isotropic version of prem ---------------------------------
	
	elseif (model_type==2) then
	
		call prem_iso
		
	!- homogeneous model plus QL6 Q ------------------------------

	elseif (model_type==3) then
		
		call homogeneous_plus_Q
	
	!- EUROPEAN BACKGROUND MODEL ---------------------------------

    elseif (model_type==4) then

    	call eumod_bg

	!- AK135 -----------------------------------------------------

	elseif (model_type==7) then

		call ak135

	endif

	!======================================================================
	! save to files
	!======================================================================

	!==============================================================
	! open files, unformatted
	!==============================================================

	call int2str(my_rank,mrc)

	open(unit=11,file='../MODELS/rhoinv'//mrc(1:len_trim(mrc)),action='write',form='unformatted')
	open(unit=12,file='../MODELS/mu'//mrc(1:len_trim(mrc)),action='write',form='unformatted')
	open(unit=13,file='../MODELS/lambda'//mrc(1:len_trim(mrc)),action='write',form='unformatted')
	open(unit=14,file='../MODELS/A'//mrc(1:len_trim(mrc)),action='write',form='unformatted')
	open(unit=15,file='../MODELS/B'//mrc(1:len_trim(mrc)),action='write',form='unformatted')
	open(unit=16,file='../MODELS/C'//mrc(1:len_trim(mrc)),action='write',form='unformatted')

	if (is_diss==1) then	
		open(unit=17,file='../MODELS/Q'//mrc(1:len_trim(mrc)),action='write',form='unformatted')
	endif

	!==============================================================
	! write parameters to files, formatted or unformatted
	!==============================================================

	write(11) rhoinv(:,:,:,:,:,:)
	write(12) mu(:,:,:,:,:,:)
	write(13) lambda(:,:,:,:,:,:)
	write(14) A(:,:,:,:,:,:)
	write(15) B(:,:,:,:,:,:)
	write(16) C(:,:,:,:,:,:)

	if (is_diss==1) then
		write(17) Q(:,:,:,:,:,:)
	endif

	close(unit=11)
	close(unit=12)
	close(unit=13)
	close(unit=14)
	close(unit=15)
	close(unit=16)

	if (is_diss==1) then
		close(unit=17)
	endif

	!======================================================================
	! make a profile
	!======================================================================


	open(unit=11,file='../MODELS/prof_rhoinv'//mrc(1:len_trim(mrc)),action='write')
	open(unit=12,file='../MODELS/prof_lambda'//mrc(1:len_trim(mrc)),action='write')
	open(unit=13,file='../MODELS/prof_mu'//mrc(1:len_trim(mrc)),action='write')
	open(unit=21,file='../MODELS/prof_a'//mrc(1:len_trim(mrc)),action='write')
	open(unit=22,file='../MODELS/prof_b'//mrc(1:len_trim(mrc)),action='write')
	open(unit=23,file='../MODELS/prof_c'//mrc(1:len_trim(mrc)),action='write')

	if (is_diss==1) then
		open(unit=24,file='../MODELS/prof_Q'//mrc(1:len_trim(mrc)),action='write')
	endif

        do i=0,nz_loc
        	do k=0,lpd

			write(11,*) rhoinv(0,0,i,0,0,k), z(i,k)
			write(12,*) lambda(0,0,i,0,0,k), z(i,k)
			write(13,*) mu(0,0,i,0,0,k), z(i,k)
			write(21,*) A(0,0,i,0,0,k), z(i,k)
			write(22,*) B(0,0,i,0,0,k), z(i,k)
			write(23,*) C(0,0,i,0,0,k), z(i,k)

			if (is_diss==1) then
				write(24,*) Q(0,0,i,0,0,k), z(i,k)
			endif

		enddo
	enddo

	close(unit=11)
	close(unit=12)
	close(unit=13)
	close(unit=21)
	close(unit=22)
	close(unit=23)

	if (is_diss==1) then
		close(unit=24)
	endif

	!=======================================================================
	! clean up
	!=======================================================================

	deallocate(rhoinv)
	deallocate(mu)
	deallocate(lambda)

	deallocate(A)
	deallocate(B)
	deallocate(C)

	deallocate(Q)

    deallocate(rho)
    deallocate(cp)
    deallocate(cs)

	deallocate(min_index)
	deallocate(max_index)

    deallocate(x)
    deallocate(y)
    deallocate(z)

	deallocate(XX)
	deallocate(YY)
	deallocate(ZZ)

	close(unit=99)

	call MPI_Finalize(ierr)

end program generate_models



!==============================================================================
! subroutine for integer to string conversion
!==============================================================================

subroutine int2str(value, string)
implicit none

	integer, intent(in) :: value
	character(len=*), intent(inout) :: string

	character(len=10) :: c
	integer :: k, n, new_value, is
	real :: e


	e=1e9
	is=0

	if (value==0) then
		string(1:1)='0'
		string(2:10)=' '
	else

		new_value=value

		do k=1,10
			c(k:k)=char(floor(new_value/e)+48)

			if ((floor(new_value/e)==0) .and. (is==0)) then
				n=k
			else
				is=1
			endif

			!write(*,*) c
			new_value=new_value-e*floor(new_value/e)
			e=e/10
			string(k:k)=' '

		enddo

		string(1:10-n)=c(n+1:10)

	endif

	if (len(string)>10) then
		string(11:len(string))=' '
	endif

return
end subroutine int2str
