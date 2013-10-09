!==============================================================================
! read and distribute input
!==============================================================================

subroutine ses3d_input
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! variables
	!======================================================================

    character(len=15) :: junk

	integer :: pp, dummy
    integer :: status(MPI_STATUS_SIZE)
    integer :: nx_min_loc, nx_max_loc
    integer :: ny_min_loc, ny_max_loc
    integer :: nz_min_loc, nz_max_loc
    integer :: ix_multi_loc, iy_multi_loc, iz_multi_loc
    integer :: i, k

    real :: xmin_loc, xmax_loc
    real :: ymin_loc, ymax_loc
    real :: zmin_loc, zmax_loc

	!======================================================================
	! read box file, only process 0, and send to other processes
	!======================================================================

	write(99,*) ' '
	write(99,*) 'begin input'
	write(99,*) '------------------------------------------------------------'
	write(99,*) 'read boxfile ***********************************************'

	if (my_rank==0) then

		write(*,*) '----------------------------------------'
		write(*,*) 'begin input'
		write(*,*) 'process 0 reading boxfile'
		open(unit=15,file='../MODELS/MODELS/boxfile',status='old',action='read')

		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk

		write(99,*) '- parallelisation parameters -------------------------------'

		read(15,*) pp
		write(99,*) 'total number of processors:', pp
		read(15,*) px
		write(99,*) 'number of processors in theta direction:', px
		read(15,*) py
		write(99,*) 'number of processors in phi direction:', py
		read(15,*) pz
		write(99,*) 'number of processors in z direction:', pz

		if (pp/=p) then

			write(*,*) 'incorrect number of processes'
			write(*,*) 'inconsistency in boxfile and s_mpirun'
			write(99,*) 'ERROR: incorrect number of processes'
			write(99,*) 'ERROR: inconsistency in boxfile and s_mpirun'

		endif

		write(*,*) 'number of processors:', pp
		read(15,*) junk

		do k=1,p,1

			read(15,*) dummy
			read(15,*) ix_multi_loc, iy_multi_loc, iz_multi_loc
			read(15,*) nx_min_loc, nx_max_loc
			read(15,*) ny_min_loc, ny_max_loc
			read(15,*) nz_min_loc, nz_max_loc
			read(15,*) xmin_loc, xmax_loc
			read(15,*) ymin_loc, ymax_loc
			read(15,*) zmin_loc, zmax_loc
			read(15,*) junk

			if (k==1) then

				ix_multi=ix_multi_loc
				iy_multi=iy_multi_loc
				iz_multi=iz_multi_loc

				nx=nx_max_loc-nx_min_loc
				ny=ny_max_loc-ny_min_loc
				nz=nz_max_loc-nz_min_loc

				xmin=xmin_loc
				xmax=xmax_loc
				ymin=ymin_loc
				ymax=ymax_loc
				zmin=zmin_loc
				zmax=zmax_loc

				write(99,*)
				write(99,*) '- geometrical parameters ----------------------------------'
				write(99,*)
				write(99,*) 'itheta_multi=',ix_multi,' iphi_multi=',iy_multi,' iz_multi=',iz_multi
				write(99,*) 'n_theta=',nx,' n_phi=',ny,' n_z=',nz
				write(99,*) 'theta_min=',xmin*180/pi,' theta_max=',xmax*180/pi
				write(99,*) 'phi_min=',ymin*180/pi,' phi_max=',ymax*180/pi
				write(99,*) 'z_min=',zmin,' z_max=',zmax

			else

				call MPI_Send(ix_multi_loc, 1, MPI_INTEGER, k-1, 1, MPI_COMM_WORLD, ierr)
				call MPI_Send(iy_multi_loc, 1, MPI_INTEGER, k-1, 2, MPI_COMM_WORLD, ierr)
				call MPI_Send(iz_multi_loc, 1, MPI_INTEGER, k-1, 3, MPI_COMM_WORLD, ierr)
				call MPI_Send(nx_max_loc-nx_min_loc, 1, MPI_INTEGER, k-1, 4, MPI_COMM_WORLD, ierr)
				call MPI_Send(ny_max_loc-ny_min_loc, 1, MPI_INTEGER, k-1, 5, MPI_COMM_WORLD, ierr)
				call MPI_Send(nz_max_loc-nz_min_loc, 1, MPI_INTEGER, k-1, 6, MPI_COMM_WORLD, ierr)
				call MPI_Send(xmin_loc, 1, MPI_REAL, k-1, 7, MPI_COMM_WORLD, ierr)
				call MPI_Send(xmax_loc, 1, MPI_REAL, k-1, 8, MPI_COMM_WORLD, ierr)
				call MPI_Send(ymin_loc, 1, MPI_REAL, k-1, 9, MPI_COMM_WORLD, ierr)
				call MPI_Send(ymax_loc, 1, MPI_REAL, k-1, 10, MPI_COMM_WORLD, ierr)
				call MPI_Send(zmin_loc, 1, MPI_REAL, k-1, 11, MPI_COMM_WORLD, ierr)
				call MPI_Send(zmax_loc, 1, MPI_REAL, k-1, 12, MPI_COMM_WORLD, ierr)

			endif

		enddo

		close(unit=15)

	else

		call MPI_Recv(ix_multi, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(iy_multi, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(iz_multi, 1, MPI_INTEGER, 0, 3, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(nx, 1, MPI_INTEGER, 0, 4, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(ny, 1, MPI_INTEGER, 0, 5, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(nz, 1, MPI_INTEGER, 0, 6, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(xmin, 1, MPI_REAL, 0, 7, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(xmax, 1, MPI_REAL, 0, 8, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(ymin, 1, MPI_REAL, 0, 9, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(ymax, 1, MPI_REAL, 0, 10, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(zmin, 1, MPI_REAL, 0, 11, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(zmax, 1, MPI_REAL, 0, 12, MPI_COMM_WORLD, status, ierr)

		write(99,*) '- geometrical parameters ----------------------------------'
		write(99,*) 'itheta_multi=',ix_multi,' iphi_multi=',iy_multi,' iz_multi=',iz_multi
		write(99,*) 'n_theta=',nx,' n_phi=',ny,' n_z=',nz
		write(99,*) 'theta_min=',xmin*180/pi,' theta_max=',xmax*180/pi
		write(99,*) 'phi_min=',ymin*180/pi,' phi_max=',ymax*180/pi
		write(99,*) 'z_min=',zmin,' z_max=',zmax

	endif

	call mpi_bcast(p,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(px,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(py,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(pz,1,mpi_integer,0,mpi_comm_world,ierr)

	write(99,*) '- parallelisation parameters -------------------------------'
	write(99,*) 'total number of processors:', p
	write(99,*) 'number of processors in theta direction:', px
	write(99,*) 'number of processors in phi direction:', py
	write(99,*) 'number of processors in z direction:', pz
	write(99,*) '------------------------------------------------------------'

	!======================================================================
	! read adequate event, setup and relax files, only process 0
	!======================================================================

	if (my_rank==0) then

		!================================================================================================
		! read event file
		!================================================================================================

		write(*,*) 'process 0 reading event_'//trim(event_indices(i_events))//' file'
		write(99,*) 'read event file **********************************************'
		open (unit=15,file='../INPUT/event_'//trim(event_indices(i_events)),status='old',action='read')

		!- time stepping parameters ----------------------------------

		write(99,*) '- time stepping parameters ---------------------------------'
		read(15,*)junk
		read(15,*)nt
		write(99,*) 'number of time steps: nt=', nt
		read(15,*)dt
		write(99,*) 'time increment: dt=', dt

		!- source parameters -----------------------------------------

		write(99,*) '- source parameters ----------------------------------------'
		read(15,*)junk
		read(15,*)xxs
		write(99,*) 'theta-coordinate of the source: xxs=',xxs
		read(15,*)yys
		write(99,*) 'phi-coordinate of the source: yys=',yys
		read(15,*)zzs
		write(99,*) 'z-coordinate of the source (depth): zzs=',zzs
		read(15,*)source_type
		write(99,*) 'source_type: ',source_type
                read(15,*)MOM_xx
                write(99,*) 'Mxx: ', MOM_xx
                read(15,*)MOM_yy
                write(99,*) 'Myy: ', MOM_yy
                read(15,*)MOM_zz
                write(99,*) 'Mzz: ', MOM_zz
                read(15,*)MOM_xy
                write(99,*) 'Mxy: ', MOM_xy
                read(15,*)MOM_xz
                write(99,*) 'Mxz: ', MOM_xz
                read(15,*) MOM_yz
                write(99,*) 'Myz: ', MOM_yz

		!- output directory ------------------------------------------

		write(99,*) '- output directory -----------------------------------------'
		read(15,*)junk
		read(15,'(A100)')ofd
		write(99,*) 'output field directory: ', ofd
        !- Create the directory if it does not exist.
        call system('mkdir -p '//ofd)

		!- output flags ----------------------------------------------

		write(99,*) '- output flags ---------------------------------------------'
		read(15,*)junk
                read(15,*)ssamp
		write(99,*) 'output rate: ssamp=',ssamp
		read(15,*)output_displacement
		write(99,*) 'output_displcement=', output_displacement

		close(unit=15)

		!================================================================================================
		! read setup file
		!================================================================================================

		write(*,*) 'process 0 reading setup file'
		write(99,*) 'read setup file **********************************************'
		open (unit=15,file='../INPUT/setup',status='old',action='read')

		!- model geometry and dissipation -----------------

		write(99,*) '- model geometry, anisotropy and dissipation ---------------'
		read(15,*)junk
		read(15,*)xmin_global
		write(99,*) 'global minimum theta-extension: xmin=',xmin_global
		read(15,*)xmax_global
		write(99,*) 'global maximum theta-extension: xmax=',xmax_global
		read(15,*)ymin_global
		write(99,*) 'global minimum phi-extension: ymin=',ymin_global
		read(15,*)ymax_global
		write(99,*) 'global maximum phi-extension: ymax=',ymax_global
		read(15,*)zmin_global
		write(99,*) 'global minimum z-extension: zmin=',zmin_global
		read(15,*)zmax_global
		write(99,*) 'global maximum z-extension: zmax=',zmax_global
		read(15,*)is_diss
		write(99,*) 'dissipation on: ', is_diss
		read(15,*)integer_dummy
		write(99,*) 'model type:', integer_dummy

		!- computational setup, parallelisation -----------------------

		write(99,*) '- computational setup, parallelisation ---------------------'
		read(15,*)junk
		read(15,*) integer_dummy
		write(99,*) 'total number of elements in x-direction:', integer_dummy
		read(15,*) integer_dummy
		write(99,*) 'total number of elements in y-direction:', integer_dummy
		read(15,*) integer_dummy
		write(99,*) 'total number of elements in z-direction:', integer_dummy
		read(15,*) integer_dummy
		write(99,*) 'LAGRANGE polynomial degree:', integer_dummy
		read(15,*) integer_dummy
		write(99,*) 'number of processors in x direction:', integer_dummy
		read(15,*) integer_dummy
		write(99,*) 'number of processors in y direction:', integer_dummy
		read(15,*) integer_dummy
		write(99,*) 'number of processors in z direction:', integer_dummy

		!- adjoint flags ---------------------------------------------

		write(99,*) '- adjoint flags --------------------------------------------'
		read(15,*)junk
		read(15,*)adjoint_flag
		write(99,*) 'adjoint_flag=', adjoint_flag
		read(15,*)samp_ad
		write(99,*) 'forward field storing rate:', samp_ad
		read(15,'(A100)') ffd
		ffd=ffd(1:len_trim(ffd))//trim(event_indices(i_events))//'/'
		write(99,*) 'forward field directory: ', ffd
        !- Create the directory if it does not exist.
        call system('mkdir -p '//ffd)

		write(99,*) '------------------------------------------------------------'

		close(unit=15)

                xxs=xxs*pi/180
                yys=yys*pi/180

                xmin_global=xmin_global*pi/180
                xmax_global=xmax_global*pi/180
                ymin_global=ymin_global*pi/180
                ymax_global=ymax_global*pi/180

		!================================================================================================
		! read relax file
		!================================================================================================

		write(*,*) 'process 0 reading relax file'
		write(99,*) 'read relax file **********************************************'
		open (unit=15,file='../INPUT/relax',status='old',action='read')

		!- relaxation times

		write(99,*) '- relaxation times -----------------------------------------'
		read(15,*) junk
		read(15,*) tau_p(1)
		read(15,*) tau_p(2)
		read(15,*) tau_p(3)
		write(99,*) tau_p(1), tau_p(2), tau_p(3)

		!- weights of relaxation mechnisms

		write(99,*) '- weights of relaxation mechanisms -------------------------'
		read(15,*) junk
		read(15,*) D_p(1)
		read(15,*) D_p(2)
		read(15,*) D_p(3)
		write(99,*) D_p(1), D_p(2), D_p(3)

		sum_D_p=sum(D_p)

		write(99,*) '------------------------------------------------------------'

		close(unit=15)

	endif

	!=====================================================================
	! broadcast input parameters to the other processes
	!=====================================================================

	call mpi_bcast(nt, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(dt, 1, mpi_real, 0, mpi_comm_world, ierr)

	call mpi_bcast(xxs, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(yys, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zzs, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(source_type, 1, mpi_integer, 0, mpi_comm_world, ierr)

    call mpi_bcast(MOM_xx, 1, mpi_real, 0, mpi_comm_world, ierr)
    call mpi_bcast(MOM_yy, 1, mpi_real, 0, mpi_comm_world, ierr)
    call mpi_bcast(MOM_zz, 1, mpi_real, 0, mpi_comm_world, ierr)
    call mpi_bcast(MOM_xy, 1, mpi_real, 0, mpi_comm_world, ierr)
    call mpi_bcast(MOM_xz, 1, mpi_real, 0, mpi_comm_world, ierr)
    call mpi_bcast(MOM_yz, 1, mpi_real, 0, mpi_comm_world, ierr)

	call mpi_bcast(xmin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(xmax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(ymin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(ymax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zmin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zmax_global, 1, mpi_real, 0, mpi_comm_world, ierr)

	call mpi_bcast(is_diss, 1, mpi_integer, 0, mpi_comm_world, ierr)

	call mpi_bcast(ofd,100,mpi_character,0,mpi_comm_world,ierr)

	call mpi_bcast(ssamp, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(output_displacement, 1, mpi_integer, 0, mpi_comm_world, ierr)

	call mpi_bcast(adjoint_flag,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(samp_ad, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(ffd,100,mpi_character,0,mpi_comm_world,ierr)

	call mpi_bcast(tau_p,3,mpi_real,0,mpi_comm_world,ierr)
	call mpi_bcast(D_p,3,mpi_real,0,mpi_comm_world,ierr)

	!=============================================================================
	! read source time function, adequate process only
	!=============================================================================

	if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

		if ((xmin .le. xxs) .and. (xmax .ge. xxs) .AND. &
		    (ymin .le. yys) .and. (ymax .ge. yys) .and. &
		    (zmin .le. zmax_global-zzs) .and. (zmax .ge. zmax_global-zzs)) then

			write(99,*) 'process ', my_rank, 'reading source time function'
                	write(*,*) 'process', my_rank, 'speaking ...'
			is_source=1.0

			open(unit=15,file='../INPUT/stf',status='old',action='read')

			do i=1, nt, 1
				read(15,*) so(i)
			end do

			close(unit=15)

		else

			write(99,*) 'source not located in this processor box'
			is_source=0.0

		endif

		call mpi_reduce(is_source*my_rank,source_processor,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

		if (my_rank==0) then

                	call int2str(int(source_processor),junk)

			write(*,*) 'process', junk, 'reading source time function'
                	write(99,*) 'process', source_processor, 'reading source time function'

		endif

	endif

	!=============================================================================
	! Read structural information (rho, mu, lambda, A, B, C), formatted or unformatted
	!=============================================================================

	call int2str(my_rank,junk)

	if (my_rank==0) then
		write(*,*) "reading structural information"
	endif

	write(99,*) '------------------------------------------------------------'
	write(99,*) 'read structural input **************************************'
	write(99,*) '------------------------------------------------------------'

	if (is_diss==1) then
		write(99,*) 'Dissipation on'
	else
		write(99,*) 'Dissipation off'
	endif

	open(unit=15,file='../MODELS/MODELS/rhoinv'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')
	open(unit=16,file='../MODELS/MODELS/mu'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')
	open(unit=17,file='../MODELS/MODELS/lambda'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')
	open(unit=21,file='../MODELS/MODELS/A'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')
	open(unit=22,file='../MODELS/MODELS/B'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')
	open(unit=23,file='../MODELS/MODELS/C'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')

	read(15) rhoinv(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	read(16) mu(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	read(17) lambda(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	read(21) A(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	read(22) B(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	read(23) C(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)

	close(unit=15)
	close(unit=16)
	close(unit=17)
	close(unit=21)
	close(unit=22)
	close(unit=23)

	if (is_diss==1) then

		open (unit=21,file='../MODELS/MODELS/Q'//junk(1:len_trim(junk)),status='old',action='read',form='unformatted')
		read(21) QQ(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
		close(unit=21)

		tau=2/(pi*QQ)

	endif

	write(99,*) 'minimum rho: ',minval(1/rhoinv(:,:,:,:,:,:)), ' maximum rho: ',maxval(1/rhoinv(:,:,:,:,:,:))
	write(99,*) 'minimum lambda: ',minval(lambda(:,:,:,:,:,:)), ' maximum lambda: ',maxval(lambda(:,:,:,:,:,:))
	write(99,*) 'minimum mu: ',minval(mu(:,:,:,:,:,:)), ' maximum mu: ', maxval(mu(:,:,:,:,:,:))
	write(99,*) 'minimum A: ',minval(A(:,:,:,:,:,:)), ' maximum A: ', maxval(A(:,:,:,:,:,:))
	write(99,*) 'minimum B: ',minval(B(:,:,:,:,:,:)), ' maximum B: ', maxval(B(:,:,:,:,:,:))
	write(99,*) 'minimum C: ',minval(C(:,:,:,:,:,:)), ' maximum C: ', maxval(C(:,:,:,:,:,:))

	if (is_diss==1) then

	       write(99,*) 'minimum tau: ',minval(tau(:,:,:,:,:,:)), ' maximum tau: ', maxval(tau(:,:,:,:,:,:))

	endif

	write(99,*) '------------------------------------------------------------'

	!======================================================================
	! read receiver locations, process  0 only
	!======================================================================

	if ((adjoint_flag==1) .or. (adjoint_flag==0)) then

		write(99,*) 'read receiver locations ************************************'

		if (my_rank==0) then

			write(*,*) 'process 0 reading receiver locations'
			write(99,*) '- receiver locations (colat [deg], lon [deg], depth [m] ----'

			open(unit=10, file='../INPUT/recfile_'//trim(event_indices(i_events)),status='old',action='read')

			read(10,*) nr_global	!number of receivers

			do i=1,nr_global

				read(10,*) station_name_global(i)
				read(10,*) recloc_global(1,i), recloc_global(2,i), recloc_global(3,i)
				write(99,*) station_name_global(i), recloc_global(1,i), recloc_global(2,i), recloc_global(3,i)

                        	recloc_global(1,i)=recloc_global(1,i)*pi/180
                        	recloc_global(2,i)=recloc_global(2,i)*pi/180

			enddo

			close(10)

			write(99,*) '------------------------------------------------------------'

		endif

		call mpi_bcast(station_name_global, 12*maxnr, mpi_character, 0, mpi_comm_world, ierr)
		call mpi_bcast(recloc_global, 3*maxnr, mpi_real, 0, mpi_comm_world, ierr)
		call mpi_bcast(nr_global,1,mpi_integer,0,mpi_comm_world,ierr)

	endif

	!=============================================================================
        ! read saving vector
        !=============================================================================

	if (adjoint_flag==2) then

		call int2str(my_rank,junk)

        	write(99,*) 'read saving vector'


		open(unit=10,file=ffd(1:len_trim(ffd))//'saving_vector_'//junk(1:len_trim(junk)),action='read')

        	do k=1,nt

           		read(10,*) saving_vector(k)

        	enddo

        	close(unit=10)

	endif

	!=============================================================================
        ! read adjoint source locations, process 0 only, and broadcast to others
        !=============================================================================

	if (adjoint_flag==2) then

        	if (my_rank==0) then

			write(99,*) 'read adjoint source locations *******************************'

           		!- read adjoint source locations ------------------------------------------

           		open(unit=10,file='../ADJOINT/'//trim(event_indices(i_events))//'/ad_srcfile',action='read')

			read(10,*) nr_adsrc_global

			do k=1,nr_adsrc_global

              			read(10,*) ad_srcloc_global(1,k), ad_srcloc_global(2,k), ad_srcloc_global(3,k)

              			ad_srcloc_global(1,k)=ad_srcloc_global(1,k)*pi/180
              			ad_srcloc_global(2,k)=ad_srcloc_global(2,k)*pi/180

              			write(99,*) ad_srcloc_global(1,k)*180/pi, ad_srcloc_global(2,k)*180/pi, ad_srcloc_global(3,k)
              			write(*,*) ad_srcloc_global(1,k)*180/pi, ad_srcloc_global(2,k)*180/pi, ad_srcloc_global(3,k)

           		enddo

			write(99,*) 'number of adjoint sources: ', nr_adsrc_global
           		write(*,*) 'number of adjoint sources: ', nr_adsrc_global

           	close(unit=10)

        	endif

        	call mpi_bcast(nr_adsrc_global, 1, mpi_integer, 0, mpi_comm_world, ierr)
        	call mpi_bcast(ad_srcloc_global, 3*maxnr, mpi_real, 0, mpi_comm_world, ierr)

	endif

	!=============================================================================
        ! finish
        !=============================================================================

	write(99,*) 'end input'
	write(99,*) '------------------------------------------------------------'

end subroutine ses3d_input
