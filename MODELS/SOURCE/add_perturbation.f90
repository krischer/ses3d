!================================================================================
!= add 3D variations to existing Earth model
!= last modified on 18 September 2012 by Andreas Fichtner
!================================================================================

program ses3d_add_perturbation
use parameters
use variables
implicit none
include 'mpif.h'

	!-----------------------------------------------------------------------
	!- local variables
	!-----------------------------------------------------------------------

	character(len=20) :: string_rank
	character(len=15) :: junk    ! dummy variable 
    
	integer :: pp, dummy
	integer :: status(MPI_STATUS_SIZE)
    
	integer :: nx_min_loc, nx_max_loc
	integer :: ny_min_loc, ny_max_loc
	integer :: nz_min_loc, nz_max_loc
    
	integer :: ix_multi_loc, iy_multi_loc, iz_multi_loc
    
	integer :: i, j, k, l, m, n, nsubvol, isubvol
	
	integer :: nbx, nby, nbz
	integer :: start_x, start_y, start_z, end_x, end_y, end_z

	real :: xmin_loc, xmax_loc
	real :: ymin_loc, ymax_loc
	real :: zmin_loc, zmax_loc
	
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: cox, coy, coz
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: csh, csv
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: d_mu, d_b, d_lambda
	
	real, allocatable, dimension(:) :: bx, by, bz
	real, allocatable, dimension(:,:,:) :: d_rho, d_beta, d_beta_sh, d_beta_sv, d_alpha
	
	!-----------------------------------------------------------------------
	!- initializations
	!-----------------------------------------------------------------------
	
	call mpi_init(ierr)
	call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
	call mpi_comm_size(mpi_comm_world, p, ierr)
	
	call int2str(my_rank, string_rank)
		
	!-----------------------------------------------------------------------
	! read box file, only process 0, and send to other processes
	!-----------------------------------------------------------------------
	
	if (my_rank==0) then
		
		write(*,*) '----------------------------------------'
		write(*,*) 'begin input'
		write(*,*) 'process 0 reading boxfile' 
		open(unit=15,file='../MODELS/boxfile',status='old',action='read')
	
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
	
		read(15,*) pp
		read(15,*) px
		read(15,*) py
		read(15,*) pz
		
		if (pp/=p) then
			
			write(*,*) 'incorrect number of processes'
			write(*,*) 'inconsistency in boxfile and s_mpirun'
			
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
		
	endif

	call mpi_bcast(p,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(px,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(py,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(pz,1,mpi_integer,0,mpi_comm_world,ierr)
	
	!-----------------------------------------------------------------------
	! read Par_1 file, only process 0, and send to others
	!-----------------------------------------------------------------------

	if (my_rank==0) then

		write(*,*) 'process 0 reading setup file'
		
		open(unit=15,file='../../INPUT/setup',status='old',action='read') 

		read(15,*)junk
		read(15,*)xmin_global
		read(15,*)xmax_global
		read(15,*)ymin_global
		read(15,*)ymax_global
		read(15,*)zmin_global
		read(15,*)zmax_global
		read(15,*)is_diss

		close(15)

                xmin_global=xmin_global*pi/180
                xmax_global=xmax_global*pi/180
                ymin_global=ymin_global*pi/180
                ymax_global=ymax_global*pi/180
                
	endif
	
	call mpi_bcast(xmin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(xmax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(ymin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(ymax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zmin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zmax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(is_diss, 1, mpi_integer, 0, mpi_comm_world, ierr)
	
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

	!-----------------------------------------------------------------------
	!- make coordinate lines
	!-----------------------------------------------------------------------
	
	dx=(xmax-xmin)/(nx+1)		! width of one element in x direction
	dy=(ymax-ymin)/(ny+1)		! width of one element in y direction
	dz=(zmax-zmin)/(nz+1)		! width of one element in z direction
	
	do i=0,nx
		do n=0,lpd
			cox(i,:,:,n,:,:)=xmin+i*dx+0.5*(1+knots(n))*dx
		enddo
	enddo

	do j=0,ny
		do n=0,lpd
			coy(:,j,:,:,n,:)=ymin+j*dy+0.5*(1+knots(n))*dy
		enddo
	enddo
	
	do k=0,nz
		do n=0,lpd
			coz(:,:,k,:,:,n)=zmax-k*dz-0.5*(1+knots(n))*dz
		enddo
	enddo

	!-----------------------------------------------------------------------
	! read structural information from existing files (rho, mu, lambda, A, B, C)
	!-----------------------------------------------------------------------

	if (my_rank==0) then
		write(*,*) "reading structural information"
	endif
	
	open(unit=15,file='../MODELS/rhoinv'//string_rank(1:len_trim(string_rank)),status='old',action='read',form='unformatted')
	open(unit=16,file='../MODELS/mu'//string_rank(1:len_trim(string_rank)),status='old',action='read',form='unformatted')
	open(unit=17,file='../MODELS/lambda'//string_rank(1:len_trim(string_rank)),status='old',action='read',form='unformatted')
	open(unit=21,file='../MODELS/A'//string_rank(1:len_trim(string_rank)),status='old',action='read',form='unformatted')
	open(unit=22,file='../MODELS/B'//string_rank(1:len_trim(string_rank)),status='old',action='read',form='unformatted')
	open(unit=23,file='../MODELS/C'//string_rank(1:len_trim(string_rank)),status='old',action='read',form='unformatted')
			
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

	!-----------------------------------------------------------------------
	!- make S velocities, P velocity and density
	!-----------------------------------------------------------------------

	csh=sqrt(mu*rhoinv)
	csv=sqrt((mu+B)*rhoinv)
	cp=sqrt((lambda+2*mu)*rhoinv)
	rho(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)=1/rhoinv(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)

	!-----------------------------------------------------------------------
	!- read number of subvolumes
	!-----------------------------------------------------------------------
	
	open(unit=55,file='../MODELS_3D/block_x',action='read',status='old');
	open(unit=56,file='../MODELS_3D/block_y',action='read',status='old');
	open(unit=57,file='../MODELS_3D/block_z',action='read',status='old');
	
	read(55,*) nsubvol
	read(56,*) nsubvol
	read(57,*) nsubvol

	!-----------------------------------------------------------------------
	!- open model perturbations
	!-----------------------------------------------------------------------

	open(unit=14,file='../MODELS_3D/dvsh',action='read',status='old')
	open(unit=15,file='../MODELS_3D/dvsv',action='read',status='old')
	open(unit=17,file='../MODELS_3D/dvp',action='read',status='old')	
	open(unit=18,file='../MODELS_3D/drho',action='read',status='old')	

	read(14,*) nsubvol
	read(15,*) nsubvol
	read(17,*) nsubvol
	read(18,*) nsubvol

	!-----------------------------------------------------------------------
	!- loop over the subvolumes
	!-----------------------------------------------------------------------

	do isubvol=1,nsubvol

	  !- make coordinate lines ---------------------------------------------

	  read(55,*) nbx	
	  read(56,*) nby
	  read(57,*) nbz
	
	  allocate(bx(1:nbx))
	  allocate(by(1:nby))
	  allocate(bz(1:nbz))
	
	  do i=1,nbx
	    read(55,*) bx(i)
	    bx(i)=bx(i)*pi/180
	  enddo
	
	  do i=1,nby
	    read(56,*) by(i)
	    by(i)=by(i)*pi/180
	  enddo
	
	  do i=1,nbz
	    read(57,*) bz(i)
	    bz(i)=bz(i)*1000
	  enddo

	  !- allocate model perturbations --------------------------------------

	  allocate(d_rho(1:(nbx-1),1:(nby-1),1:(nbz-1)))
	  allocate(d_alpha(1:(nbx-1),1:(nby-1),1:(nbz-1)))
	  allocate(d_beta_sh(1:(nbx-1),1:(nby-1),1:(nbz-1)))
	  allocate(d_beta_sv(1:(nbx-1),1:(nby-1),1:(nbz-1)))

	  !- read model perturbations ------------------------------------------

	  read(14,*) dummy
	  read(15,*) dummy
	  read(17,*) dummy
	  read(18,*) dummy

	  do i=1,nbx-1
	  do j=1,nby-1
	  do k=1,nbz-1
	
	      read(14,*) d_beta_sh(i,j,k)
	      read(15,*) d_beta_sv(i,j,k)
	      read(17,*) d_alpha(i,j,k)
	      read(18,*) d_rho(i,j,k)
		
	  enddo
	  enddo
	  enddo

	  d_beta_sh=1000*d_beta_sh		! convert from km/s to m/s
	  d_beta_sv=1000*d_beta_sv		! convert from km/s to m/s
	  d_alpha=1000*d_alpha			! convert from km/s to m/s
	  d_rho=1000*d_rho                      ! convert from g/cm^3 to kg/m^3
	
	  !- add perturbation to the current model -----------------------------
	
	  start_x=1
	  start_y=1
	  start_z=1
	
	  end_x=nbx-1
	  end_y=nby-1
	  end_z=nbz-1
	
	  do i=1,nbx-1
		if (bx(i)<xmin) then
			start_x=i
		endif
	  enddo
	  do i=nbx-1,1,-1
		if (bx(i)>xmax) then
			end_x=i
		endif
	  enddo
	
	  do i=1,nby-1
		if (by(i)<ymin) then
			start_y=i
		endif
	  enddo
	  do i=nby-1,1,-1
		if (by(i)>ymax) then
			end_y=i
		endif
	  enddo
	
	  do i=1,nbz-1
		if (bz(i)<zmin) then
			start_z=i
		endif
	  enddo
	  do i=nbz-1,1,-1
		if (bz(i)>zmax) then
			end_z=i
		endif
	  enddo
	
	  do i=start_x,end_x
	  do j=start_y,end_y
	  do k=start_z,end_z
	  
	    where ((cox>=bx(i)) .and. (cox<bx(i+1)) .and. &
		   (coy>=by(j)) .and. (coy<by(j+1)) .and. &
		   (coz>bz(k)) .and. (coz<=bz(k+1)))

		      csv=csv+d_beta_sv(i,j,k)
		      csh=csh+d_beta_sh(i,j,k)
		      cp=cp+d_alpha(i,j,k)
		      rho=rho+d_rho(i,j,k)
	
	    end where

	  enddo
	  enddo
	  enddo

	  !- deallocate ------------------------------------------------------------

	  deallocate(bx)
	  deallocate(by)
	  deallocate(bz)

	  deallocate(d_beta_sh)
	  deallocate(d_beta_sv)
	  deallocate(d_alpha)
	  deallocate(d_rho)

	enddo  !- loop over subvolumes

	!- close files ---------------------------------------------------------

	close(unit=14)
	close(unit=15)
	close(unit=17)
	close(unit=18)

	close(unit=55)
	close(unit=56)
	close(unit=57)

	!-----------------------------------------------------------------------
	!- recompute elastic parameters
	!-----------------------------------------------------------------------
  
	mu=rho*csh*csh
	B=rho*csv*csv-mu
	lambda=rho*cp*cp-2*mu

	!-----------------------------------------------------------------------
	!- write new model to files
	!-----------------------------------------------------------------------

	rhoinv(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)=1/rho(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)		

	open(unit=11,file='../MODELS/rhoinv'//string_rank(1:len_trim(string_rank)),action='write',form='unformatted')
	open(unit=12,file='../MODELS/mu'//string_rank(1:len_trim(string_rank)),action='write',form='unformatted')
	open(unit=13,file='../MODELS/lambda'//string_rank(1:len_trim(string_rank)),action='write',form='unformatted')
	open(unit=14,file='../MODELS/A'//string_rank(1:len_trim(string_rank)),action='write',form='unformatted')
	open(unit=15,file='../MODELS/B'//string_rank(1:len_trim(string_rank)),action='write',form='unformatted')
	open(unit=16,file='../MODELS/C'//string_rank(1:len_trim(string_rank)),action='write',form='unformatted')
	
	write(11) rhoinv(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	write(12) mu(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	write(13) lambda(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	write(14) A(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	write(15) B(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
	write(16) C(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)

	close(UNIT=11)
	close(UNIT=12)
	close(UNIT=13)
	close(UNIT=14)
	close(UNIT=15)
	close(UNIT=16)

	!-----------------------------------------------------------------------
	!- clean up
	!-----------------------------------------------------------------------
	
	call mpi_finalize(ierr)

end program ses3d_add_perturbation




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

