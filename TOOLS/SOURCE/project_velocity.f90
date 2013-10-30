!===============================================================================
! function for the integration of LAGRANGE polynomials
!===============================================================================

real function int_lag(idx,deg,lim_left,lim_right)
implicit none

	!-----------------------------------------------------------------------
	!- idx=index of LAGRANGE polynomial
	!- deg=degree of LAGRANGE polynomial (lpd)
	!- lim_left=left integration limit
	!- lim_right=right integration limit
	!-----------------------------------------------------------------------

	!-----------------------------------------------------------------------
	!- variables
	!-----------------------------------------------------------------------

	integer, intent(in) :: idx, deg
	real, intent(in) :: lim_left, lim_right

	integer :: i, j, k

	real :: norm
	real, dimension(1:deg) :: z0
	real, dimension(0:deg) :: a, knots

	!-----------------------------------------------------------------------
	!- initialisations
	!-----------------------------------------------------------------------

	if (deg==2) then

		knots(0)=-1.0
		knots(1)=0.0
		knots(2)=1.0

	elseif (deg==3) then

		knots(0)=-1.0
		knots(1)=-0.4472135954999579
		knots(2)=0.4472135954999579
		knots(3)=1.0

	elseif (deg==4) then

		knots(0)=-1.0
		knots(1)=-0.6546536707079772
		knots(2)=0.0
		knots(3)=0.6546536707079772
		knots(4)=1.0

	elseif (deg==5) then

		knots(0)=-1.0
		knots(1)=-0.7650553239294647
		knots(2)=-0.2852315164806451
		knots(3)=0.2852315164806451
		knots(4)=0.7650553239294647
		knots(5)=1.0

	elseif (deg==6) then

		knots(0)=-1.0
		knots(1)=-0.8302238962785670
		knots(2)=-0.4688487934707142
		knots(3)=0.0
		knots(4)=0.4688487934707142
		knots(5)=0.8302238962785670
		knots(6)=1.0

	elseif (deg==7) then

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
	!- computations
	!-----------------------------------------------------------------------

	!- make vector of negative roots ---------------------------------------

	j=1
	do i=0,deg
		if (i .ne. idx) then
			z0(j)=-knots(i)
			j=j+1
		endif
	enddo

	!- compute polynomial coefficients -------------------------------------

	if (deg==1) then

		a(1)=1
		a(0)=z0(1)

	elseif (deg==2) then

		a(2)=1
		a(1)=z0(1)+z0(2)
		a(0)=z0(1)*z0(2)

	elseif (deg>2) then

		!- initialisation of the iteration at degree 2

		a(2)=1
		a(1)=z0(1)+z0(2)
		a(0)=z0(1)*z0(2)

		!- successively compute coefficients for higher polynomials

		do i=3,deg

			!- compute the remaining coefficients

			do k=(i-2),1,-1
				a(k)=a(k)*z0(i)+a(k-1)
			enddo

			!- set the trivial coefficients for degree i >= 3

			a(i)=1
			a(i-1)=0
			a(0)=1
			do k=1,i
				a(i-1)=a(i-1)+z0(k)
				a(0)=a(0)*z0(k)
			enddo


		enddo

	endif

	!- compute constant factor ---------------------------------------------

	norm=1.0

	do i=0,deg
		if (i .ne. idx) then
			norm=norm*(knots(idx)-knots(i))
		endif
	enddo

	!- evaluate the integral -----------------------------------------------

	int_lag=0

	do i=0,deg

		int_lag=int_lag+(lim_right**(i+1)-lim_left**(i+1))*a(i)/(i+1)

	enddo

	int_lag=int_lag/norm

end function int_lag

!==============================================================================
! function for integer to string conversion
!==============================================================================

character(len=10) function int2str(value)
implicit none

	integer, intent(in) :: value
	character(len=12) :: format_string

	format_string='(i'//achar(floor(log10(real(value)))+1+48)//')'

	write(int2str,format_string)  value

end function int2str

!===============================================================================
! projection onto the basis functions
!===============================================================================

program make_gradient
use parameters
use variables
implicit none

	!-0-==========================================================================
	! local variables
	!=============================================================================

	character(len=5) :: comp
	character(len=10) :: int2str
	character(len=60) :: junk, fn, cit
	character(len=140) :: fn_grad, fn_output

	integer :: nx_min_loc, nx_max_loc
	integer :: ny_min_loc, ny_max_loc
	integer :: nz_min_loc, nz_max_loc
	integer :: ix_multi_loc, iy_multi_loc, iz_multi_loc
	integer :: rank, ind

	integer :: nsubvol, isubvol

	integer :: i, j, k, l, m, n, indx, indy, indz

	integer :: ix_min, ix_max, iy_min, iy_max, iz_min, iz_max
	real :: res(1:2), dbx, dby, dbz

	real :: xlima, xlimb, ylima, ylimb, zlima, zlimb
	real :: intx, inty, intz
	real :: int_lag

	integer :: nbx, nby, nbz
	real, allocatable, dimension(:) :: bxco, byco, bzco
	real, allocatable, dimension(:,:,:) :: gradient_csv, gradient_csh, gradient_cp, gradient_rho

	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: u_csv, u_csh, u_cp, u_rho

	!-0-==========================================================================
	! initialisations
	!=============================================================================

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

	!-0-==========================================================================
	! read external input
	!=============================================================================

	write(*,*) 'iteration: '
	read(*,*) it

	write(*,*) 'directory for sensitivity densities: '
	read(*,*) fn_grad
	write(*,*) 'confirm directory: ', fn_grad

	write(*,*) 'directory for output: '
	read(*,*) fn_output
	write(*,*) 'confirm directory: ', fn_output

	cit=int2str(it)

	!-0-==========================================================================
	! open gradient file and read header block coordinates
	!=============================================================================

	open(unit=120,file=fn_output(1:len_trim(fn_output))//'gradient_csh',action='write')
	open(unit=130,file=fn_output(1:len_trim(fn_output))//'gradient_csv',action='write')
	open(unit=140,file=fn_output(1:len_trim(fn_output))//'gradient_cp',action='write')
	open(unit=150,file=fn_output(1:len_trim(fn_output))//'gradient_rho',action='write')

	open(unit=71,file='../MODELS/MODELS_3D/block_x',status='old',action='read')
	open(unit=72,file='../MODELS/MODELS_3D/block_y',status='old',action='read')
	open(unit=73,file='../MODELS/MODELS_3D/block_z',status='old',action='read')

	read(71,*) nsubvol
	read(72,*) nsubvol
	read(73,*) nsubvol

	write(120,*) nsubvol
	write(130,*) nsubvol
	write(140,*) nsubvol
	write(150,*) nsubvol

	!-0-==========================================================================
	! loop over subvolumes
	!=============================================================================

	do isubvol=1,nsubvol

		!-1-==================================================================
		! read individual subvolumes
		!=====================================================================

		read(71,*) nbx
		read(72,*) nby
		read(73,*) nbz

		allocate(gradient_csh(1:nbx-1,1:nby-1,1:nbz-1))
		allocate(gradient_csv(1:nbx-1,1:nby-1,1:nbz-1))
		allocate(gradient_cp(1:nbx-1,1:nby-1,1:nbz-1))
		allocate(gradient_rho(1:nbx-1,1:nby-1,1:nbz-1))

		allocate(bxco(1:nbx))
		allocate(byco(1:nby))
		allocate(bzco(1:nbz))

		gradient_csh(:,:,:)=0.0
		gradient_csv(:,:,:)=0.0
		gradient_cp(:,:,:)=0.0
		gradient_rho(:,:,:)=0.0

		do i=1,nbx
			read(71,*) bxco(i)
		enddo
		do j=1,nby
			read(72,*) byco(j)
		enddo
		do k=1,nbz
			read(73,*) bzco(k)
		enddo

		dbx=bxco(2)-bxco(1)
		dby=byco(2)-byco(1)
		dbz=bzco(2)-bzco(1)

		bxco=bxco*pi/180
		byco=byco*pi/180
		bzco=bzco*1000

		!-1-==================================================================
		! read boxfile header
		!=====================================================================

		write(*,*) '----------------------------------'

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

		read(15,*) p
		read(15,*) px
		read(15,*) py
		read(15,*) pz

		!-1-==================================================================
		! read individual box information
		!=====================================================================

		read(15,*) junk

		do ind=1,p

		      write(*,*) ind

		      !-2-============================================================
		      ! read boxfile entries for one processor
		      !===============================================================

		      read(15,*) rank
		      read(15,*) ix_multi_loc, iy_multi_loc, iz_multi_loc
		      read(15,*) nx_min_loc, nx_max_loc
		      read(15,*) ny_min_loc, ny_max_loc
		      read(15,*) nz_min_loc, nz_max_loc
		      read(15,*) xmin, xmax
		      read(15,*) ymin, ymax
		      read(15,*) zmin, zmax
		      read(15,*) junk

		      nx=nx_max_loc-nx_min_loc
		      ny=ny_max_loc-ny_min_loc
		      nz=nz_max_loc-nz_min_loc

		      junk=int2str(rank-1)

		      !-2-============================================================
		      ! make coordinates
		      !===============================================================

		      dx=(xmax-xmin)/(nx+1)		! width of one element in x direction
		      dy=(ymax-ymin)/(ny+1)		! width of one element in y direction
		      dz=(zmax-zmin)/(nz+1)		! width of one element in z direction

		      do i=0,nx
		      do n=0,lpd
			x(i,n)=xmin+i*dx+0.5*(1+knots(n))*dx
			sin_theta(i,:,:,n,:,:)=sin(x(i,n))
		      enddo
		      enddo

		      do j=0,ny
		      do n=0,lpd
			y(j,n)=ymin+j*dy+0.5*(1+knots(n))*dy
		      enddo
		      enddo

		      do k=0,nz
		      do n=0,lpd
			z(k,n)=zmax-k*dz-0.5*(1+knots(n))*dz
		      enddo
		      enddo

		      !-2-============================================================
		      ! load sensitivity kernels
		      !===============================================================

		      !- csv

		      fn=fn_grad(1:len_trim(fn_grad))//'grad_csv_'//junk(1:len_trim(junk))//'_'//cit(1:len_trim(cit))

		      write(*,*) 'open file ', fn
		      open(unit=10,file=fn,action='read',form='unformatted')
		      read(10) u_csv
		      close(unit=10)

		      !- csh

		      fn=fn_grad(1:len_trim(fn_grad))//'grad_csh_'//junk(1:len_trim(junk))//'_'//cit(1:len_trim(cit))

		      write(*,*) 'open file ', fn
		      open(unit=10,file=fn,action='read',form='unformatted')
		      read(10) u_csh
		      close(unit=10)

		      !- cp

		      fn=fn_grad(1:len_trim(fn_grad))//'grad_cp_'//junk(1:len_trim(junk))//'_'//cit(1:len_trim(cit))

		      write(*,*) 'open file ', fn
		      open(unit=10,file=fn,action='read',form='unformatted')
		      read(10) u_cp
		      close(unit=10)

		      !- rho

		      fn=fn_grad(1:len_trim(fn_grad))//'grad_rho_'//junk(1:len_trim(junk))//'_'//cit(1:len_trim(cit))

		      write(*,*) 'open file ', fn
		      open(unit=10,file=fn,action='read',form='unformatted')
		      read(10) u_rho
		      close(unit=10)

		      !-2-============================================================
		      ! loop over individual blocks
		      !===============================================================

		      do indx=1,nbx-1
		      do indy=1,nby-1
		      do indz=1,nbz-1

			    !-3-======================================================
			    !- check if processor overlaps with this model block
			    !=========================================================

			    if ((bxco(indx+1)>=xmin) .and. (bxco(indx)<=xmax) .and. &
				(byco(indy+1)>=ymin) .and. (byco(indy)<=ymax) .and. &
				(bzco(indz+1)>=zmin) .and. (bzco(indz)<=zmax)) then

				  !-4-================================================
				  ! find indices of integration limits
				  !===================================================

				  !- x-direction

				  if (bxco(indx)>xmin) then
				    res=minloc(abs(x-bxco(indx)))
				    ix_min=res(1)-1
				    if ((x(ix_min,lpd)==bxco(indx))) then
				      ix_min=ix_min+1
				    endif
				  else				! eventuell redundant
				    ix_min=0
				  endif

				  if (bxco(indx+1)<xmax) then
				    res=minloc(abs(x-bxco(indx+1)))
				    ix_max=res(1)-1
				    if (x(ix_max,0)==bxco(indx+1)) then
				      ix_max=ix_max-1
				    endif
				  else				! eventuell redundant
				    ix_max=nx
				  endif

				  !- y-direction

				  if (byco(indy)>ymin) then
				    res=minloc(abs(y-byco(indy)))
				    iy_min=res(1)-1
				    if ((y(iy_min,lpd)==byco(indy))) then
				      iy_min=iy_min+1
				    endif
				  else
				    iy_min=0
				  endif

				  if (byco(indy+1)<ymax) then
				    res=minloc(abs(y-byco(indy+1)))
				    iy_max=res(1)-1
				    if (y(iy_max,0)==byco(indy+1)) then
				      iy_max=iy_max-1
				    endif
				  else
				    iy_max=ny
				  endif

				  !- z-direction

				  if (bzco(indz)>zmin) then
				    res=minloc(abs(z-bzco(indz)))
				    iz_max=res(1)-1
				  else
				    iz_max=nz
				  endif

				  if (bzco(indz+1)<zmax) then
				    res=minloc(abs(z-bzco(indz+1)))
				    iz_min=res(1)-1
				  else
				    iz_min=0
				  endif

				  !-4-================================================
				  ! perform integration
				  !===================================================

				  do i=ix_min,ix_max	! loop over elements in this box in x direction

				    !- integration limits in x direction ---------------

				    if (bxco(indx)>x(i,0)) then
				      xlima=(2*bxco(indx)-x(i,lpd)-x(i,0))/dx
				    else
				      xlima=-1.0
				    endif

				    if (bxco(indx+1)<x(i,lpd)) then
				      xlimb=(2*bxco(indx+1)-x(i,lpd)-x(i,0))/dx
				    else
				      xlimb=1.0
				    endif

				 do j=iy_min,iy_max	! loop over elements in this box in y direction

				    !- integration limits in y direction ---------------

				    if (byco(indy)>y(j,0)) then
				      ylima=(2*byco(indy)-y(j,lpd)-y(j,0))/dy
				    else
				      ylima=-1.0
				    endif

				    if (byco(indy+1)<y(j,lpd)) then
				      ylimb=(2*byco(indy+1)-y(j,lpd)-y(j,0))/dy
				    else
				      ylimb=1.0
				    endif

				do k=iz_min,iz_max	! loop over elements in this box in z direction

				    !- integration limits in z direction ---------------

				    if (bzco(indz)>z(k,lpd)) then
				      zlimb=(z(k,0)+z(k,lpd)-2*bzco(indz))/dz
				    else
				      zlimb=1.0
				    endif

				    if (bzco(indz+1)<z(k,0)) then
				      zlima=(z(k,0)+z(k,lpd)-2*bzco(indz+1))/dz
				    else
				      zlima=-1.0
				    endif

				!- integration ---------------------------------------------

				do l=0,lpd
				  intx=int_lag(l,lpd,xlima,xlimb)
				do m=0,lpd
				  inty=int_lag(m,lpd,ylima,ylimb)
				do n=0,lpd
				  intz=int_lag(n,lpd,zlima,zlimb)

				  gradient_csv(indx,indy,indz)=gradient_csv(indx,indy,indz)+u_csv(i,j,k,l,m,n)*intx*inty*intz/(dbx*dby*dbz)
				  gradient_csh(indx,indy,indz)=gradient_csh(indx,indy,indz)+u_csh(i,j,k,l,m,n)*intx*inty*intz/(dbx*dby*dbz)
				  gradient_cp(indx,indy,indz)=gradient_cp(indx,indy,indz)+u_cp(i,j,k,l,m,n)*intx*inty*intz/(dbx*dby*dbz)
				  gradient_rho(indx,indy,indz)=gradient_rho(indx,indy,indz)+u_rho(i,j,k,l,m,n)*intx*inty*intz/(dbx*dby*dbz)

				enddo
				enddo
				enddo

				enddo
				enddo
				enddo

			  endif	! if box is in this processor

		      enddo	! loop over boxes in z direction
		      enddo	! loop over boxes in y direction
		      enddo	! loop over boxes in x direction

		enddo	! loop over processors

		close(unit=15)	! close boxfile

		!-1-==================================================================
		! write output
		!=====================================================================

		write(120,*) (nbx-1)*(nby-1)*(nbz-1)
		write(130,*) (nbx-1)*(nby-1)*(nbz-1)
		write(140,*) (nbx-1)*(nby-1)*(nbz-1)
		write(150,*) (nbx-1)*(nby-1)*(nbz-1)

		do i=1,nbx-1
		do j=1,nby-1
		do k=1,nbz-1

		  write(120,*) gradient_csh(i,j,k)
		  write(130,*) gradient_csv(i,j,k)
		  write(140,*) gradient_cp(i,j,k)
		  write(150,*) gradient_rho(i,j,k)

		enddo
		enddo
		enddo

		!-1-==================================================================
		! clean up
		!=====================================================================

		deallocate(gradient_csh)
		deallocate(gradient_csv)
		deallocate(gradient_cp)
		deallocate(gradient_rho)

		deallocate(bxco)
		deallocate(byco)
		deallocate(bzco)

	enddo	!- loop over subvolumes

	!-0-==========================================================================
        ! clean up
        !=============================================================================

	close(unit=120)
	close(unit=130)
	close(unit=140)
	close(unit=150)

end program make_gradient
