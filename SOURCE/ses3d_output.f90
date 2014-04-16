!==============================================================================
! store field values at lower order collocation points
!==============================================================================

subroutine ses3d_store(field,name,n,deg,tag)
use parameters
use variables
implicit none
include 'mpif.h'

        ! field - field variable
        ! name  - string containing the name of the field variable
        ! n     - iteration number
        ! deg   - polynomial degree of the stored field

        !======================================================================
        ! local variables
        !======================================================================

        integer, intent(in) :: n, deg, tag
        integer :: i, j, k, ip, iq, ir

        real :: cknots(0:7)
		real :: lgll, lgll_dummy
		real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(in) :: field
        real, dimension(0:nx,0:ny,0:nz,0:deg,0:deg,0:deg) :: cfield

        character(len=10) :: dummy, ns
        character(len=1000) :: fn
		character(len=*), intent(in) :: name

        !======================================================================
        ! compress only if deg<lpd
        !======================================================================

        if (deg<lpd) then

 		   !======================================================================
           ! determine collocation points (knots) and integration weights
           !======================================================================

           if (deg==2) then

			   cknots(0)=-1.0
			   cknots(1)=0.0
			   cknots(2)=1.0

           elseif (deg==3) then

			   cknots(0)=-1.0
			   cknots(1)=-0.4472135954999579
			   cknots(2)=0.4472135954999579
			   cknots(3)=1.0

           elseif (deg==4) then

			   cknots(0)=-1.0
			   cknots(1)=-0.6546536707079772
			   cknots(2)=0.0
			   cknots(3)=0.6546536707079772
			   cknots(4)=1.0

           elseif (deg==5) then

			   cknots(0)=-1.0
			   cknots(1)=-0.7650553239294647
			   cknots(2)=-0.2852315164806451
			   cknots(3)=0.2852315164806451
			   cknots(4)=0.7650553239294647
			   cknots(5)=1.0

           elseif (deg==6) then

			   cknots(0)=-1.0
			   cknots(1)=-0.8302238962785670
			   cknots(2)=-0.4688487934707142
			   cknots(3)=0.0
			   cknots(4)=0.4688487934707142
			   cknots(5)=0.8302238962785670
			   cknots(6)=1.0

           elseif (deg==7) then

			   cknots(0)=-1.0
			   cknots(1)=-0.8717401485096066
			   cknots(2)=-0.5917001814331423
			   cknots(3)=-0.2092992179024789
			   cknots(4)=0.2092992179024789
			   cknots(5)=0.5917001814331423
			   cknots(6)=0.8717401485096066
			   cknots(7)=1.0

           endif

           !======================================================================
           ! generate lower degree interpolant
           !======================================================================

           cfield(:,:,:,:,:,:)=0.0

           do i=0,deg
           do j=0,deg
           do k=0,deg

              	do ip=0,lpd
                do iq=0,lpd
				do ir=0,lpd

					lgll_dummy=lgll(lpd,ip,cknots(i))*lgll(lpd,iq,cknots(j))*lgll(lpd,ir,cknots(k))
                	cfield(0:nx,0:ny,0:nz,i,j,k)=cfield(0:nx,0:ny,0:nz,i,j,k)+field(0:nx,0:ny,0:nz,ip,iq,ir)*lgll_dummy

               	enddo
                enddo
              	enddo

           enddo
           enddo
           enddo

        endif

        !======================================================================
        ! store field variable
        !======================================================================

        call int2str(my_rank,dummy)

        if (n==1) then
			
			fn=ffd(1:len_trim(ffd))//name(1:len_trim(name))//'_'//dummy(1:len_trim(dummy))

           	write(99,*) fn

           	open(unit=tag,file=fn,action='write',form='unformatted')

        endif

        if (deg<lpd) then

           	write(tag) cfield

        else

           	write(tag) field

        endif

        if (n==nt) then

           close(unit=tag)

        endif

end subroutine ses3d_store


!==============================================================================
! read field values at lower order collocation points and interpolate
! to higher order collocation points
!==============================================================================

subroutine ses3d_restore(field,name,idx,deg,tag)
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! local variables
	!======================================================================

    integer, intent(in) :: idx, deg, tag
    integer :: i, j, k, il, im, in

	real :: lgll, lgll_dummy
   	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(inout) :: field
	real, dimension(0:nx,0:ny,0:nz,0:deg,0:deg,0:deg) :: cfield

  	character(len=*), intent(in) :: name
   	character(len=10) :: dummy
  	character(len=1000) :: fn

	!======================================================================
	! read lower degree polynomial collocation point values
	!======================================================================

 	call int2str(my_rank,dummy)

	!- open file in first iteration

  	if (idx==1) then

	   	fn=ffd(1:len_trim(ffd))//name(1:len_trim(name))//'_'//dummy(1:len_trim(dummy))

      	open(unit=tag,file=fn,action='read',form='unformatted',position='append')
       	backspace(unit=tag)

  	endif

 	!- read and check if interpolation is necessary

  	if (deg<lpd) then
		
		read(tag) cfield
        backspace(unit=tag)
	   	if (idx<nt) then
			backspace(unit=tag)
	   	endif

	else

       	read(tag) field
       	backspace(unit=tag)
		if (idx<nt) then
          	backspace(unit=tag)
		endif

	endif

 	!- close file in last iteration

    if (idx==nt) then

  	  	close(unit=tag)

  	endif

	!======================================================================
	! interpolate to higher polynomial collocation point values if necessary
	!======================================================================

  	if (deg<lpd) then

           field=0.0

           do i=0,lpd
           do j=0,lpd
           do k=0,lpd

              	do il=0,deg
               	do im=0,deg
                do in=0,deg

						lgll_dummy=lgll(deg,il,knots(i))*lgll(deg,im,knots(j))*lgll(deg,in,knots(k))
                       	field(0:nx,0:ny,0:nz,i,j,k)=field(0:nx,0:ny,0:nz,i,j,k)+cfield(0:nx,0:ny,0:nz,il,im,in)*lgll_dummy

                enddo
                enddo
              	enddo

         enddo
         enddo
         enddo

 	endif

end subroutine ses3d_restore


subroutine ses3d_output(n)
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! local variables
	!======================================================================

	integer, intent(in) :: n
   	integer :: i, j, k
	character(len=10) :: dummy, ns
	character(len=1000) :: fn
	character(len=12) :: sname

	!======================================================================
	! write displacement fields and preconditioner to files
	!======================================================================

	call int2str(my_rank,dummy)
	call int2str(n,ns)

	if (output_displacement==1) then

		fn=ofd(1:len_trim(ofd))//'vx_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
		open(unit=10,file=fn,action='write',form='unformatted')

		fn=ofd(1:len_trim(ofd))//'vy_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
		open(unit=11,file=fn,action='write',form='unformatted')

		fn=ofd(1:len_trim(ofd))//'vz_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
      	open(unit=12,file=fn,action='write',form='unformatted')

		write(10) vx
		write(11) vy
		write(12) vz

		close(unit=10)
		close(unit=11)
		close(unit=12)

	endif

	!======================================================================
	! write Frechet derivatives
	!======================================================================

	if (adjoint_flag==2) then

		!- Kernels for elastic parameters.=============================

		fn=ofd(1:len_trim(ofd))//'grad_rho_'//dummy(1:len_trim(dummy))
		open(unit=10,file=fn,action='write',form='unformatted')

		fn=ofd(1:len_trim(ofd))//'grad_cp_'//dummy(1:len_trim(dummy))
    	open(unit=12,file=fn,action='write',form='unformatted')

		fn=ofd(1:len_trim(ofd))//'grad_csh_'//dummy(1:len_trim(dummy))
		open(unit=14,file=fn,action='write',form='unformatted')

		fn=ofd(1:len_trim(ofd))//'grad_csv_'//dummy(1:len_trim(dummy))
		open(unit=15,file=fn,action='write',form='unformatted')

		write(10) grad_rho
		write(12) grad_cp
		write(14) grad_csh
		write(15) grad_csv

		close(unit=10)
		close(unit=12)
		close(unit=14)
		close(unit=15)

		!- Kernels for visco-elastic parameters. ======================

		if (is_diss==1) then

		   	fn=ofd(1:len_trim(ofd))//'grad_Q_mu_'//dummy(1:len_trim(dummy))
		   	open(unit=10,file=fn,action='write',form='unformatted')

			fn=ofd(1:len_trim(ofd))//'grad_alpha_mu_'//dummy(1:len_trim(dummy))
		   	open(unit=11,file=fn,action='write',form='unformatted')

			fn=ofd(1:len_trim(ofd))//'grad_Q_kappa_'//dummy(1:len_trim(dummy))
		   	open(unit=12,file=fn,action='write',form='unformatted')

			fn=ofd(1:len_trim(ofd))//'grad_alpha_kappa_'//dummy(1:len_trim(dummy))
		   	open(unit=13,file=fn,action='write',form='unformatted')

			write(10) grad_Q_mu
			write(11) grad_alpha_mu
			write(12) grad_Q_kappa
			write(13) grad_alpha_kappa

			close(unit=10)
			close(unit=11)
			close(unit=12)
			close(unit=13)

		endif


	endif

  	!======================================================================
   	! write saving vector to files if forward adjoint simulation
   	!======================================================================

	if (adjoint_flag==1) then

		fn=ffd(1:len_trim(ffd))//'saving_vector_'//dummy(1:len_trim(dummy))
        open(unit=10,file=fn,action='write')

        do k=1,nt

           	write(10,*) saving_vector(k)

        enddo

        close(unit=10)

	endif

	!======================================================================
	! write seismograms to files
	!======================================================================

	if ((nr>0) .and. (adjoint_flag<2)) then

	do i=1,nr

		sname=station_name_local(i)

		open(unit=51,file=ofd(1:len_trim(ofd))//sname//'.x',action='write')
		open(unit=52,file=ofd(1:len_trim(ofd))//sname//'.y',action='write')
		open(unit=53,file=ofd(1:len_trim(ofd))//sname//'.z',action='write')

		write(51,*) 'theta component seismograms'
		write(52,*) 'phi component seismograms'
		write(53,*) 'r component seismograms'

		write(51,*) 'nt=', nt
		write(51,*) 'dt=', dt

		write(52,*) 'nt=', nt
		write(52,*) 'dt=', dt

		write(53,*) 'nt=', nt
		write(53,*) 'dt=', dt

		write(51,*) 'receiver location (colat [deg],lon [deg],depth [m])'
		write(52,*) 'receiver location (colat [deg],lon [deg],depth [m])'
		write(53,*) 'receiver location (colat [deg],lon [deg],depth [m])'

		write(51,*) 'x=', recloc(1,i)*180/pi, 'y=', recloc(2,i)*180/pi, 'z=', recloc(3,i)
     	write(52,*) 'x=', recloc(1,i)*180/pi, 'y=', recloc(2,i)*180/pi, 'z=', recloc(3,i)
		write(53,*) 'x=', recloc(1,i)*180/pi, 'y=', recloc(2,i)*180/pi, 'z=', recloc(3,i)

		write(51,*) 'source location (colat [deg],lon [deg],depth [m])'
		write(52,*) 'source location (colat [deg],lon [deg],depth [m])'
		write(53,*) 'source location (colat [deg],lon [deg],depth [m])'

		write(51,*) 'x=', xxs*180/pi, 'y=', yys*180/pi, 'z=', zzs
		write(52,*) 'x=', xxs*180/pi, 'y=', yys*180/pi, 'z=', zzs
		write(53,*) 'x=', xxs*180/pi, 'y=', yys*180/pi, 'z=', zzs

		do j=1,nt

			write(51,*) seismogram_x(i,j)
			write(52,*) seismogram_y(i,j)
			write(53,*) seismogram_z(i,j)

		enddo

		close(unit=51)
		close(unit=52)
		close(unit=53)

	enddo

	endif

end subroutine ses3d_output


subroutine record_seismograms
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! local variables
	!======================================================================

 	integer :: idx
   	integer :: i, j, k
   	real :: dummy
	real :: lgll

    !======================================================================
    ! interpolate displacement field
    !======================================================================

	if (nr>0) then

		do idx=1,nr

			do i=0,lpd
			do j=0,lpd
			do k=0,lpd

				dummy=lgll(lpd,i,recloc_std(1,idx))*lgll(lpd,j,recloc_std(2,idx))*lgll(lpd,k,recloc_std(3,idx))

				seismogram_x(idx,it)=seismogram_x(idx,it)+vx(rx(idx),ry(idx),rz(idx),i,j,k)*dummy
	                       	seismogram_y(idx,it)=seismogram_y(idx,it)+vy(rx(idx),ry(idx),rz(idx),i,j,k)*dummy
	                       	seismogram_z(idx,it)=seismogram_z(idx,it)+vz(rx(idx),ry(idx),rz(idx),i,j,k)*dummy

			enddo
			enddo
			enddo

		enddo

	endif

end subroutine record_seismograms
