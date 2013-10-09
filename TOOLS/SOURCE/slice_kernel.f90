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
	! read setup file
	!======================================================================

	open (unit=15,file='../INPUT/setup',status='OLD',action='READ')

        read(15,*)junk
	read(15,*)xmin_global
	write(*,*) 'global minimum theta-extension: theta_min=',xmin_global
	read(15,*)xmax_global
	write(*,*) 'global maximum theta-extension: theta_max=',xmax_global
	read(15,*)ymin_global
	write(*,*) 'global minimum phi-extension: phi_min=',ymin_global
	read(15,*)ymax_global
	write(*,*) 'global maximum phi-extension: phi_max=',ymax_global
	read(15,*)zmin_global
	write(*,*) 'global minimum z-extension: zmin=',zmin_global
	read(15,*)zmax_global
	write(*,*) 'global maximum z-extension: zmax=',zmax_global

	read(15,*)is_diss

	close(15)

	!======================================================================
	! read external input
	!======================================================================

	write(*,*) 'plane (theta=1, phi=2, r=3): '
	read(*,*) plane
	write(*,*) 'parameter (1=rho, 2=cp, 3=cs, 4=csh, 5=csv): '
	read(*,*) comp
	write(*,*) 'coordinate value (in metres for r and degrees of theta or phi): '
	read(*,*) val
	write(*,*) 'iteration: '
	read(*,*) it

        if ((plane==1) .or. (plane==2)) then

           val=val*pi/180

        endif

	!======================================================================
	! open and write slice logfile
	!======================================================================

	if (plane==1) then

		fn=dir(1:len_trim(dir))//'slice_log_x'

	elseif (plane==2) then

		fn=dir(1:len_trim(dir))//'slice_log_y'

	elseif (plane==3) then

		fn=dir(1:len_trim(dir))//'slice_log_z'

	endif

	open(unit=93,file=fn,action='write')

	write(93,*) 'geometry: colat_min, colat_max, lon_min, lon_max'
	write(93,*) xmin_global, xmax_global, ymin_global, ymax_global
	write(93,*) 'anisotropy, dissipation'
	write(93,*) is_diss
	write(93,*) '==================================================='

	!======================================================================
	! read boxfile
	!======================================================================

	write(*,*) '----------------------------------'

	open(unit=15,file='../MODELS/MODELS/boxfile',status='old',action='read ')

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

	!======================================================================
	! write slice logfile
	!======================================================================

	write(93,*) 'number of processors:'

	if (plane==1) then
		write(93,*) py*pz
	elseif (plane==2) then
		write(93,*) px*pz
	elseif (plane==3) then
		write(93,*) px*py
	endif

	write(93,*) 'processor indices:'

	!======================================================================
	! read individual box information
	!======================================================================

	read(15,*) junk

	do ind=1,p

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

		call int2str(rank-1,junk)
		call int2str(it,ns)

		index_n=0
		index_l=0

		!==============================================================
		! x=const plane
		!==============================================================

		if ((plane==1) .and. (val<=xmax) .and. (val>xmin)) then

			write(93,*) junk

			!======================================================
			! allocate slice
			!======================================================

			allocate(slice(0:ny,0:nz,0:lpd,0:lpd),stat=status)

			!======================================================
			! x-coordinate line
			!======================================================

			dx=(xmax-xmin)/(nx+1)		! width of one element in x direction

			do i=0,nx
				do n=0,lpd

					x(i,n)=xmin+i*dx+0.5*(1+knots(n))*dx

					if (abs(val-x(i,n))<abs(val-x(index_n,index_l))) then
						index_n=i
						index_l=n
					endif

				enddo
			enddo

			write(*,*) index_n, index_l, x(index_n,index_l)

			!======================================================
			! open files and transfer values
			!======================================================

			if (comp==1) then
				fn=dir(1:len_trim(dir))//'grad_rho_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==2) then
				fn=dir(1:len_trim(dir))//'grad_cp_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==3) then
				fn=dir(1:len_trim(dir))//'grad_cs_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==4) then
				fn=dir(1:len_trim(dir))//'grad_csh_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==5) then
				fn=dir(1:len_trim(dir))//'grad_csv_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			endif

			write(*,*) 'opening ', fn
			open(unit=10,file=fn,action='read',form='unformatted')
			read(10) u
			close(unit=10)

			slice(0:ny,0:nz,0:lpd,0:lpd)=u(index_n,0:ny,0:nz,index_l,0:lpd,0:lpd)

			!==============================================================
			! write slices to files
			!==============================================================

			if (comp==1) then
				fn=dir(1:len_trim(dir))//'slice_x_grad_rho_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==2) then
				fn=dir(1:len_trim(dir))//'slice_x_grad_cp_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==3) then
				fn=dir(1:len_trim(dir))//'slice_x_grad_cs_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==4) then
				fn=dir(1:len_trim(dir))//'slice_x_grad_csh_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==5) then
				fn=dir(1:len_trim(dir))//'slice_x_grad_csv_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			endif

			open(unit=10,file=fn,action='write')

			do ind_z=0,nz
				do lz=0,lpd
					do ind_y=0,ny
						do ly=0,lpd
							write(10,*) slice(ind_y,ind_z,ly,lz)
						enddo
					enddo
				enddo
			enddo

			close(unit=10)

			deallocate(slice,stat=status)

		!==============================================================
		! y=const plane
		!==============================================================

		elseif ((plane==2) .and. (val<=ymax) .and. (val>ymin)) then

			write(93,*) junk

			!======================================================
			! allocate slice
			!======================================================

			allocate(slice(0:nx,0:nz,0:lpd,0:lpd),stat=status)

			!======================================================
			! y-coordinate line
			!======================================================

			dy=(ymax-ymin)/(ny+1)		! width of one element in y direction

			do j=0,ny
				do n=0,lpd

					y(j,n)=ymin+j*dy+0.5*(1+knots(n))*dy

					if (abs(val-y(j,n))<abs(val-y(index_n,index_l))) then

						index_n=j
						index_l=n

					endif

				enddo
			enddo

			write(*,*) index_n, index_l, y(index_n,index_l)

			!======================================================
			! open files and transfer values
			!======================================================

			if (comp==1) then
				fn=dir(1:len_trim(dir))//'grad_rho_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==2) then
				fn=dir(1:len_trim(dir))//'grad_cp_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==3) then
				fn=dir(1:len_trim(dir))//'grad_cs_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==4) then
				fn=dir(1:len_trim(dir))//'grad_csh_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==5) then
				fn=dir(1:len_trim(dir))//'grad_csv_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			endif

			write(*,*) 'opening ', fn
			open(unit=10,file=fn,action='read',form='unformatted')
			read(10) u
			close(unit=10)

			slice(0:nx,0:nz,0:lpd,0:lpd)=u(0:nx,index_n,0:nz,0:lpd,index_l,0:lpd)

			!==============================================================
			! write slices to files
			!==============================================================

			if (comp==1) then
				fn=dir(1:len_trim(dir))//'slice_y_grad_rho_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==2) then
				fn=dir(1:len_trim(dir))//'slice_y_grad_cp_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==3) then
				fn=dir(1:len_trim(dir))//'slice_y_grad_cs_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==4) then
				fn=dir(1:len_trim(dir))//'slice_y_grad_csh_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==5) then
				fn=dir(1:len_trim(dir))//'slice_y_grad_csv_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			endif

			open(unit=10,file=fn,action='write')

			do ind_z=0,nz
				do lz=0,lpd
					do ind_x=0,nx
						do lx=0,lpd
							write(10,*) slice(ind_x,ind_z,lx,lz)
						enddo
					enddo
				enddo
			enddo

			close(unit=10)

			deallocate(slice,stat=status)

		!==============================================================
		! z=const plane
		!==============================================================

		elseif ((plane==3) .and. (val<=zmax) .and. (val>zmin)) then

			write(93,*) junk

			!======================================================
			! allocate slice
			!======================================================

			allocate(slice(0:nx,0:ny,0:lpd,0:lpd),stat=status)

			!======================================================
			! z-coordinate line
			!======================================================

			dz=(zmax-zmin)/(nz+1)		! width of one element in z direction

			do k=0,nz
				do n=0,lpd

					z(k,n)=zmax-k*dz-0.5*(1+knots(n))*dz

					if (abs(val-z(k,n))<abs(val-z(index_n,index_l))) then

						index_n=k
						index_l=n

					endif
				enddo
			enddo

			write(*,*) index_n, index_l, z(index_n,index_l)

			!======================================================
			! open files and transfer values
			!======================================================

			if (comp==1) then
				fn=dir(1:len_trim(dir))//'grad_rho_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==2) then
				fn=dir(1:len_trim(dir))//'grad_cp_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==3) then
				fn=dir(1:len_trim(dir))//'grad_cs_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==4) then
				fn=dir(1:len_trim(dir))//'grad_csh_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==5) then
				fn=dir(1:len_trim(dir))//'grad_csv_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			endif

			write(*,*) 'opening ', fn
			open(unit=10,file=fn,action='read',form='unformatted')
			read(10) u(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd)
			close(unit=10)

			slice(0:nx,0:ny,0:lpd,0:lpd)=u(0:nx,0:ny,index_n,0:lpd,0:lpd,index_l)

			!==============================================================
			! write slices to files
			!==============================================================

			if (comp==1) then
				fn=dir(1:len_trim(dir))//'slice_z_grad_rho_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==2) then
				fn=dir(1:len_trim(dir))//'slice_z_grad_cp_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==3) then
				fn=dir(1:len_trim(dir))//'slice_z_grad_cs_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==4) then
				fn=dir(1:len_trim(dir))//'slice_z_grad_csh_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			elseif (comp==5) then
				fn=dir(1:len_trim(dir))//'slice_z_grad_csv_'//junk(1:len_trim(junk))//'_'//ns(1:len_trim(ns))
			endif

			open(unit=10,file=fn,action='write')

			do ind_y=0,ny
				do ly=0,lpd
					do ind_x=0,nx
						do lx=0,lpd
							write(10,*) slice(ind_x,ind_y,lx,ly)
						enddo
					enddo
				enddo
			enddo

			close(unit=10)

			deallocate(slice,stat=status)

		endif

	enddo

	close(unit=15)
	close(unit=93)

        !======================================================================
        ! clean up
        !======================================================================

        deallocate(big_slice,stat=status)

end subroutine ses3d_make_slice_grad


