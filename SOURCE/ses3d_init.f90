!==============================================================================
! subroutine for the initialisation of space matrices, mass matrix, ...
!==============================================================================
! last modified: 1 April 2010 by Andreas Fichtner
!==============================================================================

subroutine ses3d_init
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! local variables
	!======================================================================

	character(len=40) :: dummy, junk
	character(len=5) :: s_idx

	integer :: i, j, k, l, m, n, in, jn, kn, idx

    real :: delta_x, delta_y, delta_z, alpha
	real :: tol
	
	!=====================================================================

    ispml=1.0
    tol=pi/(180*100000)                !- receiver location tolerance
    
	!======================================================================
	
	call int2str(my_rank, dummy)
	
    	if (my_rank==0) then
    
	    write(*,*)'begin init '
    	    
	endif
	
	write(99,*)'begin init '
	write(99,*)'------------------------------------------------------------'
  
    !======================================================================
    ! initialise elastic parameters, memory variables and Frechet derivatives
    !======================================================================

    kappa=lambda+2*mu/3
	cs=sqrt(mu*rhoinv)
	cp=sqrt((lambda+2*mu)*rhoinv)
	
	rho=0.0
	where (rhoinv>0)
		rho=1/rhoinv
	end where

        if (is_diss==1) then

        	mu_tau=(1+sum_D_p*tau)*mu
           	Mxx=0.0; Myy=0.0; Mzz=0.0; Mxy=0.0; Mxz=0.0; Myz=0.0

        endif

	if (adjoint_flag==2) then

		grad_rho=0.0; grad_csh=0.0; grad_csv=0; grad_cp=0.0

	endif
        
	dxux=0.0; dyux=0.0; dzux=0.0; dxuy=0.0; dyuy=0.0; dzuy=0.0; dxuz=0.0; dyuz=0.0; dzuz=0.0

	vx_global=0.0;	vy_global=0.0; 	vz_global=0.0;

	sxx_pml=0.0;	syy_pml=0.0;	szz_pml=0.0;
	sxy_pml=0.0;	syx_pml=0.0;
	sxz_pml=0.0;	szx_pml=0.0;
	syz_pml=0.0;	szy_pml=0.0;

	Mxx=0.0;	Myy=0.0;	Mzz=0.0;
	Mxy=0.0;	Mxz=0.0;	Myz=0.0;

	src_xx=0.0; src_yy=0.0; src_zz=0.0; src_xy=0.0; src_yx=0.0
	src_xz=0.0; src_zx=0.0; src_yz=0.0; src_zy=0.0

	grad_rho=0.0; 	grad_cp=0.0; 	grad_csh=0.0; 	grad_csv=0.0

	!======================================================================
	! determine collocation points (knots) and integration weights
	!======================================================================
	
	if (lpd==2) then
		
		knots(0)=-1.0
		knots(1)=0.0
		knots(2)=1.0
		
		w(0)=1.0/3.0
		w(1)=4.0/3.0
		w(2)=1.0/3.0
		
	elseif (lpd==3) then
		
		knots(0)=-1.0
		knots(1)=-0.4472135954999579
		knots(2)=0.4472135954999579
		knots(3)=1.0
		
		w(0)=0.1666666666666
		w(1)=0.8333333333333
		w(2)=0.8333333333333
		w(3)=0.1666666666666
		
	elseif (lpd==4) then
		
		knots(0)=-1.0
		knots(1)=-0.6546536707079772
		knots(2)=0.0
		knots(3)=0.6546536707079772
		knots(4)=1.0
		
		w(0)=0.1000000000000
		w(1)=0.5444444444444
		w(2)=0.7111111111111
		w(3)=0.5444444444444
		w(4)=0.1000000000000
		
	elseif (lpd==5) then
		
		knots(0)=-1.0
		knots(1)=-0.7650553239294647
		knots(2)=-0.2852315164806451
		knots(3)=0.2852315164806451
		knots(4)=0.7650553239294647
		knots(5)=1.0
		
		w(0)=0.0666666666666667
		w(1)=0.3784749562978470
		w(2)=0.5548583770354862
		w(3)=0.5548583770354862
		w(4)=0.3784749562978470
		w(5)=0.0666666666666667
		
	elseif (lpd==6) then
		
		knots(0)=-1.0
		knots(1)=-0.8302238962785670
		knots(2)=-0.4688487934707142
		knots(3)=0.0
		knots(4)=0.4688487934707142
		knots(5)=0.8302238962785670
		knots(6)=1.0
		
		w(0)=0.0476190476190476
		w(1)=0.2768260473615659
		w(2)=0.4317453812098627
		w(3)=0.4876190476190476
		w(4)=0.4317453812098627
		w(5)=0.2768260473615659
		w(6)=0.0476190476190476
		
	elseif (lpd==7) then
		
		knots(0)=-1.0
		knots(1)=-0.8717401485096066
		knots(2)=-0.5917001814331423
		knots(3)=-0.2092992179024789
		knots(4)=0.2092992179024789
		knots(5)=0.5917001814331423
		knots(6)=0.8717401485096066
		knots(7)=1.0
		
		w(0)=0.0357142857142857
		w(1)=0.2107042271435061
		w(2)=0.3411226924835044
		w(3)=0.4124587946587038
		w(4)=0.4124587946587038
		w(5)=0.3411226924835044
		w(6)=0.2107042271435061
		w(7)=0.0357142857142857
		
	endif

	!======================================================================
	! initialisation of LAGRANGE polynomial derivatives
	!======================================================================
	
	do i=0,lpd
		do j=0,lpd

			call dlgll(lpd,i,j,dl(i,j))

		enddo
	enddo
	
	!======================================================================
	! initialisation of the space matrices
	!======================================================================
	
	dx=(xmax-xmin)/(nx+1)		! width of one element in x direction
	dy=(ymax-ymin)/(ny+1)		! width of one element in y direction
	dz=(zmax-zmin)/(nz+1)		! width of one element in z direction
	
	write(99,*) '- element sizes --------------------------------------------' 
	write(99,*) 'width of one element in theta direction in deg: ', dx*180/pi
	write(99,*) 'width of one element in phi direction in deg: ', dy*180/pi
	write(99,*) 'width of one element in z direction in m: ', dz

	do i=0,nx
		do n=0,lpd
			x(i,n)=xmin+i*dx+0.5*(1+knots(n))*dx
                        sin_theta(i,:,:,n,:,:)=sin(x(i,n))
                        cot_theta(i,:,:,n,:,:)=cos(x(i,n))/sin(x(i,n))
                        cos_theta(i,:,:,n,:,:)=cos(x(i,n))
			write(101,*) x(i,n)
		enddo
	enddo

	do j=0,ny
		do n=0,lpd
			y(j,n)=ymin+j*dy+0.5*(1+knots(n))*dy
			write(102,*) y(j,n)
		enddo
	enddo
	
	do k=0,nz
		do n=0,lpd
			z(k,n)=zmax-k*dz-0.5*(1+knots(n))*dz
                        r(:,:,k,:,:,n)=z(k,n)
			write(103,*) z(k,n)
		enddo
	enddo
	
	!======================================================================
	! initialisation of the JACOBIAN
	!======================================================================
	
	Jac=dx*dy*dz/8

        !======================================================================
        ! initialisation of the global and local mass matrices
        !======================================================================
	
	MM_global(:,:,:)=0.0

        do i=0,nx
        do j=0,ny
        do k=0,nz
        do l=0,lpd
        do m=0,lpd
        do n=0,lpd

		index_x=i*lpd+l	! global index in x-direction
		index_y=j*lpd+m	! global index in y-direction
		index_z=k*lpd+n	! global index in z-direction

		MM_global(index_x,index_y,index_z)=MM_global(index_x,index_y,index_z) &
			-w(l)*w(m)*w(n)*Jac*z(k,n)*z(k,n)*sin_theta(i,j,k,l,m,n)/rhoinv(i,j,k,l,m,n)

		MM(i,j,k,l,m,n)=-w(l)*w(m)*w(n)*Jac*z(k,n)*z(k,n)*sin_theta(i,j,k,l,m,n)/rhoinv(i,j,k,l,m,n)

        enddo
        enddo
        enddo
        enddo
        enddo
        enddo 

	call communicate_global_field(MM_global)

        !=======================================================================
        ! find receivers located in a particular processor box
        !=======================================================================

	if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

        	nr=0			! number of receivers in this processor box

        	write(99,*) 'receiver locations in transformed unit system --------------'


        	do n=1,nr_global	! find receivers that are located in this processor box

        	if ((recloc_global(1,n)<=xmax) .and. (recloc_global(1,n)>xmin) .and. &
               	    (recloc_global(2,n)<=ymax) .and. (recloc_global(2,n)>ymin) .and. &
               	    (zmax_global-recloc_global(3,n)<=zmax) .and. (zmax_global-recloc_global(3,n)>zmin)) then

               	write(99,*) 'receiver: ', & 
		station_name_global(n),recloc_global(1,n)*180/pi, &
		recloc_global(2,n)*180/pi,zmax_global-recloc_global(3,n)

               	nr=nr+1

               	do i=0,nx
               	do j=0,ny
               	do k=0,nz

               	!------------------------------------------------------------------------------------------------
               	!- REMARK: Due to round off errors it becomes necessary to introduce an absolute tolerance in 
               	!  order to locate a receiver within an element. This may lead to a doubling of a receiver that
               	!  is located in principle exactly on an element boundary.
               	!------------------------------------------------------------------------------------------------

               	if ((recloc_global(1,n)<=x(i,lpd)+tol) .and. (recloc_global(1,n)>=x(i,0)-tol) .and. &
                    (recloc_global(2,n)<=y(j,lpd)+tol) .and. (recloc_global(2,n)>=y(j,0)-tol) .and. &
                    (zmax_global-recloc_global(3,n)<=z(k,0)) .and. (zmax_global-recloc_global(3,n)>=z(k,lpd))) then

                           
                	! receiver coordinates in unit system

                	recloc_std(1,nr)=2*((recloc_global(1,n)-xmin)/dx-i)-1
                	recloc_std(2,nr)=2*((recloc_global(2,n)-ymin)/dy-j)-1
                	recloc_std(3,nr)=-2*((zmax_global-recloc_global(3,n)-zmax)/dz+k)-1

                	recloc(1,nr)=recloc_global(1,n)
                	recloc(2,nr)=recloc_global(2,n)
                	recloc(3,nr)=recloc_global(3,n)

			station_name_local(nr)=station_name_global(n)

                	rx(nr)=i
                	ry(nr)=j
                	rz(nr)=k
                           
			write(99,*) nz,k,z(k,:)
                	write(99,*) 'element: (', i,',',j,',',k,')'
			write(99,*) 'standard coordinates: (', recloc_std(1,nr),',',recloc_std(2,nr),',',recloc_std(3,nr),')'

               	endif	!- if receiver is in this element

                enddo	!- loop over elements
                enddo	!- loop over elements
               	enddo	!- loop over elements

            	endif	!- if the receiver is in this processor box

        	enddo	!- loop over all receivers

		write(99,*)'------------------------------------------------------------'
        	write(99,*) 'number of receivers: ', nr
		write(99,*)'------------------------------------------------------------'

        	seismogram_x(:,:)=0.0
        	seismogram_y(:,:)=0.0
        	seismogram_z(:,:)=0.0

	endif
	
     	!=========================================================================
	! point source location
	!=========================================================================
	
	if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

		if (is_source==1) then
		
			delta_x=dx
			delta_y=dy
			delta_z=dz
		
			! search for source location in x direction
		
			do i=0,nx
			do n=0,lpd
			
				if (abs(x(i,n)-xxs)<delta_x) then
					
					delta_x=abs(x(i,n)-xxs)
					isx=i
					isx_n=n
					
				endif
				
			enddo
			enddo
		
			! search for source location in y direction
		
			do j=0,ny
			do n=0,lpd
			
				if (abs(y(j,n)-yys)<delta_y) then
					
					delta_y=abs(y(j,n)-yys)
					isy=j
					isy_n=n
					
				endif
				
			enddo
			enddo
		
			! search for source location in z direction
		
			do k=0,nz
			do n=0,lpd
			
				if (abs(zmax_global-z(k,n)-zzs)<delta_z) then
					
					delta_z=abs(zmax_global-z(k,n)-zzs)
					isz=k
					isz_n=n
					
				endif
				
			enddo
			enddo

		endif	!- is_source==1

	endif
	
	!======================================================================
	! make local coordinates of the source point
	!======================================================================

	if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

		if (is_source==1) then

			zzs_loc=(z(isz,lpd)+z(isz,0)-2*(zmax_global-zzs))/(z(isz,0)-z(isz,lpd))
			yys_loc=(2*yys-y(isy,lpd)-y(isy,0))/(y(isy,lpd)-y(isy,0))		
			xxs_loc=(2*xxs-x(isx,lpd)-x(isx,0))/(x(isx,lpd)-x(isx,0))

			write(99,*) 'local source coordinates: ', xxs_loc, yys_loc, zzs_loc

		endif

	endif

	!=============================================================================
        ! make local adjoint source locations and read adjoint source time functions
        !=============================================================================

	if (adjoint_flag==2) then

	  call int2str(event_indices(i_events),junk)

	  nr_adsrc=0

	  write(99,*) '------------------------------------------------------------'
	  write(99,*) '- adjoint sources ******************************************'
	  write(99,*) '------------------------------------------------------------'

	  do k=1,nr_adsrc_global

	    !- find adjoint sources in this processor box ----------------------------

	    if ((ad_srcloc_global(1,k)<=xmax) .and. (ad_srcloc_global(1,k)>xmin) .and. &
                (ad_srcloc_global(2,k)<=ymax) .and. (ad_srcloc_global(2,k)>ymin) .and. &
                (zmax_global-ad_srcloc_global(3,k)<=zmax) .and. (zmax_global-ad_srcloc_global(3,k)>zmin)) then
              
		nr_adsrc=nr_adsrc+1

		ad_srcloc(1,nr_adsrc)=ad_srcloc_global(1,k)
		ad_srcloc(2,nr_adsrc)=ad_srcloc_global(2,k)
		ad_srcloc(3,nr_adsrc)=ad_srcloc_global(3,k)

		write(99,*) 'adjoint source ', k, ' (global), adjoint source ', nr_adsrc, ' (local)'
		write(99,*) 'position (lat,lon,depth): ', ad_srcloc(1,nr_adsrc), ad_srcloc(2,nr_adsrc), ad_srcloc(3,nr_adsrc)
		write(*,*) 'adjoint source in processor', my_rank

		call int2str(k,s_idx)

		!- open file containing the adjoint source time function -------------

		open(unit=10,file='../ADJOINT/'//junk(1:len_trim(junk))//'/ad_src_'//s_idx(1:len_trim(s_idx)),action='read')
              	
		!- read header of the adjoint source file (not used) -----------------
	
		read(10,*) dummy
		read(10,*) dummy
		read(10,*) dummy
		read(10,*) dummy

		!- read adjoint source time function components ----------------------

		do i=1,nt
                 
		  read(10,*) ad_stf_x(nr_adsrc,i), ad_stf_y(nr_adsrc,i), ad_stf_z(nr_adsrc,i)
                 		
              	enddo

		close(unit=10)

	    endif

	  enddo

	  write(99,*) '------------------------------------------------------------'

	endif

	!======================================================================
	! find indices of adjoint source locations and make local coordinates
	!======================================================================

	if (adjoint_flag==2) then

	  do idx=1,nr_adsrc

	    !- search in x-direction ------------------------------------------

	    delta_x=xmax-xmin

	    do i=0,nx
            do k=0,lpd

	      if (abs(ad_srcloc(1,idx)-x(i,k))<delta_x) then

		is(1,idx)=i
                isn(1,idx)=k
		delta_x=abs(ad_srcloc(1,idx)-x(i,k))

              endif

	    enddo
	    enddo

	    xxs_ad_loc(idx)=(2*ad_srcloc(1,idx)-x(is(1,idx),lpd)-x(is(1,idx),0))/(x(is(1,idx),lpd)-x(is(1,idx),0))

	    !- search in y-direction ------------------------------------------

	    delta_y=ymax-ymin

	    do i=0,ny
            do k=0,lpd

	      if (abs(ad_srcloc(2,idx)-y(i,k))<delta_y) then

		is(2,idx)=i
                isn(2,idx)=k
    		delta_y=abs(ad_srcloc(2,idx)-y(i,k))

	      endif

	    enddo
	    enddo

	    yys_ad_loc(idx)=(2*ad_srcloc(2,idx)-y(is(2,idx),lpd)-y(is(2,idx),0))/(y(is(2,idx),lpd)-y(is(2,idx),0))

	    !- search in z-direction ------------------------------------------

	    delta_z=zmax-zmin

	    do i=0,nz
            do k=0,lpd

	      if (abs(zmax_global-ad_srcloc(3,idx)-z(i,k))<delta_z) then

		is(3,idx)=i
                isn(3,idx)=k
		delta_z=abs(zmax_global-ad_srcloc(3,idx)-z(i,k))

	      endif

	    enddo
	    enddo

	    zzs_ad_loc(idx)=(z(is(3,idx),lpd)+z(is(3,idx),0)-2*(zmax_global-ad_srcloc(3,idx)))/(z(is(3,idx),0)-z(is(3,idx),lpd))

	  enddo

	endif

        !======================================================================
        ! PML damping profiles
        !======================================================================

        alpha=1.2

        prof_x=0.0
        prof_y=0.0
        prof_z=0.0
        
        !- upper z-boundary ---------------------------------------------------

	if (1==0) then

        if (iz_multi==pz) then

           do k=0,pml-1
              do n=0,lpd

                 delta_z=(zmax-pml*dz-z(k,n))/(pml*dz)
                 prof_z(0:nx,0:ny,k,0:lpd,0:lpd,n)=alpha*delta_z*delta_z

              enddo
           enddo

        endif

	endif

        !- lower z-boundary ---------------------------------------------------

        if (iz_multi==1) then

           do k=(nz-pml+1),nz
              do n=0,lpd

                 delta_z=(zmin+pml*dz-z(k,n))/(pml*dz)
                 prof_z(0:nx,0:ny,k,0:lpd,0:lpd,n)=alpha*delta_z*delta_z

              enddo
           enddo

        endif

        !- left x-boundary ----------------------------------------------------
                 
        if (ix_multi==1) then

           do i=0,pml-1
              do n=0,lpd

                 delta_x=(xmin+pml*dx-x(i,n))/(pml*dx)
                 prof_x(i,0:ny,0:nz,n,0:lpd,0:lpd)=alpha*delta_x*delta_x

              enddo
           enddo

        endif

        !- right x-boundary ---------------------------------------------------

        if (ix_multi==px) then

           do i=(nx-pml+1),nx
              do n=0,lpd

                 delta_x=(xmax-dx*pml-x(i,n))/(pml*dx)
                 prof_x(i,0:ny,0:nz,n,0:lpd,0:lpd)=alpha*delta_x*delta_x

              enddo
           enddo

        endif

        !- left y-boundary ----------------------------------------------------
                 
        if (iy_multi==1) then

           do i=0,pml-1
              do n=0,lpd

                 delta_y=(ymin+pml*dy-y(i,n))/(pml*dy)
                 prof_y(0:nx,i,0:nz,0:lpd,n,0:lpd)=alpha*delta_y*delta_y

              enddo
           enddo

        endif

        !- right y-boundary ---------------------------------------------------

        if (iy_multi==py) then

           do i=(ny-pml+1),ny
              do n=0,lpd

                 delta_y=(ymax-dy*pml-y(i,n))/(pml*dy)
                 prof_y(0:nx,i,0:nz,0:lpd,n,0:lpd)=alpha*delta_y*delta_y

              enddo
           enddo

        endif

        !- normalise the corners ----------------------------------------------

        prof=prof_x+prof_y+prof_z

        prof_z=2*prof_z/(1+prof/alpha)
        prof_y=2*prof_y/(1+prof/alpha)
        prof_x=2*prof_x/(1+prof/alpha)

        prof=prof_x+prof_y+prof_z
       
        taper=exp(-0.1*prof*prof)

	!- map local prof to global prof --------------------------------------

	do i=0,nx
	do j=0,ny
	do k=0,nz
	do in=0,lpd
	do jn=0,lpd
	do kn=0,lpd

		index_x=i*lpd+in	! global index in x-direction
		index_y=j*lpd+jn	! global index in y-direction
		index_z=k*lpd+kn	! global index in z-direction

		prof_global(index_x,index_y,index_z)=prof(i,j,k,in,jn,kn)
		taper_global(index_x,index_y,index_z)=taper(i,j,k,in,jn,kn)
		
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo

	!======================================================================
	! clean up
	!======================================================================
	
	if (my_rank==0) then

		write(*,*) 'end init'
	
	endif
	
	write(99,*)'------------------------------------------------------------'
	write(99,*) 'end init'
	write(99,*)'------------------------------------------------------------'
	write(99,*)'------------------------------------------------------------'
	write(99,*)'iterations: ------------------------------------------------'
	write(99,*)'------------------------------------------------------------'	


end subroutine ses3d_init
