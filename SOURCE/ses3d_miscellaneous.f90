!*****************************************************************************
! sources and communication **************************************************
!*****************************************************************************
! last modified: 17 March 2011 by Andreas Fichtner
!*****************************************************************************

!==============================================================================
! add the forward field source
!==============================================================================

subroutine add_single_force
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! variables
	!======================================================================
	
	integer :: i, in, jn, kn, idx
	real :: delta
	real :: lgll

	!======================================================================
	!- single point source, forward calculation
	!======================================================================
    
	if ((adjoint_flag==0) .or. (adjoint_flag==1)) then
	
	  if ((is_source>0) .and. (source_type<10)) then

	    do in=0,lpd
	    do jn=0,lpd
	    do kn=0,lpd

	      !- representation of the delta function at the nodes of the source-bearing element

	      delta=lgll(lpd,in,xxs_loc)*lgll(lpd,jn,yys_loc)*lgll(lpd,kn,zzs_loc)
	
	      !- add single force to stress divergence

	      if (source_type==1) then
		sx(isx,isy,isz,in,jn,kn)=sx(isx,isy,isz,in,jn,kn)-delta*so(it)
			
	      elseif (source_type==2) then
		sy(isx,isy,isz,in,jn,kn)=sy(isx,isy,isz,in,jn,kn)-delta*so(it)

	      elseif (source_type==3) then
		sz(isx,isy,isz,in,jn,kn)=sz(isx,isy,isz,in,jn,kn)-delta*so(it)
	  
	      endif

	    enddo
	    enddo
	    enddo

	  endif

	endif

	!======================================================================
	!- adjoint point sources, adjoint calculation
	!======================================================================

	if (adjoint_flag==2) then

	  !- loop through all the adjoint sources in this processor box -------

	  do idx=1,nr_adsrc

	    do in=0,lpd
	    do jn=0,lpd
	    do kn=0,lpd

	      !- representation of the delta function at the nodes of the source-bearing element

	      delta=lgll(lpd,in,xxs_ad_loc(idx))*lgll(lpd,jn,yys_ad_loc(idx))*lgll(lpd,kn,zzs_ad_loc(idx))

	      !- add single force to stress divergence

	      sx(is(1,idx),is(2,idx),is(3,idx),in,jn,kn)=sx(is(1,idx),is(2,idx),is(3,idx),in,jn,kn)-delta*ad_stf_x(idx,it)
	      sy(is(1,idx),is(2,idx),is(3,idx),in,jn,kn)=sy(is(1,idx),is(2,idx),is(3,idx),in,jn,kn)-delta*ad_stf_y(idx,it)
	      sz(is(1,idx),is(2,idx),is(3,idx),in,jn,kn)=sz(is(1,idx),is(2,idx),is(3,idx),in,jn,kn)-delta*ad_stf_z(idx,it)

	    enddo
	    enddo
	    enddo

	  enddo
  
	endif

end subroutine add_single_force


subroutine make_moment_source
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! variables
	!======================================================================
	
	integer :: i, in, jn, kn
	real :: delta
	real :: lgll

	!======================================================================
	!- make moment tensor source fields
	!======================================================================
	
	src_xx=0.0; src_yy=0.0; src_zz=0.0; src_xy=0.0; src_yx=0.0
	src_xz=0.0; src_zx=0.0; src_yz=0.0; src_zy=0.0

	if ((is_source>0) .and. (source_type==10)) then
	
		do in=0,lpd
		do jn=0,lpd
		do kn=0,lpd

			!- representation of the delta function at the nodes of the source-bearing element

			delta=lgll(lpd,in,xxs_loc)*lgll(lpd,jn,yys_loc)*lgll(lpd,kn,zzs_loc)/(w(in)*w(jn)*w(kn))
			delta=delta/(z(isz,kn)*z(isz,kn)*sin(x(isx,in))*Jac)

			!- make source field

			src_xx(isx,isy,isz,in,jn,kn)=src_xx(isx,isy,isz,in,jn,kn)-MOM_xx*delta*so(it)
			src_yy(isx,isy,isz,in,jn,kn)=src_yy(isx,isy,isz,in,jn,kn)-MOM_yy*delta*so(it)
			src_zz(isx,isy,isz,in,jn,kn)=src_zz(isx,isy,isz,in,jn,kn)-MOM_zz*delta*so(it)

			src_xy(isx,isy,isz,in,jn,kn)=src_xy(isx,isy,isz,in,jn,kn)-MOM_xy*delta*so(it)
			src_yx(isx,isy,isz,in,jn,kn)=src_yx(isx,isy,isz,in,jn,kn)-MOM_xy*delta*so(it)
			
			src_xz(isx,isy,isz,in,jn,kn)=src_xz(isx,isy,isz,in,jn,kn)-MOM_xz*delta*so(it)
			src_zx(isx,isy,isz,in,jn,kn)=src_zx(isx,isy,isz,in,jn,kn)-MOM_xz*delta*so(it)

			src_yz(isx,isy,isz,in,jn,kn)=src_yz(isx,isy,isz,in,jn,kn)-MOM_yz*delta*so(it)
			src_zy(isx,isy,isz,in,jn,kn)=src_zy(isx,isy,isz,in,jn,kn)-MOM_yz*delta*so(it)

		enddo
		enddo
		enddo

	endif

end subroutine make_moment_source

!==============================================================================
! communicate global field variables across processors
!==============================================================================

subroutine communicate_global_field(FF)
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! variables
	!======================================================================

	real, dimension(0:(nx+1)*lpd,0:(ny+1)*lpd,0:(nz+1)*lpd), intent(inout) :: FF

	integer :: mi_rec, mi_src
	integer ::num
        integer :: status(MPI_STATUS_SIZE)

        real, dimension(0:(ny+1)*lpd,0:(nz+1)*lpd) :: Bx
        real, dimension(0:(nx+1)*lpd,0:(nz+1)*lpd) :: By
        real, dimension(0:(nx+1)*lpd,0:(ny+1)*lpd) :: Bz

	!======================================================================
	! communicate in x direction (send to the right, receive from the left)
	!======================================================================

	num=((ny+1)*lpd+1)*((nz+1)*lpd+1)	! number of nodes

	mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi	! index (rank) of receiving processor
        mi_src=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-2	! index (rank) of sending processor

	!- send to the right and receive from the left

	if (ix_multi<px) then					! send if there is a processor to the right
		call mpi_send(FF((nx+1)*lpd,:,:),num,mpi_real,mi_rec,10,mpi_comm_world,ierr)
	endif

	if (ix_multi>1) then 					! receiver if there is a processor to the left
		call mpi_recv(Bx(:,:),num,mpi_real,mi_src,10,mpi_comm_world,status,ierr)
		FF(0,:,:)=FF(0,:,:)+Bx(:,:)
		call mpi_send(FF(0,:,:),num,mpi_real,mi_src,11,mpi_comm_world,ierr)
	endif

	if (ix_multi<px) then
		call mpi_recv(FF((nx+1)*lpd,:,:),num,mpi_real,mi_rec,11,mpi_comm_world,status,ierr)
	endif

	!======================================================================
	! communicate in y direction (send to the right, receive from the left)
	!======================================================================

	num=((nx+1)*lpd+1)*((nz+1)*lpd+1)

	mi_rec=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi-1
        mi_src=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi-1

	if (iy_multi<py) then
		call mpi_send(FF(:,(ny+1)*lpd,:),num,mpi_real,mi_rec,20,mpi_comm_world,ierr)
        endif

        if (iy_multi>1) then
		call mpi_recv(By(:,:),num,mpi_real,mi_src,20,mpi_comm_world,status,ierr)
		FF(:,0,:)=FF(:,0,:)+By(:,:)
		call mpi_send(FF(:,0,:),num,mpi_real,mi_src,21,mpi_comm_world,ierr)
	endif

        if (iy_multi<py) then
		call mpi_recv(FF(:,(ny+1)*lpd,:),num,mpi_real,mi_rec,21,mpi_comm_world,status,ierr)
        endif

	!======================================================================
	! communicate in z direction (send to the right, receive from the left)
	!======================================================================
	! remember that the elements are upside down in z-direction !!!
	!======================================================================	

	num=((nx+1)*lpd+1)*((ny+1)*lpd+1)
 
	mi_rec=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi-1
        mi_src=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi-1

	!- send to bottom and receive from the top

	if (iz_multi<pz) then
		call mpi_send(FF(:,:,0),num,mpi_real,mi_rec,30,mpi_comm_world,ierr)
        endif

        if (iz_multi>1) then
		call mpi_recv(Bz(:,:),num,mpi_real,mi_src,30,mpi_comm_world,status,ierr)
		FF(:,:,(nz+1)*lpd)=FF(:,:,(nz+1)*lpd)+Bz(:,:)
		call mpi_send(FF(:,:,(nz+1)*lpd),num,mpi_real,mi_src,31,mpi_comm_world,ierr)
	endif

       if (iz_multi<pz) then
		call mpi_recv(FF(:,:,0),num,mpi_real,mi_rec,31,mpi_comm_world,status,ierr)
        endif

end subroutine communicate_global_field

!==============================================================================
! communicate global vector field variables across processors
!==============================================================================

subroutine communicate_global_acceleration
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! variables
	!======================================================================

!	real, dimension(0:2,0:(nx+1)*lpd,0:(ny+1)*lpd,0:(nz+1)*lpd) :: F

	integer :: mi_rec, mi_src
	integer ::num
        integer :: status(MPI_STATUS_SIZE)

!        real, dimension(0:2,0:(ny+1)*lpd,0:(nz+1)*lpd) :: Bx
!        real, dimension(0:2,0:(nx+1)*lpd,0:(nz+1)*lpd) :: By
!        real, dimension(0:2,0:(nx+1)*lpd,0:(ny+1)*lpd) :: Bz

	!======================================================================
	! group vector into one variable
	!======================================================================

	F(0,:,:,:)=sx_global(:,:,:)
	F(1,:,:,:)=sy_global(:,:,:)
	F(2,:,:,:)=sz_global(:,:,:)

	!======================================================================
	! communicate in x direction (send to the right, receive from the left)
	!======================================================================

	num=3*((ny_max+1)*lpd+1)*((nz_max+1)*lpd+1)	! number of nodes

	mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi	! index (rank) of receiving processor
        mi_src=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-2	! index (rank) of sending processor

	!- send to the right and receive from the left

	if (ix_multi<px) then					! send if there is a processor to the right
		call mpi_send(F(:,(nx+1)*lpd,:,:),num,mpi_real,mi_rec,10,mpi_comm_world,ierr)
	endif

	if (ix_multi>1) then 					! receiver if there is a processor to the left
		call mpi_recv(Bsx(:,:,:),num,mpi_real,mi_src,10,mpi_comm_world,status,ierr)
		F(:,0,:,:)=F(:,0,:,:)+Bsx(:,:,:)
		call mpi_send(F(:,0,:,:),num,mpi_real,mi_src,11,mpi_comm_world,ierr)
	endif

	if (ix_multi<px) then
		call mpi_recv(F(:,(nx+1)*lpd,:,:),num,mpi_real,mi_rec,11,mpi_comm_world,status,ierr)
	endif

	!======================================================================
	! communicate in y direction (send to the right, receive from the left)
	!======================================================================

	num=3*((nx_max+1)*lpd+1)*((nz_max+1)*lpd+1)

	mi_rec=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi-1
        mi_src=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi-1

	if (iy_multi<py) then
		call mpi_send(F(:,:,(ny+1)*lpd,:),num,mpi_real,mi_rec,20,mpi_comm_world,ierr)
        endif

        if (iy_multi>1) then
		call mpi_recv(Bsy(:,:,:),num,mpi_real,mi_src,20,mpi_comm_world,status,ierr)
		F(:,:,0,:)=F(:,:,0,:)+Bsy(:,:,:)
		call mpi_send(F(:,:,0,:),num,mpi_real,mi_src,21,mpi_comm_world,ierr)
	endif

        if (iy_multi<py) then
		call mpi_recv(F(:,:,(ny+1)*lpd,:),num,mpi_real,mi_rec,21,mpi_comm_world,status,ierr)
        endif

	!======================================================================
	! communicate in z direction (send to the right, receive from the left)
	!======================================================================
	! remember that the elements are upside down in z-direction !!!
	!======================================================================	

	num=3*((nx_max+1)*lpd+1)*((ny_max+1)*lpd+1)
 
	mi_rec=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi-1
        mi_src=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi-1

	!- send to bottom and receive from the top

	if (iz_multi<pz) then
		call mpi_send(F(:,:,:,0),num,mpi_real,mi_rec,30,mpi_comm_world,ierr)
        endif

        if (iz_multi>1) then
		call mpi_recv(Bsz(:,:,:),num,mpi_real,mi_src,30,mpi_comm_world,status,ierr)
		F(:,:,:,(nz+1)*lpd)=F(:,:,:,(nz+1)*lpd)+Bsz(:,:,:)
		call mpi_send(F(:,:,:,(nz+1)*lpd),num,mpi_real,mi_src,31,mpi_comm_world,ierr)
	endif

	if (iz_multi<pz) then
		call mpi_recv(F(:,:,:,0),num,mpi_real,mi_rec,31,mpi_comm_world,status,ierr)
	endif

	!======================================================================
	! restore vector components
	!======================================================================

	sx_global(:,:,:)=F(0,:,:,:)
	sy_global(:,:,:)=F(1,:,:,:)
	sz_global(:,:,:)=F(2,:,:,:)

end subroutine communicate_global_acceleration

!==============================================================================
! LAGRANGE-Polynomials for GAUSS-LOBATTO collocation points
!
! n=degree of the LAGRANGE polynomials (greater than 1)
! i=index of the LAGRANGE polynomial (between 0 and n)
!==============================================================================

real function lgll(n,i,x)
implicit none

	!======================================================================
	! variables
	!======================================================================
	
	integer, intent(in) :: n, i
	real, intent(in) :: x
	
	integer :: k
	real, dimension(0:7) :: knots
	
	lgll=1.0
	
	!======================================================================
	! determine collocation points (knots)
	!======================================================================
	
	if (n==2) then
		
		knots(0)=-1.0
		knots(1)=0.0
		knots(2)=1.0
		
	elseif (n==3) then
		
		knots(0)=-1.0
		knots(1)=-0.4472135954999579
		knots(2)=0.4472135954999579
		knots(3)=1.0
		
	elseif (n==4) then
		
		knots(0)=-1.0
		knots(1)=-0.6546536707079772
		knots(2)=0.0
		knots(3)=0.6546536707079772
		knots(4)=1.0
		
	elseif (n==5) then
		
		knots(0)=-1.0
		knots(1)=-0.7650553239294647
		knots(2)=-0.2852315164806451
		knots(3)=0.2852315164806451
		knots(4)=0.7650553239294647
		knots(5)=1.0
		
	elseif (n==6) then
		
		knots(0)=-1.0
		knots(1)=-0.8302238962785670
		knots(2)=-0.4688487934707142
		knots(3)=0.0
		knots(4)=0.4688487934707142
		knots(5)=0.8302238962785670
		knots(6)=1.0
		
	elseif (n==7) then
		
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
	! compute value of the LAGRANGE polynomial
	!======================================================================
	
	do k=0,n
	
		if (k /= i) then
	
			lgll=lgll*(x-knots(k))/(knots(i)-knots(k))
		
		endif
		
	enddo
	

end function lgll


!==============================================================================
! derivatives of LAGRANGE-GAUSS-LOBATTO polynomials at the collocation points
! n=degree of the LAGRANGE polynomial
! i=index of the LAGRANGE polynomial
! j=index of the collocation point where the derivative is evaluated
! y=output=derivative of the LAGRANGE polynomial at the collocation point j
!==============================================================================

subroutine dlgll(n,i,j,y)
implicit none

	!======================================================================
	! local variables
	!======================================================================
	
	integer, intent(in) :: n, i, j
	real, intent(out) :: y
	
	integer :: k
	real, dimension(0:7) :: knots
	
	!======================================================================
	! determine collocation points (knots)
	!======================================================================
	
	if (n==2) then
		
		knots(0)=-1.0
		knots(1)=0.0
		knots(2)=1.0
		
	elseif (n==3) then
		
		knots(0)=-1.0
		knots(1)=-0.4472135954999579
		knots(2)=0.4472135954999579
		knots(3)=1.0
		
	elseif (n==4) then
		
		knots(0)=-1.0
		knots(1)=-0.6546536707079772
		knots(2)=0.0
		knots(3)=0.6546536707079772
		knots(4)=1.0
		
	elseif (n==5) then
		
		knots(0)=-1.0
		knots(1)=-0.7650553239294647
		knots(2)=-0.2852315164806451
		knots(3)=0.2852315164806451
		knots(4)=0.7650553239294647
		knots(5)=1.0
		
	elseif (n==6) then
		
		knots(0)=-1.0
		knots(1)=-0.8302238962785670
		knots(2)=-0.4688487934707142
		knots(3)=0.0
		knots(4)=0.4688487934707142
		knots(5)=0.8302238962785670
		knots(6)=1.0
		
	elseif (n==7) then
		
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
	! compute derivative of LAGRANGE polynomial i at collocation point j
	!======================================================================
	
	if (i==j) then
		
		y=0.0
	
		do k=0,n
		
			if (k /= i) then
				
				y=y+1.0/(knots(i)-knots(k))
			
			endif
			
		enddo
		
	else
		
		y=1.0/(knots(i)-knots(j))
	
		do k=0,n
		
			if ((k /= i) .and. (k /= j)) then
				
				y=y*(knots(j)-knots(k))/(knots(i)-knots(k))
			
			endif
			
		enddo
		
	endif

end subroutine dlgll


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

