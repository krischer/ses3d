subroutine homogeneous_plus_Q
use variables_gm
implicit none
include 'mpif.h'

	real :: rho_homogen, mu_homogen, lambda_homogen
	real :: A_homogen, B_homogen, C_homogen

	allocate(r(0:nz_loc,0:lpd))
	allocate(Q1d(0:nz_loc,0:lpd))

	!=============================================================
	! homogeneous model (elastic part)
	!=============================================================

	if (my_rank==0) then

		write(*,*) 'Homogeneous model'

    endif
	
	rho_homogen=1.0
	mu_homogen=0.0
	lambda_homogen=0.0
	A_homogen=0.0
	B_homogen=0.0
	C_homogen=0.0

	rhoinv(:,:,:,:,:,:)=1.0/rho_homogen
	mu(:,:,:,:,:,:)=mu_homogen
	lambda(:,:,:,:,:,:)=lambda_homogen
	A(:,:,:,:,:,:)=A_homogen
	B(:,:,:,:,:,:)=B_homogen
	C(:,:,:,:,:,:)=C_homogen
	
	!=================================================================================
	!- make Q model (smooth version of QL6)
	!=================================================================================

	r=(6371.0-z/1000.0)/271.0

	!- make 1D profile

	if (is_diss==1) then

		do i=0,nz_loc
    	do k=0,lpd

			if ((z(i,k)<=6371000.0) .and. (z(i,k)>=6100000.0)) then
				  
				Q1d(i,k)=300.0-5370.82*r(i,k)**2+14401.62*r(i,k)**3-13365.78*r(i,k)**4+4199.98*r(i,k)**5
			     
			elseif ((z(i,k)<=6100000.0) .and. (z(i,k)>=5701000.0)) then
				  
				Q1d(i,k)=165.0

        	elseif ((z(i,k)<=5701000.0) .and. (z(i,k)>=3480000.0)) then
                             
            	Q1d(i,k)=355.0

    		elseif((z(i,k)<=3480000.0) .and. (z(i,k)>=1221000.0)) then

            	Q1d(i,k)=0.0

        	else

            	Q1d(i,k)=104.0

        	endif

		enddo
		enddo
		
	else
	
		Q1d(:,:)=100000.0
		
	endif
	
	!- translate to 3D array

	do i=0,nz_loc
    	do k=0,lpd

		  	Q(:,:,i,:,:,k)=Q1d(i,k)

    	enddo
	enddo
		
end subroutine homogeneous_plus_Q
