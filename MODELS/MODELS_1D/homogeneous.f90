subroutine homogeneous
use variables_gm
implicit none
include 'mpif.h'

	!=============================================================
	! homogeneous model
	!=============================================================

	real :: rho_homogen, mu_homogen, lambda_homogen
	real :: A_homogen, B_homogen, C_homogen

	rho_homogen=1.0
	mu_homogen=0.0
	lambda_homogen=0.0
	A_homogen=0.0
	B_homogen=0.0
	C_homogen=0.0

	if (my_rank==0) then

		write(*,*) 'Homogeneous model'

        endif

	rhoinv(:,:,:,:,:,:)=1.0/rho_homogen
	mu(:,:,:,:,:,:)=mu_homogen
	lambda(:,:,:,:,:,:)=lambda_homogen
	A(:,:,:,:,:,:)=A_homogen
	B(:,:,:,:,:,:)=B_homogen
	C(:,:,:,:,:,:)=C_homogen
		
end subroutine homogeneous
