module variables_gm
implicit none

	!- parameters and indices ---------------------------------------------

	character(len=40) :: junk, mrc
	integer :: i, j, k, l, m, n, s, done
	integer :: mi, is_aniso, is_diss
	integer, allocatable, dimension(:,:) :: min_index, max_index		! minimum and maximum index boundaries
	integer :: ierr, p, my_rank
	
	real :: knots(0:7)

        !- global model setup -------------------------------------------------

	integer :: model_type, status
	integer :: nx, ny, nz, px, py, pz, lpd
	integer :: bnx, bny, bnz

	real :: z_max, z_min, y_max, y_min, x_max, x_min, dz, dy, dx

        !- local model setup --------------------------------------------------

        integer :: nx_loc, ny_loc, nz_loc
        
        real :: x_max_loc, y_max_loc, z_max_loc
        real :: x_min_loc, y_min_loc, z_min_loc

        real, allocatable, dimension(:,:) :: x, y, z
	real, allocatable, dimension(:,:,:,:,:,:) :: XX, YY, ZZ

        !- physical model -----------------------------------------------------

	real, allocatable, dimension(:,:,:,:,:,:) :: rhoinv, mu, lambda, A, B, C, cp, cs, rho, Q
        real, allocatable, dimension(:,:) :: cp1d, cs1d, rho1d, r, Q1d
			
	real, dimension(1:139) :: ak135_r, ak135_cp, ak135_cs, ak135_rho
	real, dimension(1:151) :: ak135c_r, ak135c_cp, ak135c_cs, ak135c_rho
	real :: ratio
			
        real, parameter :: pi=3.14159265358979

	!- ocean index --------------------------------------------------------

	real, allocatable, dimension(:,:,:,:) :: oci

	!- andere -------------------------------------------------------------

	integer :: integer_dummy
	real :: real_dummy


end module variables_gm
