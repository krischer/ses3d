!*****************************************************************************
!***************** parameters and global variables ***************************
!*****************************************************************************
! last modified: 1 April 2010 by Andreas Fichtner
!*****************************************************************************
!=============================================================================
! Predefined parameters
!=============================================================================

module parameters
implicit none

	integer, parameter :: nx_max=22		! number of elements in x direction per processor - 1
	integer, parameter :: ny_max=27		! number of elements in y direction per processor - 1
	integer, parameter :: nz_max=7		! number of elements in z direction per processor - 1
	
	integer, parameter :: lpd=4		! LAGRANGE polynomial degree
	
	integer, parameter :: maxnt=15000	! maximum number of time steps
	integer, parameter :: maxnr=800         ! maximum number of receivers
        integer, parameter :: pml=2
	
        integer, parameter :: nrdiss=3          ! number of relaxation mechanisms

	real, parameter :: pi=3.1415926535898

end module parameters
    
!==============================================================================
! global variables
!==============================================================================

module variables
use parameters
implicit none

        integer :: run_the_programme
        real :: ispml
        
	!======================================================================
	! model dimensions
	!======================================================================

	integer :: nx, ny, nz			! number of elements in the respective coordinate directions
	integer :: ix_multi, iy_multi, iz_multi	! multi index of the elements

	integer :: px, py, pz			! number of processors in the respective coordinate directions
	integer :: index_x, index_y, index_z
	integer :: p, my_rank, ierr

	integer :: n_events, i_events
	integer, dimension(1:1000) :: event_indices

	!======================================================================
	! physical model parameters
	!======================================================================
	
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: rhoinv, mu, lambda, kappa, mu_tau
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: A, B, C
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: tau, QQ
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: cp, cs, rho
	real, dimension(1:nrdiss) :: tau_p, D_p
	real :: sum_D_p
	
	!======================================================================
	! local displacement fields, strain fields, mass matrix and memory variables
	!======================================================================

	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: MM
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: ux, uy, uz
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: vx, vy, vz
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: ax, ay, az

        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: exx, eyy, ezz
        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: exy, exz, eyz
        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: eyx, ezx, ezy

        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: dxux, dyux, dzux
        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: dxuy, dyuy, dzuy
        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: dxuz, dyuz, dzuz

        real, dimension(1:nrdiss,0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: Mxx, Myy, Mzz
        real, dimension(1:nrdiss,0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: Mxy, Mxz, Myz

	!======================================================================
	! variables for communication
	!======================================================================

	real, dimension(0:2,0:(nx_max+1)*lpd,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: F
	real, dimension(0:2,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: Bsx
        real, dimension(0:2,0:(nx_max+1)*lpd,0:(nz_max+1)*lpd) :: Bsy
        real, dimension(0:2,0:(nx_max+1)*lpd,0:(ny_max+1)*lpd) :: Bsz

	!======================================================================
	! source fields
	!======================================================================
	
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: src_xx, src_yy, src_zz 
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: src_xy, src_yx, src_xz, src_zx, src_yz, src_zy

	!======================================================================
	! global displacement fields, strain fields, mass matrix and memory variables
	!======================================================================

	real, dimension(0:(nx_max+1)*lpd,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: MM_global
	real, dimension(0:(nx_max+1)*lpd,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: sx_global, sy_global, sz_global
	real, dimension(0:(nx_max+1)*lpd,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: ux_global, uy_global, uz_global
	real, dimension(0:(nx_max+1)*lpd,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: vx_global, vy_global, vz_global

	!======================================================================
	! weak form stress fields
	!======================================================================
	
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: sx, sy, sz
	
	!======================================================================
	! strong form stress fields
	!======================================================================
	
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: sxx,syy,szz
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: sxy,sxz,syz

        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: sxx_pml,syy_pml,szz_pml
        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: sxy_pml,syz_pml,sxz_pml
        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: syx_pml,szy_pml,szx_pml

        !======================================================================
        ! PML parameters
        !======================================================================

        real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: prof_x, prof_y, prof_z, prof, taper
	real, dimension(0:(nx_max+1)*lpd,0:(ny_max+1)*lpd,0:(nz_max+1)*lpd) :: prof_global, taper_global
        real :: E_kin, E_kin_old
	
	!======================================================================
	! geometrical parameters
	!======================================================================
	
	real, dimension(0:7) :: knots, w
	real, dimension(0:7,0:7) :: dl
	
	real :: dx, dy, dz
	real :: xmin,xmax,ymin,ymax,zmin,zmax
	real :: xmin_global, xmax_global, ymin_global, ymax_global, zmin_global, zmax_global
	
	real, dimension(0:nx_max,0:lpd) :: x
	real, dimension(0:ny_max,0:lpd) :: y
	real, dimension(0:nz_max,0:lpd) :: z
	
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: sin_theta, r, cot_theta, cos_theta

	real :: Jac
	
	!======================================================================
	! time parameters
	!======================================================================
	
	integer :: it, nt
	real :: dt

	!======================================================================
	! Source Variables
	!======================================================================

	integer :: source_type
	integer :: isx, isy, isz, isx_n, isy_n, isz_n
        
	real :: xxs, yys, zzs
	real :: xxs_loc, yys_loc, zzs_loc
	real :: is_source, source_processor
	real, dimension(1:maxnt) :: so
	
        real :: MOM_xx, MOM_yy, MOM_zz, MOM_xy, MOM_xz, MOM_yz
        
	!======================================================================
	! receiver variables
	!======================================================================
	
        integer :: nr_global, nr                   		                        	! number of receivers (first line of recfile)
	real :: recloc_global(1:3,1:maxnr), recloc(1:3,1:maxnr), recloc_std(1:3,1:maxnr)     	! receiver location x, y, depth
        integer :: rx(1:maxnr), ry(1:maxnr), rz(1:maxnr)

	character(len=12) :: station_name_global(1:maxnr), station_name_local(1:maxnr)

        !======================================================================
        ! output
        !======================================================================
        
        character(len=100) :: ofd                                                                       ! output file directory
        character(len=100) :: ffd                                                                       ! forward field directory
	character(len=100) :: dir

	integer :: output_displacement
	integer :: ssamp
        
        real :: seismogram_x(1:maxnr,1:maxnt)
        real :: seismogram_y(1:maxnr,1:maxnt)
        real :: seismogram_z(1:maxnr,1:maxnt)
	
	!======================================================================
	! other variables
	!======================================================================
	
	integer :: is_diss, integer_dummy

        !======================================================================
        ! adjoint computations
        !======================================================================

	integer :: adjoint_flag
        integer :: samp_ad
        integer :: saving_vector(1:maxnt)
	integer :: nr_adsrc_global, nr_adsrc

	integer, dimension(1:3,1:maxnr) :: is, isn

	real :: ad_srcloc_global(3,maxnr), ad_srcloc(3,maxnr)
	real, dimension(1:maxnr) :: xxs_ad_loc, yys_ad_loc, zzs_ad_loc
	real, dimension(1:maxnr,1:maxnt) :: ad_stf_x, ad_stf_y, ad_stf_z

	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: grad_rho, grad_cp, grad_csh, grad_csv
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: vx_fw, vy_fw, vz_fw
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: exx_fw, eyy_fw, ezz_fw
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: exy_fw, exz_fw, eyz_fw

	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: div_u_fw, div_u_rw, e_ddot_e

end module variables
!==============================================================================
