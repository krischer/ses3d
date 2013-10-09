!*****************************************************************************
! computation of Frechet kernels *********************************************
!*****************************************************************************
! last modified: 1 April 2010 by Andreas Fichtner
!*****************************************************************************

subroutine ses3d_grad
use parameters
use variables
implicit none
include 'mpif.h'

        !======================================================================
	! The sensitivity densities assume that all measurements are made on velocity fields !!!
	!======================================================================

        if (saving_vector(nt+1-it)>0) then

		!==============================================================
        	! read dynamic fields
        	!==============================================================

           	!- read velocity fields ---------------------------------------

           	call ses3d_restore(vx_fw,'vx',it,4,1021)
           	call ses3d_restore(vy_fw,'vy',it,4,1022)
           	call ses3d_restore(vz_fw,'vz',it,4,1023)

		!- read strain rates ------------------------------------------

		call ses3d_restore(exx_fw,'exx',it,4,2021)
		call ses3d_restore(eyy_fw,'eyy',it,4,2022)
		call ses3d_restore(ezz_fw,'ezz',it,4,2023)

		call ses3d_restore(exy_fw,'exy',it,4,2024)
		call ses3d_restore(exz_fw,'exz',it,4,2025)
		call ses3d_restore(eyz_fw,'eyz',it,4,2026)

		!==============================================================
        	! compute Frechet kernels
        	!==============================================================

		div_u_fw=exx_fw+eyy_fw+ezz_fw
		div_u_rw=dxux+dyuy+dzuz

		e_ddot_e=dxux*exx_fw+dyuy*eyy_fw+dzuz*ezz_fw+(dxuy+dyux)*exy_fw+(dxuz+dzux)*exz_fw+(dyuz+dzuy)*eyz_fw

                grad_cp=grad_cp+samp_ad*dt*(2*rho*cp*div_u_fw*div_u_rw)/Jac

                grad_rho=grad_rho+samp_ad*dt*(vx_fw*vx+vy_fw*vy+vz_fw*vz)/Jac
		grad_rho=grad_rho+samp_ad*dt*((cp*cp-2*cs*cs)*div_u_fw*div_u_rw+2*cs*cs*e_ddot_e)/Jac

		grad_csh=grad_csh+samp_ad*dt*(4*rho*cs*(e_ddot_e-div_u_rw*div_u_fw-exz_fw*(dxuz+dzux)-eyz_fw*(dyuz+dzuy)))/Jac
		grad_csv=grad_csv+samp_ad*dt*(4*rho*cs*(exz_fw*(dxuz+dzux)+eyz_fw*(dyuz+dzuy)))/Jac

        endif

end subroutine ses3d_grad
