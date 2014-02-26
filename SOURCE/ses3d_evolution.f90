!*****************************************************************************
! time evolution *************************************************************
!*****************************************************************************
! last modified: 17 March 2011 by Andreas Fichtner
!*****************************************************************************

subroutine ses3d_evolution
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! local variables
	!======================================================================

	integer :: i,j,k,q, in, jn, kn
	real :: dummy_1, dummy_2, dummy_3
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: L, M_dummy, Mdx, Mdy, Mdz
	real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: r_sin_theta
	real, dimension(0:nx_max,0:ny_max,0:nz_max) :: dummy_x, dummy_y, dummy_z

	!======================================================================
	! map global v to local versions for computations of local stresses
	!======================================================================

	do i=0,nx
	do in=0,lpd

	  index_x=i*lpd+in	! global index in x-direction

	do j=0,ny
	do jn=0,lpd

	  index_y=j*lpd+jn	! global index in y-direction

	do k=0,nz
	do kn=0,lpd

	  index_z=k*lpd+kn	! global index in z-direction

	  vx(i,j,k,in,jn,kn)=vx_global(index_x,index_y,index_z)
	  vy(i,j,k,in,jn,kn)=vy_global(index_x,index_y,index_z)
	  vz(i,j,k,in,jn,kn)=vz_global(index_x,index_y,index_z)

	enddo
	enddo
	enddo
	enddo
	enddo
	enddo

        !======================================================================
        ! initialisations
        !======================================================================

        exx=0.0; eyy=0.0; ezz=0.0; exy=0.0; exz=0.0; eyz=0.0; eyx=0.0; ezx=0.0; ezy=0.0

        !----------------------------------------------------------------------
        !- derivative proxies of velocity field
        !----------------------------------------------------------------------

       	do i=0,lpd
           do q=0,lpd

           	dummy_3=2.0*dl(q,i)/dz
            	dummy_2=2.0*dl(q,i)/dy
            	dummy_1=2.0*dl(q,i)/dx

		!----------------------------------------------------------------
	      	!- velocity field -----------------------------------------------
	      	!----------------------------------------------------------------

            	ezz(:,:,:,:,:,i)=ezz(:,:,:,:,:,i)+vz(:,:,:,:,:,q)*dummy_3		! (grad v)_(r r)
            	ezy(:,:,:,:,:,i)=ezy(:,:,:,:,:,i)+vy(:,:,:,:,:,q)*dummy_3		! (grad v)_(r phi)
            	ezx(:,:,:,:,:,i)=ezx(:,:,:,:,:,i)+vx(:,:,:,:,:,q)*dummy_3		! (grad v)_(r theta)

            	eyy(:,:,:,:,i,:)=eyy(:,:,:,:,i,:)+vy(:,:,:,:,q,:)*dummy_2		! (grad v)_(phi phi)
	      	eyx(:,:,:,:,i,:)=eyx(:,:,:,:,i,:)+vx(:,:,:,:,q,:)*dummy_2		! (grad v)_(phi theta)
	      	eyz(:,:,:,:,i,:)=eyz(:,:,:,:,i,:)+vz(:,:,:,:,q,:)*dummy_2		! (grad v)_(phi r)

            	exx(:,:,:,i,:,:)=exx(:,:,:,i,:,:)+vx(:,:,:,q,:,:)*dummy_1		! (grad v)_(theta theta)
	      	exy(:,:,:,i,:,:)=exy(:,:,:,i,:,:)+vy(:,:,:,q,:,:)*dummy_1		! (grad v)_(theta phi)
	      	exz(:,:,:,i,:,:)=exz(:,:,:,i,:,:)+vz(:,:,:,q,:,:)*dummy_1		! (grad v)_(theta r)

          enddo
        enddo

	!-------------------------------------------------------------------------------------------------
        !- complete velocity gradients in spherical coordinates ------------------------------------------
        !-------------------------------------------------------------------------------------------------

	r_sin_theta=r*sin_theta

	ezz=-ezz					! (grad v)_(r r)
	eyy=eyy/(r_sin_theta)+(vz+vx*cot_theta)/r	! (grad v)_(phi phi)
	exx=(vz+exx)/r					! (grad v)_(theta theta)
	ezx=-ezx					! (grad v)_(r theta)
	exz=(exz-vx)/r					! (grad v)_(theta r)
	ezy=-ezy					! (grad v)_(r phi)
	eyz=-vy/r+eyz/(r_sin_theta)			! (grad v)_(phi r)
	exy=exy/r					! (grad v)_(theta phi)
	eyx=-vy*cot_theta/r+eyx/(r_sin_theta)		! (grad v)_(phi theta)

	!-------------------------------------------------------------------------------------------------
        !- integrate velocity gradients to displacement gradients ----------------------------------------
        !-------------------------------------------------------------------------------------------------

	dzuz=dzuz+dt*ezz
	dyuy=dyuy+dt*eyy
	dxux=dxux+dt*exx
	dzux=dzux+dt*ezx
	dxuz=dxuz+dt*exz
	dzuy=dzuy+dt*ezy
	dyuz=dyuz+dt*eyz
	dxuy=dxuy+dt*exy
	dyux=dyux+dt*eyx

	!-------------------------------------------------------------------------------------------------
	!- build apml strain rates -----------------------------------------------------------------------
	!-------------------------------------------------------------------------------------------------

        exx=exx+(prof_z+prof_y)*dxux*ispml
        eyy=eyy+(prof_x+prof_z)*dyuy*ispml
        ezz=ezz+(prof_x+prof_y)*dzuz*ispml

        exy=exy+(prof_y+prof_z)*dxuy*ispml
        eyx=eyx+(prof_x+prof_z)*dyux*ispml

        exz=exz+(prof_y+prof_z)*dxuz*ispml
        ezx=ezx+(prof_x+prof_y)*dzux*ispml

        eyz=eyz+(prof_x+prof_z)*dyuz*ispml
        ezy=ezy+(prof_x+prof_y)*dzuy*ispml

	exy=(exy+eyx)/2.0
	exz=(exz+ezx)/2.0
	eyz=(eyz+ezy)/2.0

	!======================================================================
	! make strong form stress rates
	!======================================================================

        ! - diagonal stress rates ---------------------------------------------

        if (is_diss==1) then	!- dissipation on

        	Mdx=0.0
           	Mdy=0.0
           	Mdz=0.0

           	do k=1,nrdiss

              		Mdx=Mdx+Mxx(k,:,:,:,:,:,:)
              		Mdy=Mdy+Myy(k,:,:,:,:,:,:)
              		Mdz=Mdz+Mzz(k,:,:,:,:,:,:)

           	enddo

           	M_dummy=Mdx+Mdy+Mdz

           	L=(kappa-2.0*mu_tau/3.0)*(exx+eyy+ezz)-2.0*mu*tau*M_dummy/3.0
              	sxx=L+2.0*mu_tau*exx+2.0*mu*tau*Mdx+C*ezz+A*(eyy+exx)
              	syy=L+2.0*mu_tau*eyy+2.0*mu*tau*Mdy+C*ezz+A*(eyy+exx)
              	szz=L+2.0*mu_tau*ezz+2.0*mu*tau*Mdz+C*(eyy+exx)

        else			!- dissipation off

              	L=lambda*(exx+eyy+ezz)
              	sxx=L+2*mu*exx+C*ezz+A*(eyy+exx)
              	syy=L+2*mu*eyy+C*ezz+A*(eyy+exx)
              	szz=L+2*mu*ezz+C*(eyy+exx)

        endif

	!- off-diagonal stress rates -------------------------------------------------

        if (is_diss==1) then		!- dissipation on

           	Mdx=0.0
           	Mdy=0.0
           	Mdz=0.0

           	do k=1,nrdiss

              		Mdx=Mdx+Mxy(k,:,:,:,:,:,:)
              		Mdy=Mdy+Mxz(k,:,:,:,:,:,:)
              		Mdz=Mdz+Myz(k,:,:,:,:,:,:)

           	enddo

              	sxy=2.0*mu_tau*exy+2.0*mu*tau*Mdx
              	sxz=2.0*mu_tau*exz+2.0*mu*tau*Mdy+2.0*B*exz
              	syz=2.0*mu_tau*eyz+2.0*mu*tau*Mdz+2.0*B*eyz

        else				!- dissipation off

              	sxy=2.0*mu*exy
              	sxz=2.0*mu*exz+2.0*B*exz
              	syz=2.0*mu*eyz+2.0*B*eyz

        endif

        !======================================================================
        ! march pml stresses (integrate stress rates)
        !======================================================================

        sxx_pml=sxx_pml+dt*(sxx-ispml*prof_x*sxx_pml)
        syy_pml=syy_pml+dt*(syy-ispml*prof_y*syy_pml)
        szz_pml=szz_pml+dt*(szz-ispml*prof_z*szz_pml)

        sxy_pml=sxy_pml+dt*(sxy-ispml*prof_x*sxy_pml)
        syx_pml=syx_pml+dt*(sxy-ispml*prof_y*syx_pml)

        sxz_pml=sxz_pml+dt*(sxz-ispml*prof_x*sxz_pml)
        szx_pml=szx_pml+dt*(sxz-ispml*prof_z*szx_pml)

        syz_pml=syz_pml+dt*(syz-ispml*prof_y*syz_pml)
        szy_pml=szy_pml+dt*(syz-ispml*prof_z*szy_pml)

	!======================================================================
	! make weak form of stress divergence
	!======================================================================

	!- make moment tensor source if there is one

	if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

	    call make_moment_source

	endif

	!- compute stress divergence (with moment tensor added to the stress tensor)

	sx=0.0; sy=0.0; sz=0.0

	do i=0,lpd
	do j=0,lpd
	do k=0,lpd

        	dummy_1=(2/dx)*Jac*w(j)*w(k)
        	dummy_2=(2/dy)*Jac*w(i)*w(k)
        	dummy_3=-(2/dz)*Jac*w(i)*w(j)

		do q=0,lpd

        		dummy_x=w(q)*dl(i,q)*r(:,:,:,q,j,k)*sin_theta(:,:,:,q,j,k)*dummy_1
        		dummy_y=w(q)*dl(j,q)*r(:,:,:,i,q,k)*dummy_2
        		dummy_z=w(q)*dl(k,q)*r(:,:,:,i,j,q)*r(:,:,:,i,j,q)*sin_theta(:,:,:,i,j,q)*dummy_3

        		sx(:,:,:,i,j,k)=sx(:,:,:,i,j,k) &
                       		+(szx_pml(:,:,:,i,j,q)+src_zx(:,:,:,i,j,q))*dummy_z &
                       		+(syx_pml(:,:,:,i,q,k)+src_yx(:,:,:,i,q,k))*dummy_y &
                       		+(sxx_pml(:,:,:,q,j,k)+src_xx(:,:,:,q,j,k))*dummy_x

        		sy(:,:,:,i,j,k)=sy(:,:,:,i,j,k) &
                       		+(szy_pml(:,:,:,i,j,q)+src_zy(:,:,:,i,j,q))*dummy_z &
                       		+(syy_pml(:,:,:,i,q,k)+src_yy(:,:,:,i,q,k))*dummy_y &
                       		+(sxy_pml(:,:,:,q,j,k)+src_xy(:,:,:,q,j,k))*dummy_x

        		sz(:,:,:,i,j,k)=sz(:,:,:,i,j,k) &
                       		+(szz_pml(:,:,:,i,j,q)+src_zz(:,:,:,i,j,q))*dummy_z &
                       		+(syz_pml(:,:,:,i,q,k)+src_yz(:,:,:,i,q,k))*dummy_y &
                       		+(sxz_pml(:,:,:,q,j,k)+src_xz(:,:,:,q,j,k))*dummy_x

		enddo

	enddo
	enddo
	enddo

	!======================================================================
	! add external forces single force
	!======================================================================

	call add_single_force

	!======================================================================
	! map local force vectors to global force vectors
	!======================================================================

	sx_global(:,:,:)=0.0
	sy_global(:,:,:)=0.0
	sz_global(:,:,:)=0.0

	do i=0,nx
	do in=0,lpd

	  index_x=i*lpd+in	! global index in x-direction

	do j=0,ny
	do jn=0,lpd

	  index_y=j*lpd+jn	! global index in y-direction

	do k=0,nz
	do kn=0,lpd

	  index_z=k*lpd+kn	! global index in z-direction

	  sx_global(index_x,index_y,index_z)=sx_global(index_x,index_y,index_z)+sx(i,j,k,in,jn,kn)
	  sy_global(index_x,index_y,index_z)=sy_global(index_x,index_y,index_z)+sy(i,j,k,in,jn,kn)
	  sz_global(index_x,index_y,index_z)=sz_global(index_x,index_y,index_z)+sz(i,j,k,in,jn,kn)

	enddo
	enddo
	enddo
	enddo
	enddo
	enddo

	call communicate_global_acceleration

        !======================================================================
        ! make accelerations / divide weak stresses by mass matrix
        !======================================================================

	sx_global=sx_global/MM_global
	sy_global=sy_global/MM_global
	sz_global=sz_global/MM_global

	!======================================================================
	! time extrapolation of the displacement field
	!======================================================================

	vx_global=vx_global+dt*(sx_global-ispml*prof_global*vx_global)
	vy_global=vy_global+dt*(sy_global-ispml*prof_global*vy_global)
	vz_global=vz_global+dt*(sz_global-ispml*prof_global*vz_global)

        !======================================================================
        ! time extrapolation of memory variables
        !======================================================================

        if (is_diss==1) then

           do k=1,nrdiss

              L=dt*D_p(k)/tau_p(k)

              Mxx(k,:,:,:,:,:,:)=Mxx(k,:,:,:,:,:,:)-dt*Mxx(k,:,:,:,:,:,:)/tau_p(k) &
                                -L(:,:,:,:,:,:)*exx(:,:,:,:,:,:)

              Myy(k,:,:,:,:,:,:)=Myy(k,:,:,:,:,:,:)-dt*Myy(k,:,:,:,:,:,:)/tau_p(k) &
                                -L(:,:,:,:,:,:)*eyy(:,:,:,:,:,:)

              Mzz(k,:,:,:,:,:,:)=Mzz(k,:,:,:,:,:,:)-dt*Mzz(k,:,:,:,:,:,:)/tau_p(k) &
                                -L(:,:,:,:,:,:)*ezz(:,:,:,:,:,:)

              Mxy(k,:,:,:,:,:,:)=Mxy(k,:,:,:,:,:,:)-dt*Mxy(k,:,:,:,:,:,:)/tau_p(k) &
                                -L(:,:,:,:,:,:)*exy(:,:,:,:,:,:)

              Mxz(k,:,:,:,:,:,:)=Mxz(k,:,:,:,:,:,:)-dt*Mxz(k,:,:,:,:,:,:)/tau_p(k) &
                                -L(:,:,:,:,:,:)*exz(:,:,:,:,:,:)

              Myz(k,:,:,:,:,:,:)=Myz(k,:,:,:,:,:,:)-dt*Myz(k,:,:,:,:,:,:)/tau_p(k) &
                                -L(:,:,:,:,:,:)*eyz(:,:,:,:,:,:)

           enddo

        endif

	!=======================================================================
	! taper boundaries when pmls are turned off
	!=======================================================================

	vx_global=taper_global*vx_global; vy_global=taper_global*vy_global; vz_global=taper_global*vz_global

      	sxx_pml=sxx_pml*taper; syy_pml=syy_pml*taper; szz_pml=szz_pml*taper
       	sxy_pml=sxy_pml*taper; syx_pml=syx_pml*taper
       	szx_pml=szx_pml*taper; sxz_pml=sxz_pml*taper
       	syz_pml=syz_pml*taper; szy_pml=szy_pml*taper

end subroutine ses3d_evolution
