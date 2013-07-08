
!*****************************************************************************
!************** main program for spherical spectral elements *****************
!*****************************************************************************

program ses3d_main
use parameters
use variables
implicit none
include 'mpif.h'

        character(len=8) :: date1,date2
	character(len=10) :: time1,time2  
	character(len=20) :: dummy
	real :: vx_max, vx_min

	!======================================================================
	! initialisations
	!======================================================================
        
	call mpi_init(ierr)
	call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
	call mpi_comm_size(mpi_comm_world, p, ierr)
 
	!======================================================================
	! read master Par_master file and communicate content
	!======================================================================

	if (my_rank==0) then

	  open(unit=99,file='../INPUT/event_list',action='read')

	  read(99,*) n_events

	  do it=1,n_events
	    read(99,*) event_indices(it)
	  enddo

	  close(unit=99)

	endif

	call mpi_bcast(n_events,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(event_indices,1000,mpi_integer,0,mpi_comm_world,ierr)

	!======================================================================
	! open and write logfile
	!======================================================================
	
	call int2str(my_rank, dummy)
    
	open(unit=99, file='../DATA/LOGFILES/logfile'//dummy(1:len_trim(dummy)), action='WRITE')
    
	write(99,*) '---------------------------------------'
	write(99,*) 'logfile process ', my_rank
	write(99,*) '---------------------------------------'
    
	if (my_rank==0) then
	
		call date_and_time(date1, time1)
		
		write(*,*) '------------------------------------------------------------------'
		write(*,*) 'starting date: ', date1(7:8), '. ', date1(5:6), '. ', date1(5:5), date1(4:4)
		write(*,*) 'starting time: ', time1(1:2), ':', time1(3:4), ',', time1(5:6)
		write(*,*) '------------------------------------------------------------------'
		
	endif
    
	!======================================================================
	! start loop over events
	!======================================================================

	do i_events=1,n_events

	    !======================================================================
	    ! read parameters from files 'Par' and 'recfile' and 'boxfile
	    !======================================================================

	    call ses3d_input
	
	    !======================================================================
	    ! various initializations (model space, source, receivers, boundaries)
	    !======================================================================

	    call ses3d_init
  
	    !======================================================================
	    ! start time evolution
	    !======================================================================

	    do it=1,nt

		!==============================================================
		! forward time stepping
		!==============================================================
		
		write(99,*) "iteration", it, "ispml=", ispml
		
		call ses3d_evolution

		if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

			call record_seismograms

		endif
		
		call mpi_reduce(maxval(vx),vx_max,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
		call mpi_reduce(minval(vx),vx_min,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)

		!==============================================================
		! (screen) output
		!==============================================================

		if (my_rank==0) then

			if ((adjoint_flag==0) .or. (adjoint_flag==1)) then

			write(*,*) 'iteration ',it ,':   ', vx_min, '< vx <', vx_max, 'ispml=',ispml

			elseif (adjoint_flag==2) then

   			write(*,*) 'adjoint iteration ',it ,':   ', vx_min, '< vx <', vx_max, 'ispml=',ispml

			endif

		endif
		
                if (mod(it,ssamp)==0) then

                   call ses3d_output(it)
                   
                endif

		!==============================================================
		! save forward wavefield
		!==============================================================

		if (adjoint_flag==1) then

                	if ((mod(it,samp_ad)==0) .or. (it==1) .or. (it==nt)) then

				!- store strain rate --------------------------
		       
				call ses3d_store(vx,'vx',it,lpd,1021)
				call ses3d_store(vy,'vy',it,lpd,1022)
				call ses3d_store(vz,'vz',it,lpd,1023)

				call ses3d_store(dxux,'exx',it,lpd,1041)
				call ses3d_store(dyuy,'eyy',it,lpd,1042)
				call ses3d_store(dzuz,'ezz',it,lpd,1043)
		
				call ses3d_store((dxuy+dyux)/2,'exy',it,lpd,1044)
				call ses3d_store((dxuz+dzux)/2,'exz',it,lpd,1045)
				call ses3d_store((dyuz+dzuy)/2,'eyz',it,lpd,1046)

				!- document saving vector ---------------------
				
                        	saving_vector(it)=it

                	endif

		endif

		!==============================================================
		! compute Frechet derivatives
		!==============================================================

		if (adjoint_flag==2) then

			call ses3d_grad

		endif

	    enddo	! end of iteration, do it=1,nt
   
	    call ses3d_output(nt)

	enddo		! end of loop over events

	!======================================================================
	! date and time output
	!======================================================================

	if (my_rank==0) then
	
		call date_and_time(date2, time2)
		
		write(*,*) '------------------------------------------------------------------'
		write(*,*) 'starting date: ', date1(7:8), '. ', date1(5:6), '. ', date1(5:5), date1(4:4)
		write(*,*) 'starting time: ', time1(1:2), ':', time1(3:4), ',', time1(5:6)
		write(*,*) '------------------------------------------------------------------'
		write(*,*) 'finishing date: ', date2(7:8), '. ', date2(5:6), '. ', date2(5:5), date2(4:4)
		write(*,*) 'finishing time: ', time2(1:2), ':', time2(3:4), ',', time2(5:6)
		write(*,*) '------------------------------------------------------------------'
		
	endif

	close(99)
	call mpi_finalize(ierr)

end program ses3d_main
