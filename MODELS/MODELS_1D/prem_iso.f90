!================================================================================
!= Isotropic version of PREM. Dziewonski & Anderson, PEPI 1981.
!= Last modified on 5 February 2014 by Andreas Fichtner.
!================================================================================

subroutine prem_iso
use variables_gm
implicit none
include 'mpif.h'
		
	!=============================================================
	! EUROPEAN 1D BACKGROUND MODEL
	!=============================================================
		
	if (my_rank==0) then

		write(*,*) 'EUMOD 1D background'

    	endif

    	!- allocate memory

    	allocate(cp1d(0:nz_loc,0:lpd))
    	allocate(cs1d(0:nz_loc,0:lpd))
    	allocate(rho1d(0:nz_loc,0:lpd))
    	allocate(r(0:nz_loc,0:lpd))
	allocate(Q1d(0:nz_loc,0:lpd))

	!=================================================================================
	!- make elastic model ------------------------------------------------------------
	!=================================================================================

	r=z/6371e3

    	do i=0,nz_loc
    		do k=0,lpd

			!- crust
		
			if ((z(i,k)<=6371000.0) .and. (z(i,k)>=6356000.0)) then	! 0-15km
				  
				rho1d(i,k)=2.60
				cp1d(i,k)=5.80
				cs1d(i,k)=3.20
			     
			elseif ((z(i,k)<=6356000.0) .and. (z(i,k)>=6346600.0)) then	! 15-24.4km
				  
				rho1d(i,k)=2.90
				cp1d(i,k)=6.80
				cs1d(i,k)=3.90
		       
            		!- LID

        		elseif ((z(i,k)<=6346600.0) .and. (z(i,k)>=6291000.0)) then	! 24-80 km
                             
        			rho1d(i,k)=2.6910+0.6924*r(i,k)
            			cp1d(i,k)=4.1875+3.9382*r(i,k)
               			cs1d(i,k)=2.1519+2.3481*r(i,k)

        		!- LVZ
        	
        		elseif ((z(i,k)<=6291000.0) .and. (z(i,k)>=6151000.0)) then	! 80-220 km
        	
        			rho1d(i,k)=2.6910+0.6924*r(i,k)
            			cp1d(i,k)=4.1875+3.9382*r(i,k)
               			cs1d(i,k)=2.1519+2.3481*r(i,k)

            		!- Transition zone

        		elseif((z(i,k)<=6151000.0) .and. (z(i,k)>=5971000.0)) then	! 220-400 km

				rho1d(i,k)=7.1089-3.8045*r(i,k)
            			cp1d(i,k)=20.3926-12.2569*r(i,k)
            			cs1d(i,k)=8.9496-4.4597*r(i,k)

			elseif((z(i,k)<=5971000.0) .and. (z(i,k)>=5771000.0)) then	! 400-600 km

            			rho1d(i,k)=11.2494-8.0298*r(i,k)
            			cp1d(i,k)=39.7027-32.6166*r(i,k)
            			cs1d(i,k)=22.3512-18.5856*r(i,k)

        		elseif((z(i,k)<=5771000.0) .and. (z(i,k)>=5701000.0)) then	! 600-670 km

            			rho1d(i,k)=5.3197-1.4836*r(i,k)
               			cp1d(i,k)=19.0957-9.8672*r(i,k)
               			cs1d(i,k)=9.9839-4.9324*r(i,k)

        		!- Lower mantle

        		elseif((z(i,k)<=5701000.0) .and. (z(i,k)>=5600000.0)) then	! 670-771 km

        			rho1d(i,k)=7.9565-6.4761*r(i,k)+5.5283*r(i,k)*r(i,k)-3.0807*r(i,k)*r(i,k)*r(i,k)
            			cp1d(i,k)=29.2766-23.6026*r(i,k)+5.5242*r(i,k)*r(i,k)-2.5514*r(i,k)*r(i,k)*r(i,k)
               			cs1d(i,k)=22.3459-17.2473*r(i,k)-2.0834*r(i,k)*r(i,k)+0.9783*r(i,k)*r(i,k)*r(i,k)

          		elseif((z(i,k)<=5600000.0) .and. (z(i,k)>=3630000.0)) then	! 771-2741 km

               			rho1d(i,k)=7.9565-6.4761*r(i,k)+5.5283*r(i,k)*r(i,k)-3.0807*r(i,k)*r(i,k)*r(i,k)
               			cp1d(i,k)=24.9520-40.4673*r(i,k)+51.4832*r(i,k)*r(i,k)-26.6419*r(i,k)*r(i,k)*r(i,k)
               			cs1d(i,k)=11.1671-13.7818*r(i,k)+17.4575*r(i,k)*r(i,k)-9.2777*r(i,k)*r(i,k)*r(i,k)

    			elseif((z(i,k)<=3630000.0) .and. (z(i,k)>=3480000.0)) then	! 2741-2756 km

             			rho1d(i,k)=7.9565-6.4761*r(i,k)+5.5283*r(i,k)*r(i,k)-3.0807*r(i,k)*r(i,k)*r(i,k)
               			cp1d(i,k)=15.3891-5.3181*r(i,k)+5.5242*r(i,k)*r(i,k)-2.5514*r(i,k)*r(i,k)*r(i,k)
               			cs1d(i,k)=6.9254+1.4672*r(i,k)-2.0834*r(i,k)*r(i,k)+0.9783*r(i,k)*r(i,k)*r(i,k)

        		!- Outer core

        		elseif((z(i,k)<=3480000.0) .and. (z(i,k)>=1221500.0)) then

        			rho1d(i,k)=12.5815-1.2638*r(i,k)-3.6426*r(i,k)*r(i,k)-5.5281*r(i,k)*r(i,k)*r(i,k)
        			cp1d(i,k)=11.0487-4.0362*r(i,k)+4.8023*r(i,k)*r(i,k)-13.5732*r(i,k)*r(i,k)*r(i,k)
            			cs1d(i,k)=0.0

        		!- Inner core

        		elseif(z(i,k)<=1221500.0) then

            			rho1d(i,k)=13.0885-8.8381*r(i,k)*r(i,k)
            			cp1d(i,k)=11.2622-6.3640*r(i,k)*r(i,k)
        			cs1d(i,k)=3.6678-4.4475*r(i,k)*r(i,k)

        		endif

		enddo

    	enddo

	!=================================================================================
	!- make Q model 
	!=================================================================================

	r=(6371.0-z/1000.0)/271.0

	if (is_diss==1) then

		do i=0,nz_loc
    		do k=0,lpd

				!- crust
		
				if ((z(i,k)<=6371000.0) .and. (z(i,k)>=6356000.0)) then	! 0-15km
				  
					Q1d(i,k)=600.0	
			     
				elseif ((z(i,k)<=6356000.0) .and. (z(i,k)>=6346600.0)) then	! 15-24.4km
				  
					Q1d(i,k)=600.0	
		       
            			!- LID

        			elseif ((z(i,k)<=6346600.0) .and. (z(i,k)>=6291000.0)) then	! 24-80 km
                             
        				Q1d(i,k)=600.0	

        			!- LVZ
        	
        			elseif ((z(i,k)<=6291000.0) .and. (z(i,k)>=6151000.0)) then	! 80-220 km
        	
        				Q1d(i,k)=80.0	

            			!- Transition zone

        			elseif((z(i,k)<=6151000.0) .and. (z(i,k)>=5971000.0)) then	! 220-400 km

					Q1d(i,k)=143.0	

				elseif((z(i,k)<=5971000.0) .and. (z(i,k)>=5771000.0)) then	! 400-600 km

            				Q1d(i,k)=143.0	

        			elseif((z(i,k)<=5771000.0) .and. (z(i,k)>=5701000.0)) then	! 600-670 km

            				Q1d(i,k)=143.0		

        			!- Lower mantle

        			elseif((z(i,k)<=5701000.0) .and. (z(i,k)>=5600000.0)) then	! 670-771 km

        				Q1d(i,k)=312.0	

          			elseif((z(i,k)<=5600000.0) .and. (z(i,k)>=3630000.0)) then	! 771-2741 km

               				Q1d(i,k)=312.0	

    				elseif((z(i,k)<=3630000.0) .and. (z(i,k)>=3480000.0)) then	! 2741-2756 km

             				Q1d(i,k)=312.0	

        			!- Outer core

        			elseif((z(i,k)<=3480000.0) .and. (z(i,k)>=1221500.0)) then

        				Q1d(i,k)=100000.0	
        			
        			!- Inner core

        			elseif(z(i,k)<=1221500.0) then

            				Q1d(i,k)=84.6	

        			endif
			
			enddo
		enddo
		
	else
	
		Q1d(:,:)=100000.0
		
	endif

	!============================================================================
	! translate to 3D arrays
	!============================================================================

	rho1d=1000.0*rho1d
	cp1d=1000.0*cp1d
	cs1d=1000.0*cs1d

	do i=0,nz_loc
    	do k=0,lpd

    		rhoinv(:,:,i,:,:,k)=1/rho1d(i,k)
        	mu(:,:,i,:,:,k)=cs1d(i,k)*cs1d(i,k)*rho1d(i,k)
        	lambda(:,:,i,:,:,k)=(cp1d(i,k)*cp1d(i,k)-2*cs1d(i,k)*cs1d(i,k))*rho1d(i,k)
		Q(:,:,i,:,:,k)=Q1d(i,k)

    	enddo
	enddo

	!- initialise anisotropic parameters

	A(:,:,:,:,:,:)=0.0
	B(:,:,:,:,:,:)=0.0
	C(:,:,:,:,:,:)=0.0

	!- clean up

	deallocate(cp1d)
	deallocate(cs1d)
	deallocate(rho1d)
	deallocate(Q1d)
	deallocate(r)
                  
end subroutine prem_iso
