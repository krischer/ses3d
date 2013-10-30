program ses3d_make_slice
use parameters
use variables
implicit none

	integer :: field_type

	write(*,*) 'field type (displacement=1, model=2, gradient=3): '
	read(*,*) field_type
	write(*,*) 'directory: '
	read(*,*) dir

	if (field_type==1) then

		call ses3d_make_slice_displacement

	elseif (field_type==2) then

		call ses3d_make_slice_model

	endif

end program ses3d_make_slice


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
