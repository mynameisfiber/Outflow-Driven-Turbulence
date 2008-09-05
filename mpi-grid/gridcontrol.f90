program gridcontrol
include "mpif.h"

integer :: n=126, localn=42
integer, DIMENSION(3) :: dims, offset
integer :: rank, myid, cartrank=0
integer :: ierr
integer, DIMENSION(3) :: coords

call MPI_INIT(ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, myid, ierr)

dims = (/ n/localn,n/localn,n/localn /)
call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims, &
		(/ .TRUE.,.TRUE.,.TRUE. /),.FALSE.,COMM_CART,ierr)
call MPI_CART_MAP(MPI_COMM_WORLD,3,dims, &
	     (/ .TRUE.,.TRUE.,.TRUE. /),cartrank,ierr)

call MPI_CART_COORDS(COMM_CART,cartrank,3, &
	        coords,ierr)
offset = coords * localn

print*,"procrank=",rank,"cartrank=",cartrank,"(x,y,z) = ",offset

call MPI_FINALIZE(ierr)

return
end
