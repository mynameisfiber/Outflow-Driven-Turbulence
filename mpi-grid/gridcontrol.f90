program gridcontrol
include "mpif.h"

integer :: maxx,maxy,maxz
integer, DIMENSION(3) :: dims = (/ 3,3,2 /)
integer :: rank, myid, cartrank=0
integer :: ierr
integer, DIMENSION(3) :: coords

call MPI_INIT(ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, myid, ierr)

call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims, &
		(/ .TRUE.,.TRUE.,.TRUE. /),.FALSE.,COMM_CART,ierr)

call MPI_CART_MAP(MPI_COMM_WORLD,3,dims, &
	     (/ .TRUE.,.TRUE.,.TRUE. /),cartrank,ierr)


call MPI_CART_COORDS(COMM_CART,cartrank,3, &
	        coords,ierr)

print*,"procrank=",rank,"cartrank=",cartrank
print*,"    (x,y,z) = ",coords

call MPI_FINALIZE(ierr)

return
end
