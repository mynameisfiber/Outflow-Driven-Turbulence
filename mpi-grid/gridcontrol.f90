program gridcontrol
include "mpif.h"


! Variable Definitions:
!     procs - number of processes
!     globaln - global length
!     n - length of local node
!     ghost - number of ghost cells
!     dims - # of nodes for each dimension
!     u - test grid
!     ghost - number of ghost zones
!     offset - coordinate in node space
!     coords - origin of local grid in the global grid
!     node - node ID in grid
!     ierr - MPI error code

integer, parameter :: procs = 8, n=16, globaln=(procs)**(1.0/3)*n, ghost=3
integer, parameter, DIMENSION(3) :: dims = (/ 2,2,2 /) 
real, DIMENSION(n,n,n,4) :: u
integer, DIMENSION(3) :: offset, coords
integer :: node
integer :: ierr

!Check that all variables make sense
IF (procs**(1.0/3) .ne. int(procs**(1.0/3))) THEN
  PRINT*,"# Procs must be a perfect cube"
END IF

!Initialize MPI
call MPI_INIT(ierr)

!Create the cartesian grid and find the local ID in it
call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,(/ .TRUE.,.TRUE.,.TRUE. /), &
     .FALSE.,COMM_CART,ierr)
call MPI_CART_MAP(MPI_COMM_WORLD,3,dims,(/ .TRUE.,.TRUE.,.TRUE. /), &
     node,ierr)

!Collects the local coordinates in node space and calculates the offset
call MPI_CART_COORDS(COMM_CART,node,3,offset,ierr)
coords = offset * n

print*,"node=",node,"(x,y,z) = ",coords

!Initialize the grid
u = node

!Now test out the boundry conditions
call boundary(u,n,ghost,node)

if (node .eq. 0) call output(u,n)


call MPI_FINALIZE(ierr)
return


CONTAINS

  SUBROUTINE boundary(u,n,ghost,node)
    integer :: n, node, ghost, dim, ndown, nup, count, ierr, tmp
    integer, dimension(3) :: sds, sde, sus, sue, rd, ru
    real, DIMENSION(n,n,n,4) :: u
    
    rd = (/ ghost, n, n /)
    ru = n-rd+1
    sds = (/ ghost+1, 1, 1 /)
    sde = (/ 2*ghost, n, n/)
    sus = n - sde + 1
    sue = n - sds + 1
    
    count = n*n*ghost*4
    call MPI_Barrier (COMM_CART,ierr)
    if (node .eq. 0)then
      print*,"rd",rd
      print*,"ru",ru
      print*,"sds",sds
      print*,"sde",sde
      print*,"sus",sus
      print*,"sue",sue
      print*,"size(ru)",size(u(ru(1):n,ru(2):n,ru(3):n,:))
      print*,"size(rd)",size(u(1:rd(1),1:rd(2),1:rd(3),:))
      print*,"size(sd)",size(u(sds(1):sde(1),sds(2):sde(2),sds(3):sde(3),:))
      print*,"size(su)",size(u(sus(1):sue(1),sus(2):sue(2),sus(3):sue(3),:))
      print*,"count",count
    end if
    DO dim=0,2
      !Find neighboors
      call MPI_Cart_shift (COMM_CART, dim, 1, &
           ndown, tmp, ierr)
      call MPI_Cart_shift (COMM_CART, dim, -1, &
           nup, tmp, ierr)
                          
      !First we send to the down
      print*,"D",node,"->",ndown
      call MPI_Sendrecv(u(sds(1):sde(1),sds(2):sde(2),sds(3):sde(3),:), count,&
                        MPI_REAL, ndown, node, u(ru(1):n,ru(2):n,ru(3):n,:), &
                        count, MPI_REAL, nup, nup, COMM_CART, &
                        MPI_STATUS_IGNORE, ierr)
      
      !Now we send to the up
      print*,"U",node,"->",nup
      call MPI_Sendrecv(u(sus(1):sue(1),sus(2):sue(2),sus(3):sue(3),:), count,&
                        MPI_REAL, nup, node, u(1:rd(1),1:rd(2),1:rd(3),:), &
                        count, MPI_REAL, ndown,ndown,COMM_CART, &
                        MPI_STATUS_IGNORE,ierr)
                        
      rd = cshift(rd,1)
      ru = cshift(ru,1)
      sds = cshift(sds,1)
      sde = cshift(sde,1)
      sus = cshift(sus,1)
      sue = cshift(sue,1)
    END DO
  END SUBROUTINE boundary
  
  SUBROUTINE output(u,n)
    integer :: n
    real, dimension(n,n,n,4) :: u
    
    OPEN(UNIT=1, FILE="output-test")
    DO i = 1,n
      DO j = 1,n
        DO k = 1, n
          write(1,900) u(i,j,k,1) , u(i,j,k,2),  u(i,j,k,3), u(i,j,k,4); 
          900 format(E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
               !i,j,k, rho, rho vx, rho vy, rho vz. 
        END DO !k
      END DO !j
    END DO !i
    CLOSE(1)
  END SUBROUTINE output
end
